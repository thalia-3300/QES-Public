#include "WindProfilerShapefile.h"
#include "WINDSGeneralData.h"
#include <algorithm>
#include <iostream>
#include <limits>
#include <stdexcept>

WindProfilerShapefile::WindProfilerShapefile(const WINDSInputData *wid, const qes::Domain &domain)
  : grid(*wid->metParams->windShapefile)
{
  const auto &input = *wid->metParams->windShapefile;
  const auto &sim = *wid->simParams;
  grid.validateCoverage(input.originX, input.originY,
                        input.originX + (domain.nx()-1)*double(domain.dx()),
                        input.originY + (domain.ny()-1)*double(domain.dy()));
  if (!sim.demFile.empty())
    grid.validateDEM(sim.demFile, input.originX, input.originY, sim.originFlag,
                     sim.DEMDistancex, sim.DEMDistancey, sim.UTMx, sim.UTMy);
  if (wid->buildingsParams && !wid->buildingsParams->shpFile.empty())
    grid.validateVectorCRS(wid->buildingsParams->shpFile);
  std::cout << "[QES-WINDS]\t Wind shapefile: " << grid.nx() << " x " << grid.ny()
            << " points; reference height " << input.height << " m AGL; one timestep.\n";
}

void WindProfilerShapefile::interpolateWindProfile(const WINDSInputData *wid, WINDSGeneralData *wgd)
{
  const auto &input = *wid->metParams->windShapefile;
  const auto &d = wgd->domain;
  const int nx = d.nx(), ny = d.ny(), nz = d.nz();
  std::fill(wgd->u0.begin(), wgd->u0.end(), 0);
  std::fill(wgd->v0.begin(), wgd->v0.end(), 0);
  std::fill(wgd->w0.begin(), wgd->w0.end(), 0);
  const auto ground = [&](int i, int j) {
    return wgd->terrain[d.cell2d(std::clamp(i, 0, nx-2), std::clamp(j, 0, ny-2))];
  };
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      // Sample at staggered u/v faces, not cell centres. The extra tangential
      // padding row/column copies the last physical centre (not extrapolation).
      const double ux = input.originX + i*double(d.dx());
      const double uy = input.originY + (std::min(j, ny-2)+0.5)*d.dy();
      const double vx = input.originX + (std::min(i, nx-2)+0.5)*d.dx();
      const double vy = input.originY + j*double(d.dy());
      const double u = grid.sample(ux, uy)[0], v = grid.sample(vx, vy)[1];
      const double zu = 0.5*(ground(i-1,j) + ground(i,j));
      const double zv = 0.5*(ground(i,j-1) + ground(i,j));
      for (int k = 1; k < nz-1; ++k) {
        const auto face = d.face(i,j,k);
        const double up = u * input.profileScale(d.z[k] - zu);
        const double vp = v * input.profileScale(d.z[k] - zv);
        if (!std::isfinite(up) || !std::isfinite(vp)
            || std::abs(up) > std::numeric_limits<float>::max()
            || std::abs(vp) > std::numeric_limits<float>::max())
          throw std::runtime_error("windShapefile: vertical profile exceeds finite QES velocity range");
        wgd->u0[face] = up;
        wgd->v0[face] = vp;
      }
      wgd->u0[d.face(i,j,nz-1)] = wgd->u0[d.face(i,j,nz-2)];
      wgd->v0[d.face(i,j,nz-1)] = wgd->v0[d.face(i,j,nz-2)];
    }
  }
}
