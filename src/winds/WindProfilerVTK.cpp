#include "WindProfilerVTK.h"
#include "WindVTKGrid.h"
#include "WINDSGeneralData.h"
#include <algorithm>
#include <cstdio>
#include <stdexcept>

namespace {
bool solid(const WINDSGeneralData &w,int i,int j,int k)
{
  const auto &d=w.domain;
  if (i<0 || j<0 || k<0 || i>=d.nx()-1 || j>=d.ny()-1 || k>=d.nz()-1) return false;
  const int flag=w.icellflag[d.cell(i,j,k)];
  return flag==0 || flag==2;
}
void padding(const WINDSGeneralData &w,std::vector<float> &a,int component)
{
  const auto &d=w.domain;
  for (int k=0;k<d.nz();++k) for (int j=0;j<d.ny();++j) for (int i=0;i<d.nx();++i) {
    if (k==0) { a[d.face(i,j,k)]=0; continue; } // below-ground ghost layer
    const int ii=component==0 ? i : std::min(i,d.nx()-2);
    const int jj=component==1 ? j : std::min(j,d.ny()-2);
    const int kk=component==2 ? k : std::min(k,d.nz()-2);
    if (i!=ii || j!=jj || k!=kk) a[d.face(i,j,k)]=a[d.face(ii,jj,kk)];
  }
}
}
void WindProfilerVTK::interpolateWindProfile(const WINDSInputData *wid,WINDSGeneralData *w)
{
  const auto &input=*wid->metParams->windVTK;
  const auto &grid=*input.field;
  const auto &d=w->domain;
  std::vector<float> *arrays[]={&w->u0,&w->v0,&w->w0};
  std::vector<int> *masks[]={&w->vtk_u_mask,&w->vtk_v_mask,&w->vtk_w_mask};
  for (int c=0;c<3;++c) {
    auto &a=*arrays[c]; auto &mask=*masks[c];
    std::fill(a.begin(),a.end(),0); std::fill(mask.begin(),mask.end(),2);
    size_t sampled=0,blocked=0; double maxCast=0,minValue=0,maxValue=0;
    for (int k=1;k<d.nz()-(c==2?0:1);++k)
      for (int j=0;j<d.ny()-(c==1?0:1);++j)
        for (int i=0;i<d.nx()-(c==0?0:1);++i) {
          const auto id=d.face(i,j,k);
          const double x=input.originX+(i+(c==0?0:.5))*d.dx();
          const double y=input.originY+(j+(c==1?0:.5))*d.dy();
          const double z=input.originZ+(c==2?d.z_face[k]:d.z[k]);
          // A physical face touching a QES stair/building solid has zero flow.
          // Also mask faces below the continuous VTK surface (stair aliasing).
          const bool blockedFace=solid(*w,i,j,k) || solid(*w,i-(c==0),j-(c==1),k-(c==2))
            || z<=grid.ground(x,y);
          if (blockedFace) { mask[id]=0; ++blocked; continue; }
          const double value=grid.sample(x,y,z)[c];
          a[id]=value; mask[id]=1; ++sampled;
          maxCast=std::max(maxCast,std::abs(value-double(a[id])));
          minValue=std::min(minValue,value); maxValue=std::max(maxValue,value);
        }
    padding(*w,a,c);
    std::printf("[windVTK] %c: sampled=%zu, solid/terrain=%zu, range=[%.9g, %.9g] m/s; max float rounding=%.9g m/s\n",
                "uvw"[c],sampled,blocked,minValue,maxValue,maxCast);
    if (!sampled) throw std::runtime_error("windVTK: no covered air faces for a velocity component");
  }
  w->u0_imported=w->u0; w->v0_imported=w->v0; w->w0_imported=w->w0;
}
void WindProfilerVTK::enforceMasks(WINDSGeneralData *w)
{
  std::vector<float> *a[]={&w->u0,&w->v0,&w->w0};
  const std::vector<int> *m[]={&w->vtk_u_mask,&w->vtk_v_mask,&w->vtk_w_mask};
  for (int c=0;c<3;++c) {
    for (size_t id=0;id<a[c]->size();++id) if ((*m[c])[id]==0) (*a[c])[id]=0;
    padding(*w,*a[c],c);
  }
}
void WindProfilerVTK::finishVelocity(WINDSGeneralData *w)
{
  // Solid faces were constrained by the same geometry used for solver coefficients.
  // Only update padding here; changing physical faces after the solve would
  // invalidate its mass-conservation equation.
  padding(*w,w->u,0); padding(*w,w->v,1); padding(*w,w->w,2);
}
void WindProfilerVTK::reportChange(const WINDSGeneralData *w)
{
  const std::vector<float> *a[]={&w->u0,&w->v0,&w->w0};
  const std::vector<float> *b[]={&w->u0_imported,&w->v0_imported,&w->w0_imported};
  const std::vector<int> *m[]={&w->vtk_u_mask,&w->vtk_v_mask,&w->vtk_w_mask};
  for (int c=0;c<3;++c) {
    double maximum=0,squares=0; size_t n=0,changed=0;
    for (size_t id=0;id<a[c]->size();++id) if ((*m[c])[id]!=2) {
      const double delta=double((*a[c])[id])-(*b[c])[id];
      if (!std::isfinite(delta)) throw std::runtime_error("windVTK: non-finite parameterised velocity");
      maximum=std::max(maximum,std::abs(delta)); squares+=delta*delta; ++n; changed+=delta!=0;
    }
    std::printf("[windVTK] Parameterisations %c: changed=%zu/%zu physical faces, max delta=%.9g, RMS delta=%.9g m/s\n",
                "uvw"[c],changed,n,maximum,n?std::sqrt(squares/n):0);
  }
}
void WindProfilerVTK::reportFinal(const WINDSGeneralData *w)
{
  const auto &d=w->domain;
  double maxAll=0,sumAll=0,maxBoundary=0; size_t n=0;
  for (int k=1;k<d.nz()-1;++k) for (int j=0;j<d.ny()-1;++j) for (int i=0;i<d.nx()-1;++i) {
    if (solid(*w,i,j,k)) continue;
    const auto f=d.face(i,j,k);
    const double div=(double(w->u[f+1])-w->u[f])/d.dx()
      +(double(w->v[f+d.nx()])-w->v[f])/d.dy()
      +(double(w->w[f+d.nx()*d.ny()])-w->w[f])/d.dz_array[k];
    if (!std::isfinite(div)) throw std::runtime_error("windVTK: non-finite final divergence");
    maxAll=std::max(maxAll,std::abs(div)); sumAll+=div*div; ++n;
    if (i==0 || j==0 || i==d.nx()-2 || j==d.ny()-2 || k==d.nz()-2)
      maxBoundary=std::max(maxBoundary,std::abs(div));
  }
  std::printf("[windVTK] Final all-air divergence: cells=%zu, max=%.9g, RMS=%.9g s^-1; prescribed boundary/top max=%.9g s^-1\n",
              n,maxAll,n?std::sqrt(sumAll/n):0,maxBoundary);
}
