#include "WindVTKInput.h"
#include "WindVTKGrid.h"
#include "WINDSInputData.h"
#include "util/QESFileSystemHandler.h"
#include "util/QEStime.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

void WindVTKInput::parseValues()
{
  const std::vector<std::string> allowed={"file","terrainFile","velocityUnits","timeStamp","epsg",
    "originX","originY","originZ","siteZ0","reciprocal","terrainTolerance"};
  for (const auto &entry:tree) {
    if (entry.first=="<xmlcomment>") continue;
    if (std::find(allowed.begin(),allowed.end(),entry.first)==allowed.end() || tree.count(entry.first)!=1)
      throw ParseException("windVTK: unknown or repeated option "+entry.first);
  }
  parsePrimitive(true,file,"file");
  parsePrimitive(true,terrainFile,"terrainFile");
  parsePrimitive(true,velocityUnits,"velocityUnits");
  parsePrimitive(true,timeStamp,"timeStamp");
  parsePrimitive(true,epsg,"epsg");
  parsePrimitive(true,originX,"originX");
  parsePrimitive(true,originY,"originY");
  parsePrimitive(true,originZ,"originZ");
  parsePrimitive(false,siteZ0,"siteZ0");
  parsePrimitive(false,reciprocal,"reciprocal");
  parsePrimitive(false,terrainTolerance,"terrainTolerance");
  if (file.empty() || terrainFile.empty() || velocityUnits!="m/s" || epsg<=0
      || !std::isfinite(originX) || !std::isfinite(originY) || !std::isfinite(originZ)
      || !std::isfinite(siteZ0) || siteZ0<=0 || !std::isfinite(reciprocal)
      || !std::isfinite(terrainTolerance) || terrainTolerance<0)
    throw ParseException("windVTK: files, projected EPSG, m/s units, finite origins/stability and positive roughness required");
  QEStime checked(timeStamp);
  file=QESfs::get_absolute_path(file);
  terrainFile=QESfs::get_absolute_path(terrainFile);
}

void WindVTKInput::prepare(WINDSInputData *wid)
{
  auto &s=*wid->simParams;
  const auto &m=*wid->metParams;
  if (m.windShapefile || !m.sensors.empty() || !m.sensorName.empty() || s.wrfInputData || wid->hrrrInput)
    throw ParseException("windVTK cannot be combined with other wind input modes");
  if (!s.demFile.empty() || s.DTE_heightField)
    throw ParseException("windVTK: use terrainFile for the reference DEM; terrain is built from the VTK bottom layer");
  if (s.totalTimeIncrements>1 || s.meshTypeFlag!=0 || s.verticalStretching!=0
      || s.domainRotation!=0 || s.halo_x!=0 || s.halo_y!=0 || m.z0_domain_flag!=0
      || s.readCoefficientsFlag!=0 || s.logLawFlag!=0 || s.wrfCoupling
      || s.DEMDistancex!=0 || s.DEMDistancey!=0)
    throw ParseException("windVTK requires one timestep, stair terrain, uniform grid/roughness, zero rotation/halo/DEM offsets, fresh coefficients, no WRF coupling and logLawFlag=0");
  for (int c=0;c<3;++c)
    if (s.domain[c]<2 || !std::isfinite(s.grid[c]) || s.grid[c]<=0)
      throw ParseException("windVTK requires at least two cells per axis and positive finite spacing");
  if (s.originFlag!=0 && s.originFlag!=1)
    throw ParseException("windVTK requires originFlag 0 or 1; originX/Y/Z define the physical origin");
  // UTM fields are legacy float fields used by the building importer. Reject
  // conflicting user coordinates, and require exact float-representable origins
  // when that importer must subtract them from a projected building shapefile.
  if ((s.UTMx!=0 && s.UTMx!=float(originX)) || (s.UTMy!=0 && s.UTMy!=float(originY)))
    throw ParseException("windVTK: simulation UTMx/UTMy disagree with originX/Y");
  if (s.UTMZone!=0) {
    const int zone=(epsg>=32601 && epsg<=32660) ? epsg-32600 : (epsg>=32701 && epsg<=32760) ? epsg-32700 : 0;
    if (s.UTMZone!=zone) throw ParseException("windVTK: UTMZone disagrees with EPSG");
  }
  field=std::make_shared<WindVTKGrid>(file);
  const double xmax=originX+s.domain[0]*double(s.grid[0]);
  const double ymax=originY+s.domain[1]*double(s.grid[1]);
  field->ground(originX,originY); field->ground(xmax,ymax); // strict horizontal coverage
  field->validateTerrain(terrainFile,epsg,terrainTolerance,originX,originY,xmax,ymax);
  if (wid->buildingsParams && !wid->buildingsParams->shpFile.empty()) {
    field->validateVectorCRS(wid->buildingsParams->shpFile,epsg);
    if (originX==0 || originY==0 || double(float(originX))!=originX || double(float(originY))!=originY)
      throw ParseException("windVTK: building shapefile requires nonzero, exactly float-representable originX/Y (legacy building coordinates)");
  }
  s.UTMx=originX; s.UTMy=originY;
  // Check extrema on all source break lines clipped to the QES footprint.
  // Bilinear surfaces attain their extrema at these rectangle vertices.
  std::vector<double> cx={originX,xmax},cy={originY,ymax};
  for (double x:field->x()) if (x>originX && x<xmax) cx.push_back(x);
  for (double y:field->y()) if (y>originY && y<ymax) cy.push_back(y);
  const double ceiling=originZ+s.domain[2]*double(s.grid[2]);
  for (double y:cy) for (double x:cx) {
    if (field->ground(x,y)<originZ || field->ground(x,y)>=ceiling)
      throw ParseException("windVTK: vertical domain must start at/below terrain and end above all terrain");
    if (ceiling>field->top(x,y))
      throw ParseException("windVTK: missing vertical coverage at QES domain top");
  }
  const int nx=s.domain[0]+1,ny=s.domain[1]+1;
  std::vector<double> heights(size_t(nx)*ny);
  for (int j=0;j<ny;++j) for (int i=0;i<nx;++i)
    heights[i+size_t(nx)*j]=localGround(i*double(s.grid[0]),j*double(s.grid[1]));
  s.DTE_heightField=new DTEHeightField(heights,{nx,ny,1},{s.grid[0],s.grid[1],s.grid[2]},0,0);
  s.DTE_mesh=new Mesh(s.DTE_heightField->getTris());
  // Do not call setDomain(): the explicit originZ must not be re-zeroed.
}
double WindVTKInput::localGround(double x,double y) const
{
  return field->ground(originX+x,originY+y)-originZ;
}
