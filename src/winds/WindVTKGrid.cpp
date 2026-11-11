#include "WindVTKGrid.h"
#include <gdal_priv.h>
#include <ogrsf_frmts.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>

namespace {
[[noreturn]] void fail(const std::string &s) { throw std::runtime_error("windVTK: " + s); }
using Dataset = std::unique_ptr<GDALDataset, decltype(&GDALClose)>;
std::string line(std::ifstream &f)
{
  std::string s;
  while (std::getline(f, s)) {
    if (!s.empty() && s.back() == '\r') s.pop_back();
    if (!s.empty()) return s;
  }
  fail("truncated VTK header");
}
// Legacy VTK binary numbers are big endian, independent of the host.
double number(std::ifstream &f, const std::string &type)
{
  const int n = type == "double" ? 8 : type == "float" ? 4 : 0;
  if (!n) fail("only float/double VTK arrays are supported");
  unsigned char b[8];
  if (!f.read(reinterpret_cast<char *>(b), n)) fail("truncated binary VTK array");
  uint64_t bits = 0;
  for (int i=0; i<n; ++i) bits = (bits << 8) | b[i];
  double value;
  if (n == 8) std::memcpy(&value, &bits, 8);
  else { uint32_t small = bits; float v; std::memcpy(&v, &small, 4); value = v; }
  if (!std::isfinite(value)) fail("non-finite VTK coordinate or velocity");
  return value;
}
OGRSpatialReference reference(int epsg)
{
  OGRSpatialReference s;
  s.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
  if (s.importFromEPSG(epsg) != OGRERR_NONE || !s.IsProjected()
      || std::abs(s.GetLinearUnits()-1) > 1e-9)
    fail("epsg must identify a projected metre-based CRS");
  return s;
}
void sameCRS(const OGRSpatialReference *actual, OGRSpatialReference &expected)
{
  if (!actual || !actual->IsProjected() || std::abs(actual->GetLinearUnits()-1)>1e-9)
    fail("reference data requires a projected metre-based CRS");
  OGRSpatialReference copy(*actual);
  copy.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
  if (!copy.IsSame(&expected)) fail("coordinate-system mismatch with declared VTK EPSG");
}
}

WindVTKGrid::WindVTKGrid(const std::string &file)
{
  std::ifstream f(file, std::ios::binary);
  if (!f) fail("cannot open " + file);
  if (line(f).rfind("# vtk DataFile Version",0) != 0) fail("not a legacy VTK file");
  std::string title;
  if (!std::getline(f,title)) fail("truncated VTK title"); // a blank title is legal
  if (line(f) != "BINARY" || line(f) != "DATASET STRUCTURED_GRID")
    fail("requires binary legacy STRUCTURED_GRID");
  std::string tag, type, extra;
  int ni=0, nj=0, nk=0;
  std::istringstream dims(line(f));
  if (!(dims >> tag >> ni >> nj >> nk) || tag != "DIMENSIONS" || ni<2 || nj<2 || nk<2 || (dims >> extra))
    fail("invalid DIMENSIONS (at least 2 per axis required)");
  if (uint64_t(ni)*nj > std::numeric_limits<uint64_t>::max()/uint64_t(nk))
    fail("DIMENSIONS overflow");
  const uint64_t count = uint64_t(ni)*nj*nk;
  // Check the file can actually contain the arrays before allocating from its header.
  auto pos = f.tellg(); f.seekg(0,std::ios::end); auto bytes = f.tellg(); f.seekg(pos);
  if (count > uint64_t(bytes)/24 || count > std::numeric_limits<size_t>::max()/sizeof(std::array<double,3>))
    fail("DIMENSIONS exceed file size");
  size_t declared=0;
  std::istringstream points(line(f));
  if (!(points >> tag >> declared >> type) || tag != "POINTS" || declared != count || (points >> extra))
    fail("POINTS count does not match DIMENSIONS");
  xs.resize(ni); ys.resize(nj); zs.resize(count); levels=nk;
  for (int k=0;k<nk;++k) for (int j=0;j<nj;++j) for (int i=0;i<ni;++i) {
    const size_t id=i+size_t(ni)*(j+size_t(nj)*k);
    const double px=number(f,type), py=number(f,type), pz=number(f,type);
    if (j==0 && k==0) xs[i]=px;
    if (i==0 && k==0) ys[j]=py;
    if (std::abs(px-xs[i])>1e-7 || std::abs(py-ys[j])>1e-7)
      fail("x/y must be separable vertical columns; curvilinear horizontal grids unsupported");
    zs[id]=pz;
    if (k && pz<=zs[id-size_t(ni)*nj]) fail("folded or non-increasing vertical coordinates");
  }
  for (const auto *axis : {&xs,&ys})
    for (size_t i=1;i<axis->size();++i)
      if ((*axis)[i]<=(*axis)[i-1]) fail("x/y axes must be strictly increasing");
  std::istringstream data(line(f));
  if (!(data >> tag >> declared) || tag != "POINT_DATA" || declared != count || (data >> extra))
    fail("requires matching POINT_DATA");
  std::string name;
  std::istringstream vectors(line(f));
  if (!(vectors >> tag >> name >> type) || tag != "VECTORS" || name != "wind_vectors" || (vectors >> extra))
    fail("requires point-data VECTORS wind_vectors float/double");
  velocity.resize(count);
  for (auto &v : velocity) for (double &component : v) {
    component=number(f,type);
    if (std::abs(component)>std::numeric_limits<float>::max()) fail("velocity exceeds QES float range");
  }
  // Refuse unexpected arrays rather than accidentally interpreting a different layout.
  while (std::getline(f,extra))
    if (extra.find_first_not_of(" \t\r") != std::string::npos) fail("unsupported trailing VTK data");
  std::cout << "[windVTK] Loaded " << ni << " x " << nj << " x " << nk << " physical 3D points\n";
}

WindVTKGrid::Column WindVTKGrid::column(double x, double y) const
{
  if (!std::isfinite(x) || !std::isfinite(y) || x<xs.front() || x>xs.back() || y<ys.front() || y>ys.back())
    fail("missing horizontal coverage at ("+std::to_string(x)+", "+std::to_string(y)+")");
  const auto bracket=[](const std::vector<double> &a,double v) {
    return std::min(size_t(std::upper_bound(a.begin(),a.end(),v)-a.begin()-1),a.size()-2);
  };
  const size_t i=bracket(xs,x), j=bracket(ys,y), id=i+j*xs.size();
  const double a=(x-xs[i])/(xs[i+1]-xs[i]), b=(y-ys[j])/(ys[j+1]-ys[j]);
  return {{id,id+1,id+xs.size(),id+xs.size()+1}, {(1-a)*(1-b),a*(1-b),(1-a)*b,a*b}};
}
double WindVTKGrid::height(const Column &c,int k) const
{
  double z=0;
  for (int n=0;n<4;++n) z+=c.weights[n]*zs[c.ids[n]+size_t(k)*xs.size()*ys.size()];
  return z;
}
double WindVTKGrid::ground(double x,double y) const { return height(column(x,y),0); }
double WindVTKGrid::groundMaximum(double xmin,double ymin,double xmax,double ymax) const
{
  std::vector<double> cx={xmin,xmax},cy={ymin,ymax};
  for (double x:xs) if (x>xmin && x<xmax) cx.push_back(x);
  for (double y:ys) if (y>ymin && y<ymax) cy.push_back(y);
  double maximum=-std::numeric_limits<double>::infinity();
  for (double y:cy) for (double x:cx) maximum=std::max(maximum,ground(x,y));
  return maximum;
}
double WindVTKGrid::top(double x,double y) const { return height(column(x,y),levels-1); }
std::array<double,3> WindVTKGrid::sample(double x,double y,double z) const
{
  const auto c=column(x,y);
  if (!std::isfinite(z) || z<height(c,0) || z>height(c,levels-1))
    fail("missing vertical coverage at ("+std::to_string(x)+", "+std::to_string(y)+", "+std::to_string(z)+")");
  int lo=0,hi=levels-1;
  while (hi-lo>1) { const int mid=(lo+hi)/2; if (z<height(c,mid)) hi=mid; else lo=mid; }
  const double a=(z-height(c,lo))/(height(c,hi)-height(c,lo));
  std::array<double,3> result{};
  for (int n=0;n<4;++n) for (int component=0;component<3;++component)
    result[component]+=c.weights[n]*((1-a)*velocity[c.ids[n]+size_t(lo)*xs.size()*ys.size()][component]
                                   +a*velocity[c.ids[n]+size_t(hi)*xs.size()*ys.size()][component]);
  for (double v:result) if (!std::isfinite(v)) fail("non-finite interpolated velocity");
  return result;
}

void WindVTKGrid::validateTerrain(const std::string &file,int epsg,double tolerance,
                                  double xmin,double ymin,double xmax,double ymax) const
{
  GDALAllRegister();
  Dataset ds(static_cast<GDALDataset *>(GDALOpen(file.c_str(),GA_ReadOnly)),GDALClose);
  if (!ds || ds->GetRasterCount()!=1) fail("requires single-band reference terrain TIFF");
  auto expected=reference(epsg); sameCRS(ds->GetSpatialRef(),expected);
  double gt[6];
  if (ds->GetGeoTransform(gt)!=CE_None || gt[1]<=0 || gt[5]>=0 || gt[2]!=0 || gt[4]!=0)
    fail("reference terrain must be north-up with valid metre geotransform");
  for (double coordinate:gt) if (!std::isfinite(coordinate)) fail("invalid terrain geotransform");
  const int nx=ds->GetRasterXSize(),ny=ds->GetRasterYSize();
  if (nx<2 || ny<2 || xmin<gt[0]+.5*gt[1] || xmax>gt[0]+(nx-.5)*gt[1]
      || ymax>gt[3]+.5*gt[5] || ymin<gt[3]+(ny-.5)*gt[5])
    fail("QES domain lacks reference terrain pixel-centre coverage");
  auto *band=ds->GetRasterBand(1);
  int hasNodata=0; const double nodata=band->GetNoDataValue(&hasNodata);
  std::vector<double> raster(size_t(nx)*ny);
  if (band->RasterIO(GF_Read,0,0,nx,ny,raster.data(),nx,ny,GDT_Float64,0,0)!=CE_None)
    fail("cannot read reference terrain");
  for (double v:raster)
    if (!std::isfinite(v) || (hasNodata && v==nodata)) fail("invalid or nodata reference terrain");
  int scaleSet=0,offsetSet=0;
  const double scale=band->GetScale(&scaleSet),offset=band->GetOffset(&offsetSet);
  if ((scaleSet && !std::isfinite(scale)) || (offsetSet && !std::isfinite(offset)))
    fail("invalid terrain scale/offset");
  size_t checked=0; double maxError=0;
  for (double y:ys) for (double x:xs) {
    // Validate source nodes with DEM interpolation support. Outer half-pixels
    // have no bilinear support and are never extrapolated. QES is cropped inside.
    const double px=(x-gt[0])/gt[1]-.5,py=(y-gt[3])/gt[5]-.5;
    if (px<0 || py<0 || px>nx-1 || py>ny-1) continue;
    const int i=std::min(int(px),nx-2),j=std::min(int(py),ny-2);
    const double a=px-i,b=py-j;
    double z=0;
    for (int dy=0;dy<2;++dy) for (int dx=0;dx<2;++dx) {
      double v=raster[i+dx+size_t(nx)*(j+dy)];
      if (!std::isfinite(v) || (hasNodata && v==nodata)) fail("invalid or nodata reference terrain");
      v=v*(scaleSet?scale:1)+(offsetSet?offset:0);
      z+=v*(dx?a:1-a)*(dy?b:1-b);
    }
    maxError=std::max(maxError,std::abs(z-ground(x,y))); ++checked;
  }
  if (checked<4) fail("insufficient shared terrain nodes to verify registration");
  if (maxError>tolerance) fail("terrain/vertical-datum mismatch: max error "+std::to_string(maxError)+" m");
  std::cout << "[windVTK] Reference EPSG:" << epsg << "; terrain registration nodes=" << checked
            << ", max error=" << maxError << " m; VTK CRS/units explicitly asserted by XML\n";
}
void WindVTKGrid::validateVectorCRS(const std::string &file,int epsg) const
{
  GDALAllRegister();
  Dataset ds(static_cast<GDALDataset *>(GDALOpenEx(file.c_str(),GDAL_OF_VECTOR|GDAL_OF_READONLY,nullptr,nullptr,nullptr)),GDALClose);
  if (!ds || ds->GetLayerCount()<1) fail("cannot open building shapefile for CRS check");
  auto expected=reference(epsg);
  for (int i=0;i<ds->GetLayerCount();++i) sameCRS(ds->GetLayer(i)->GetSpatialRef(),expected);
}
