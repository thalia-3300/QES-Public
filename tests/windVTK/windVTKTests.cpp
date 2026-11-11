#include "winds/WindVTKGrid.h"
#include "winds/WindProfilerVTK.h"
#include "winds/WINDSGeneralData.h"
#include "winds/WINDSOutputWorkspace.h"
#include "winds/SolverFactory.h"
#include <gdal_priv.h>
#include <ogr_spatialref.h>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>

namespace fs=std::filesystem;
int checks=0;
void check(bool value,const std::string &message)
{
  ++checks;
  if (!value) throw std::runtime_error(message);
}
void rejects(const std::function<void()> &f,const std::string &expected)
{
  try { f(); } catch (const std::exception &e) {
    check(std::string(e.what()).find(expected)!=std::string::npos,"wrong rejection: "+std::string(e.what()));
    return;
  }
  throw std::runtime_error("expected rejection: "+expected);
}
constexpr double ox=636000,oy=4780000;
double ground(double x,double y) { return 101+.08*x+.06*y+.002*x*y; }
std::array<double,3> affine(double x,double y,double z)
{
  return {3+.03*x-.02*y+.005*(z-100), -1+.01*x+.025*y-.002*(z-100), .4+.003*x-.004*y+.006*(z-100)};
}
void binary(std::ofstream &f,double d,bool single)
{
  uint64_t bits;
  if (single) { float v=d; uint32_t b; std::memcpy(&b,&v,4); bits=b; }
  else std::memcpy(&bits,&d,8);
  for (int i=(single?3:7);i>=0;--i) f.put(char((bits>>(8*i))&255));
}
void vtk(const fs::path &file,bool single=false,bool invalid=false,bool folded=false,bool warped=false)
{
  std::ofstream f(file,std::ios::binary);
  f << "# vtk DataFile Version 3.0\nManufactured terrain-following affine vector field\nBINARY\n\nDATASET STRUCTURED_GRID\nDIMENSIONS 4 4 4\nPOINTS 64 " << (single?"float":"double") << '\n';
  const double levels[]={0,.15,.45,1};
  for (int k=0;k<4;++k) for (int j=0;j<4;++j) for (int i=0;i<4;++i) {
    const double g=ground(i*10,j*10),z=g+levels[k]*(200-g);
    binary(f,ox+i*10+(warped&&k==2?1:0),single); binary(f,oy+j*10,single);
    binary(f,folded&&k==2?g:z,single);
  }
  f << "\nPOINT_DATA 64\nVECTORS wind_vectors " << (single?"float":"double") << '\n';
  for (int k=0;k<4;++k) for (int j=0;j<4;++j) for (int i=0;i<4;++i) {
    const double g=ground(i*10,j*10),z=g+levels[k]*(200-g);
    auto v=affine(i*10,j*10,z);
    if (invalid && k==1 && j==1 && i==1) v[2]=std::numeric_limits<double>::quiet_NaN();
    for (double c:v) binary(f,c,single);
  }
  f << '\n';
}
void dem(const fs::path &file,int epsg=32612,double offset=0,bool nodata=false)
{
  GDALAllRegister();
  auto *driver=GetGDALDriverManager()->GetDriverByName("GTiff");
  std::unique_ptr<GDALDataset,decltype(&GDALClose)> ds(driver->Create(file.c_str(),41,41,1,GDT_Float64,nullptr),GDALClose);
  check(bool(ds),"create terrain fixture");
  double gt[]={ox-5.5,1,0,oy+35.5,0,-1}; ds->SetGeoTransform(gt);
  OGRSpatialReference s; s.importFromEPSG(epsg); ds->SetSpatialRef(&s);
  std::vector<double> values(41*41);
  for (int j=0;j<41;++j) for (int i=0;i<41;++i) values[i+41*j]=ground(i-5,35-j)+offset;
  ds->GetRasterBand(1)->SetNoDataValue(-9999);
  if (nodata) values[15+41*15]=-9999;
  check(ds->GetRasterBand(1)->RasterIO(GF_Write,0,0,41,41,values.data(),41,41,GDT_Float64,0,0)==CE_None,"write terrain");
}
std::string xml(const fs::path &dir,bool building)
{
  std::ostringstream s;
  s << "<QESWindsParameters><simulationParameters><domain>12 10 18</domain><cellSize>2 2 5</cellSize>"
       "<totalTimeIncrements>1</totalTimeIncrements><tolerance>1e-9</tolerance><maxIterations>10000</maxIterations>"
       "</simulationParameters><metParams><windVTK><file>" << (dir/"field.vtk").string() << "</file><terrainFile>"
    << (dir/"terrain.tif").string() << "</terrainFile><epsg>32612</epsg><velocityUnits>m/s</velocityUnits>"
       "<timeStamp>2020-01-01T00:00:00</timeStamp><originX>636002</originX><originY>4780002</originY>"
       "<originZ>100</originZ></windVTK></metParams>";
  if (building) s << "<buildingsParams><wallRoughness>0.01</wallRoughness><rooftopFlag>1</rooftopFlag>"
       "<upwindCavityFlag>2</upwindCavityFlag><streetCanyonFlag>1</streetCanyonFlag><streetIntersectionFlag>0</streetIntersectionFlag>"
       "<wakeFlag>2</wakeFlag><highRiseFlag>0</highRiseFlag><sidewallFlag>1</sidewallFlag>"
       "<rectangularBuilding><height>15</height><baseHeight>0</baseHeight><xStart>8</xStart><yStart>8</yStart>"
       "<length>4</length><width>4</width><buildingRotation>0</buildingRotation></rectangularBuilding></buildingsParams>";
  s << "<fileOptions><outputFlag>0</outputFlag><outputFields>all</outputFields></fileOptions></QESWindsParameters>";
  return s.str();
}
void write(const fs::path &p,const std::string &s) { std::ofstream f(p); f<<s; }
void checkFields(const WINDSGeneralData &w,bool manufactured)
{
  const auto &d=w.domain;
  const std::vector<float> *a[]={&w.u0_imported,&w.v0_imported,&w.w0_imported};
  const std::vector<int> *mask[]={&w.vtk_u_mask,&w.vtk_v_mask,&w.vtk_w_mask};
  double maxError=0; size_t nonzeroW=0;
  for (int c=0;c<3;++c) for (int k=0;k<d.nz();++k) for (int j=0;j<d.ny();++j) for (int i=0;i<d.nx();++i) {
    const auto id=d.face(i,j,k); const float value=(*a[c])[id];
    check(std::isfinite(value),"finite imported face");
    if ((*mask[c])[id]==0) check(value==0,"zero imported solid face");
    else if ((*mask[c])[id]==1) {
      if (c==2 && value!=0) ++nonzeroW;
      if (manufactured) {
        const double x=w.vtk_origin_x-ox+(i+(c==0?0:.5))*d.dx();
        const double y=w.vtk_origin_y-oy+(j+(c==1?0:.5))*d.dy();
        const double z=w.vtk_origin_z+(c==2?d.z_face[k]:d.z[k]);
        maxError=std::max(maxError,std::abs(value-affine(x,y,z)[c]));
      }
    } else {
      check((*mask[c])[id]==2,"known mask code");
      if (k==0) check(value==0,"zero ground ghost");
      else {
        const auto source=d.face(c==0?i:std::min(i,d.nx()-2),c==1?j:std::min(j,d.ny()-2),c==2?k:std::min(k,d.nz()-2));
        check(value==(*a[c])[source],"padding copies physical face");
      }
    }
  }
  check(nonzeroW>0,"preserve nonzero imported w");
  if (manufactured) check(maxError<5e-7,"affine physical-XYZ staggered interpolation error");
  if (manufactured) std::cout << "[test] Staggered affine max error=" << maxError << " m/s\n";
  std::cout << "[test] Nonzero imported w faces=" << nonzeroW << '\n';
}
void runCase(const fs::path &xmlPath,const fs::path &output,bool manufactured,bool building,int solveType)
{
  WINDSInputData wid(xmlPath.string());
  qes::Domain d(wid.simParams->domain,wid.simParams->grid);
  WINDSGeneralData w(&wid,d,solveType);
  {
    // Register output before import, matching qesWindsMain's lifecycle.
    WINDSOutputWorkspace saved(&w,output.string());
    w.resetICellFlag(); w.applyWindProfile(&wid,0,solveType);
    checkFields(w,manufactured);
    const auto u=w.u0_imported,v=w.v0_imported,z=w.w0_imported;
    w.applyParametrizations(&wid);
    check(w.u0_imported==u && w.v0_imported==v && w.w0_imported==z,"immutable pre-parameterisation snapshots");
    if (building) check(w.u0!=u || w.v0!=v || w.w0!=z,"building parameterisations change input");
    else check(w.u0==u && w.v0==v && w.w0==z,"no-building import unchanged by parameterisations");
    SolverFactory factory;
    std::unique_ptr<Solver> solver(factory.create(solveType,d,wid.simParams->tolerance));
    solver->solve(&w,wid.simParams->maxIterations);
    check(solver->diagnostics.converged,"wind solve converged");
    check(solver->diagnostics.maxDivergence<2e-6,"final solved-cell divergence <2e-6 s^-1");
    const std::vector<float> *final[]={&w.u,&w.v,&w.w};
    const std::vector<int> *mask[]={&w.vtk_u_mask,&w.vtk_v_mask,&w.vtk_w_mask};
    for (int c=0;c<3;++c) for (size_t i=0;i<final[c]->size();++i) {
      check(std::isfinite((*final[c])[i]),"finite corrected velocity");
      if ((*mask[c])[i]==0) check((*final[c])[i]==0,"zero final solid/terrain face");
    }
    // Physical top w remains prescribed; tangential/upper padding must follow
    // the corrected physical faces rather than retain stale imported values.
    for (int c=0;c<3;++c) for (int k=1;k<d.nz();++k) for (int j=0;j<d.ny();++j) for (int i=0;i<d.nx();++i) {
      const auto id=d.face(i,j,k);
      if ((*mask[c])[id]==2) {
        const auto src=d.face(c==0?i:std::min(i,d.nx()-2),c==1?j:std::min(j,d.ny()-2),c==2?k:std::min(k,d.nz()-2));
        check((*final[c])[id]==(*final[c])[src],"final padding copies corrected face");
      }
      if (c==2 && k==d.nz()-1 && (*mask[c])[id]==1)
        check(w.w[id]==w.w0[id],"prescribed physical top w face");
    }
    saved.save(w.timestamp[0]);
  }
  netCDF::NcFile nc(output.string(),netCDF::NcFile::read);
  for (const auto &entry : std::vector<std::pair<std::string,const std::vector<float> *>>{
      {"u0_imported",&w.u0_imported},{"v0_imported",&w.v0_imported},{"w0_imported",&w.w0_imported},
      {"u0",&w.u0},{"v0",&w.v0},{"w0",&w.w0},{"u",&w.u},{"v",&w.v},{"w",&w.w}}) {
    auto var=nc.getVar(entry.first); check(!var.isNull(),"saved "+entry.first);
    std::vector<float> data(entry.second->size()); var.getVar(data.data());
    check(data==*entry.second,"saved exact "+entry.first);
  }
  double origin=0; nc.getAtt("vtk_origin_z").getValues(&origin); check(origin==w.vtk_origin_z,"saved vertical datum origin");
  std::cout << "[test] Saved " << output << '\n';
}
void wakeHeightTranslation(const fs::path &dir)
{
  // Changing only the vertical coordinate datum must not change a wake.
  // Exercise the production polygonWake routine without other building
  // parameterisations or the mass-consistency solver obscuring its response.
  const auto path=dir/"wake-height.xml";
  write(path,xml(dir,true));
  struct Sample { std::vector<float> u,v,w; std::vector<int> flags; float length; };
  for (int rooftop : {1,2}) for (int highRise : {0,1}) {
    const auto sample=[&](float offset) {
      WINDSInputData wid(path.string());
      wid.buildingsParams->rooftopFlag=rooftop;
      wid.buildingsParams->highRiseFlag=highRise;
      qes::Domain d(wid.simParams->domain,wid.simParams->grid);
      WINDSGeneralData w(&wid,d,1);
      check(w.allBuildingsV.size()==1,"one wake-height fixture building");
      auto *building=w.allBuildingsV.front();
      building->base_height+=offset;
      building->height_eff+=offset;
      for (auto &z : w.domain.z) z+=offset;
      for (auto &z : w.domain.z_face) z+=offset;
      std::fill(w.u0.begin(),w.u0.end(),0.0f);
      std::fill(w.v0.begin(),w.v0.end(),2.0f);
      std::fill(w.w0.begin(),w.w0.end(),0.0f);
      building->polygonWake(&wid,&w,0);
      check(std::count(w.icellflag.begin(),w.icellflag.end(),4)
            +std::count(w.icellflag.begin(),w.icellflag.end(),5)>0,"nonempty wake-height fixture wake");
      return Sample{w.u0,w.v0,w.w0,w.icellflag,building->Lr};
    };
    const auto original=sample(0),shifted=sample(100);
    double error=0;
    for (const auto &pair : {std::make_pair(&original.u,&shifted.u),
                             std::make_pair(&original.v,&shifted.v),
                             std::make_pair(&original.w,&shifted.w)}) {
      for (size_t i=0;i<pair.first->size();++i) {
        check(std::isfinite((*pair.first)[i]) && std::isfinite((*pair.second)[i]),"finite wake-height velocities");
        error=std::max(error,double(std::abs((*pair.first)[i]-(*pair.second)[i])));
      }
    }
    std::cout << "[test] Wake vertical-datum shift: rooftop=" << rooftop
              << " highRise=" << highRise << " maxVelocityDifference=" << error
              << " Lr=" << original.length << "/" << shifted.length << '\n';
    check(error<5e-5,"wake velocities invariant to a 100 m vertical datum shift");
    check(original.flags==shifted.flags,"wake region invariant to a vertical datum shift");
    check(std::abs(original.length-shifted.length)<5e-5,"wake length invariant to a vertical datum shift");
  }
}
void manufactured(const fs::path &dir)
{
  vtk(dir/"field.vtk"); dem(dir/"terrain.tif");
  WindVTKGrid grid((dir/"field.vtk").string());
  double maxError=0;
  for (int n=0;n<200;++n) {
    const double x=.05+29.9*((n*37)%199)/199.,y=.05+29.9*((n*61)%199)/199.;
    const double z=ground(x,y)+(.001+.998*((n*17)%199)/199.)*(200-ground(x,y));
    const auto actual=grid.sample(ox+x,oy+y,z),expected=affine(x,y,z);
    for (int c=0;c<3;++c) maxError=std::max(maxError,std::abs(actual[c]-expected[c]));
  }
  check(maxError<1e-9,"physical affine interpolation across sloping layers");
  std::cout << "[test] Double-precision physical affine interpolation max error=" << maxError << " m/s\n";
  const auto bottom=grid.sample(ox,oy,ground(0,0)),top=grid.sample(ox+30,oy+30,200);
  check(std::abs(bottom[2]-affine(0,0,ground(0,0))[2])<1e-12 && std::abs(top[2]-affine(30,30,200)[2])<1e-12,"endpoint sampling");
  rejects([&]{grid.sample(ox-.01,oy+5,150);},"horizontal coverage");
  rejects([&]{grid.sample(ox+5,oy+5,200.01);},"vertical coverage");
  rejects([&]{grid.sample(ox+5,oy+5,100);},"vertical coverage");
  rejects([&]{grid.sample(ox+5,oy+5,std::numeric_limits<double>::quiet_NaN());},"vertical coverage");
  vtk(dir/"float.vtk",true); WindVTKGrid single((dir/"float.vtk").string());
  check(std::abs(single.sample(ox+15,oy+15,150)[2]-affine(15,15,150)[2])<1e-6,"big-endian float decoding");
  vtk(dir/"invalid.vtk",false,true); rejects([&]{WindVTKGrid bad((dir/"invalid.vtk").string());},"non-finite");
  vtk(dir/"folded.vtk",false,false,true); rejects([&]{WindVTKGrid bad((dir/"folded.vtk").string());},"non-increasing vertical");
  vtk(dir/"warped.vtk",false,false,false,true); rejects([&]{WindVTKGrid bad((dir/"warped.vtk").string());},"separable");
  fs::copy_file(dir/"field.vtk",dir/"truncated.vtk"); fs::resize_file(dir/"truncated.vtk",2500);
  rejects([&]{WindVTKGrid bad((dir/"truncated.vtk").string());},"truncated");
  dem(dir/"wrong-crs.tif",32613); dem(dir/"wrong-z.tif",32612,5); dem(dir/"nodata.tif",32612,0,true);
  const auto validate=[&](const std::string &file,int epsg=32612) { grid.validateTerrain((dir/file).string(),epsg,.1,ox+2,oy+2,ox+26,oy+22); };
  validate("terrain.tif");
  rejects([&]{validate("wrong-crs.tif");},"coordinate-system mismatch");
  rejects([&]{validate("wrong-z.tif");},"vertical-datum mismatch");
  rejects([&]{validate("nodata.tif");},"nodata");
  rejects([&]{validate("terrain.tif",4326);},"projected metre-based");
  const auto path=dir/"case.xml";
  write(path,xml(dir,false)); runCase(path,dir/"manufactured.nc",true,false,1);
  write(path,xml(dir,true)); runCase(path,dir/"building.nc",true,true,1);
  wakeHeightTranslation(dir);
  for (const auto &replacement:std::vector<std::pair<std::string,std::string>>{
      {"<epsg>32612</epsg>","<epsg>32613</epsg>"},
      {"<originX>636002</originX>","<originX>636020</originX>"},
      {"<originZ>100</originZ>","<originZ>120</originZ>"},
      {"<domain>12 10 18</domain>","<domain>12 10 22</domain>"},
      {"</windVTK>","<misspelled>1</misspelled></windVTK>"},
      {"</windVTK>","<epsg>32612</epsg></windVTK>"},
      {"</simulationParameters>","<meshTypeFlag>1</meshTypeFlag></simulationParameters>"},
      {"</simulationParameters>","<UTMZone>13</UTMZone></simulationParameters>"},
      {"</simulationParameters>","<WRFCoupling>true</WRFCoupling></simulationParameters>"},
      {"</simulationParameters>","<DEMDistancex>1</DEMDistancex></simulationParameters>"},
      {"</metParams>","<sensorName>unused.xml</sensorName></metParams>"}}) {
    auto bad=xml(dir,false); bad.replace(bad.find(replacement.first),replacement.first.size(),replacement.second);
    write(path,bad); rejects([&]{WINDSInputData wid(path.string());},"windVTK");
  }
}
int main(int argc,char **argv)
{
  try {
    if (argc>=3 && std::string(argv[1])=="--integration") {
      const fs::path output=argc>=4?argv[3]:(fs::temp_directory_path()/("ghosttown_vtk_"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count())+".nc")).string();
      runCase(argv[2],output,false,false,argc>=5?std::stoi(argv[4]):1);
    } else if (argc>=2 && std::string(argv[1])=="--wake-height-only") {
      const auto dir=fs::temp_directory_path()/("qes_wake_height_"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
      fs::create_directories(dir);vtk(dir/"field.vtk");dem(dir/"terrain.tif");
      wakeHeightTranslation(dir);
    } else {
      const auto dir=fs::temp_directory_path()/("qes_wind_vtk_"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
      fs::create_directories(dir); manufactured(dir);
      std::cout << "[test] Fixtures: " << dir << '\n';
    }
    std::cout << "PASS: " << checks << " checks\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "FAIL: " << e.what() << '\n'; return 1;
  }
}
