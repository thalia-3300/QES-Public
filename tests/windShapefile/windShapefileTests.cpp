#include "winds/WindShapefileGrid.h"
#include "winds/WindProfilerShapefile.h"
#include "winds/WINDSGeneralData.h"
#include "winds/SolverFactory.h"
#include "winds/Solver_CPU.h"
#include "winds/Solver_CPU_RB.h"
#include <ogrsf_frmts.h>
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>

namespace fs = std::filesystem;
using Point = std::array<double, 4>; // x,y,speed,FROM direction
int checks = 0;
void check(bool ok, const std::string &label)
{
  ++checks;
  if (!ok) throw std::runtime_error(label);
}
void near(double value, double expected, const std::string &label, double epsilon = 1e-5)
{
  check(std::isfinite(value) && std::abs(value-expected) < epsilon, label);
}
void rejects(const std::function<void()> &f, const std::string &message)
{
  try { f(); } catch (const std::exception &e) {
    check(std::string(e.what()).find(message) != std::string::npos,
          "Wrong rejection: " + std::string(e.what()) + " (expected " + message + ")");
    return;
  }
  check(false, "Expected rejection: " + message);
}

void shape(const fs::path &file, const std::vector<Point> &points, int epsg = 32755,
           bool nullSpeed = false, OGRwkbGeometryType type = wkbPoint)
{
  auto *driver = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
  std::unique_ptr<GDALDataset, decltype(&GDALClose)> ds(
    driver->Create(file.string().c_str(), 0,0,0,GDT_Unknown,nullptr), GDALClose);
  check(bool(ds), "create test shapefile");
  OGRSpatialReference srs;
  srs.importFromEPSG(epsg);
  auto *layer = ds->CreateLayer("wind", &srs, type, nullptr);
  OGRFieldDefn speed("velocity", OFTReal), direction("bearing", OFTReal);
  speed.SetWidth(24); speed.SetPrecision(12);
  direction.SetWidth(24); direction.SetPrecision(12);
  check(layer->CreateField(&speed)==OGRERR_NONE, "create speed field");
  check(layer->CreateField(&direction)==OGRERR_NONE, "create direction field");
  for (const auto &p : points) {
    std::unique_ptr<OGRFeature, decltype(&OGRFeature::DestroyFeature)> feature(
      OGRFeature::CreateFeature(layer->GetLayerDefn()), OGRFeature::DestroyFeature);
    if (!nullSpeed) feature->SetField("velocity", p[2]);
    feature->SetField("bearing", p[3]);
    if (type == wkbPoint) {
      OGRPoint geometry(p[0],p[1]); feature->SetGeometry(&geometry);
    } else {
      OGRLineString geometry;
      geometry.addPoint(p[0],p[1]); geometry.addPoint(p[0]+1,p[1]+1);
      feature->SetGeometry(&geometry);
    }
    check(layer->CreateFeature(feature.get()) == OGRERR_NONE, "create test feature");
  }
}

WindShapefileInput config(const fs::path &file)
{
  WindShapefileInput c;
  c.file = file.string(); c.height = 10; c.siteZ0 = 0.1;
  c.speedField = "velocity"; c.directionField = "bearing";
  c.speedUnits = "m/s"; c.timeStamp = "2026-09-01T00:00:00";
  c.originX = 100000; c.originY = 5000000;
  return c;
}

void synthetic(const fs::path &dir)
{
  std::vector<Point> points;
  for (int j = 0; j < 3; ++j) for (int i = 0; i < 3; ++i) {
    const double x = i*10, y = j*20;
    const double u = 2 + 0.1*x + 0.05*y, v = -3 + 0.02*x - 0.03*y;
    double direction = std::atan2(-u,-v)*180/std::acos(-1.0);
    if (direction < 0) direction += 360;
    points.push_back({100000+x,5000000+y,std::hypot(u,v),direction});
  }
  std::reverse(points.begin(), points.end()); // independent of feature ordering
  const auto file = dir/"affine.shp";
  shape(file, points);
  auto c = config(file);
  WindShapefileGrid grid(c);
  check(grid.nx()==3 && grid.ny()==3, "rectangular axes");
  for (double y : {0.,7.,20.,40.}) for (double x : {0.,3.,10.,20.}) {
    auto uv = grid.sample(100000+x,5000000+y);
    near(uv[0],2+.1*x+.05*y,"bilinear u");
    near(uv[1],-3+.02*x-.03*y,"bilinear v");
  }
  rejects([&]{grid.sample(99999,5000000);},"coverage");
  rejects([&]{grid.validateCoverage(100000,5000000,100021,5000040);},"coverage");
  rejects([&]{grid.sample(std::numeric_limits<double>::quiet_NaN(),5000000);},"coverage");
  for (auto unit : {"m/s","mps","mph","kph","km/h","kts","knots"}) {
    auto units = c; units.speedUnits = unit;
    const double factor = (units.speedUnits=="mph") ? .44704
        : (units.speedUnits=="kph" || units.speedUnits=="km/h") ? 1./3.6
        : (units.speedUnits=="kts" || units.speedUnits=="knots") ? 1852./3600 : 1.;
    WindShapefileGrid converted(units);
    near(converted.sample(100000,5000000)[0],2*factor,"speed units");
  }
  auto wrong = c; wrong.speedUnits = "bananas";
  rejects([&]{WindShapefileGrid invalid(wrong);},"speedUnits");
  wrong = c; wrong.speedField = "absent";
  rejects([&]{WindShapefileGrid invalid(wrong);},"missing speed/direction");
  wrong = c; wrong.height = 0;
  rejects([&]{wrong.validate();},"positive height");
  wrong = c; wrong.boundaryLayerFlag = 4;
  rejects([&]{wrong.validate();},"profiles supported");
  wrong = c; wrong.directionUnits = "grads";
  rejects([&]{wrong.validate();},"directionUnits");
  for (double stability : {0.,.01,-.01}) {
    c.reciprocal = stability;
    near(c.profileScale(c.height),1,"reference-height normalization");
    check(c.profileScale(30)>1 && c.profileScale(2)<1,"vertical-profile gradient");
    near(c.profileScale(-1),0,"subsurface wind");
  }
  c.reciprocal = 0;
  near(c.profileScale(20),std::log(201.)/std::log(101.),"neutral log profile");
  c.boundaryLayerFlag = 0; near(c.profileScale(20),1,"uniform profile");
  c.boundaryLayerFlag = 1;

  // Exercise actual QES face arrays, including non-square cell spacing and padding.
  SimulationParameters sim;
  MetParams met; met.windShapefile = &c;
  WINDSInputData wid; wid.simParams = &sim; wid.metParams = &met;
  qes::Domain domain(4,4,3,5,10,20);
  WINDSGeneralData wgd(domain);
  WindProfilerShapefile profiler(&wid,domain);
  profiler.interpolateWindProfile(&wid,&wgd);
  for (int j=0;j<domain.ny();++j) for (int i=0;i<domain.nx();++i) {
    const auto id = domain.face(i,j,1);
    const double u = 2+.1*(i*5)+.05*((std::min(j,3)+.5)*10);
    const double v = -3+.02*((std::min(i,3)+.5)*5)-.03*(j*10);
    near(wgd.u0[id],u,"staggered u reference height");
    near(wgd.v0[id],v,"staggered v reference height");
    near(wgd.u0[domain.face(i,j,2)],u*c.profileScale(30),"upper profile");
    near(wgd.u0[domain.face(i,j,0)],0,"bottom ghost");
    near(wgd.u0[domain.face(i,j,4)],wgd.u0[domain.face(i,j,3)],"top padding");
  }
  check(std::all_of(wgd.w0.begin(),wgd.w0.end(),[](float w){return w==0;}),"initial w zero");
  std::fill(wgd.terrain.begin(),wgd.terrain.end(),15);
  profiler.interpolateWindProfile(&wid,&wgd);
  near(wgd.u0[domain.face(1,1,1)],0,"terrain below ground");
  near(wgd.u0[domain.face(1,1,2)],3.25*c.profileScale(15),"AGL above terrain");
  rejects([&]{qes::Domain larger(5,4,3,5,10,20); WindProfilerShapefile invalid(&wid,larger);},"coverage");

  const auto bad = [&](std::vector<Point> p, const std::string &name, const std::string &message,
                       int epsg=32755, bool null=false, OGRwkbGeometryType type=wkbPoint) {
    const auto path = dir/(name+".shp"); shape(path,p,epsg,null,type);
    rejects([&]{WindShapefileGrid invalid(config(path));},message);
  };
  auto p = points; p.pop_back(); bad(p,"missing","incomplete");
  p = points; p[0]=p[1]; bad(p,"duplicate","duplicate");
  p = points; p[0][0]+=.2; bad(p,"irregular","irregular");
  p = points; p[0][2]=-1; bad(p,"negative","invalid coordinate");
  p = points; p[0][3]=361; bad(p,"angle","invalid coordinate");
  bad(points,"null","missing speed/direction value",32755,true);
  bad(points,"geographic","projected metre",4326);
  bad(points,"feet","projected metre",2227);
  bad(points,"lines","nonempty points",32755,false,wkbLineString);
  shape(dir/"no_crs.shp",points);
  fs::rename(dir/"no_crs.prj",dir/"unused.prj");
  rejects([&]{WindShapefileGrid invalid(config(dir/"no_crs.shp"));},"projected metre");
  grid.validateVectorCRS(file.string());
  rejects([&]{grid.validateVectorCRS((dir/"geographic.shp").string());},"CRS differ");
  // DEM coordinate validation without depending on private terrain files.
  const auto dem = dir/"terrain.tif";
  {
    auto *driver = GetGDALDriverManager()->GetDriverByName("GTiff");
    std::unique_ptr<GDALDataset, decltype(&GDALClose)> ds(
      driver->Create(dem.string().c_str(),6,6,1,GDT_Float32,nullptr),GDALClose);
    double transform[] = {100000,10,0,5000060,0,-10};
    ds->SetGeoTransform(transform);
    OGRSpatialReference srs; srs.importFromEPSG(32755); ds->SetSpatialRef(&srs);
  }
  grid.validateDEM(dem.string(),100000,5000000,0,0,0,0,0);
  grid.validateDEM(dem.string(),100010,5000010,1,0,0,100010,5000010);
  rejects([&]{grid.validateDEM(dem.string(),100001,5000000,0,1,0,0,0);},"pixel boundaries");
  rejects([&]{grid.validateDEM(dem.string(),100001,5000000,0,0,0,0,0);},"crop origin");
  shape(dir/"other_crs.shp",points,7854);
  WindShapefileGrid other(config(dir/"other_crs.shp"));
  rejects([&]{other.validateDEM(dem.string(),100000,5000000,0,0,0,0,0);},"CRS differ");
  // A 359/1 degree seam should point north, not average to south.
  p = {{100000,5000000,10,359},{100010,5000000,10,1},
       {100000,5000020,10,359},{100010,5000020,10,1}};
  shape(dir/"seam.shp",p);
  WindShapefileGrid seam(config(dir/"seam.shp"));
  near(seam.sample(100005,5000010)[0],0,"angle seam u");
  near(seam.sample(100005,5000010)[1],-10*std::cos(std::acos(-1.)/180),"angle seam v");
  for (auto &point : p) point[3] *= std::acos(-1.)/180;
  shape(dir/"radians.shp",p);
  c = config(dir/"radians.shp"); c.directionUnits = "radians";
  WindShapefileGrid radians(c);
  near(radians.sample(100005,5000010)[1],seam.sample(100005,5000010)[1],"radians");
  for (int angle : {0,90,180,270,360}) {
    for (auto &point : p) {point[2]=10; point[3]=angle;}
    const auto path=dir/("cardinal"+std::to_string(angle)+".shp"); shape(path,p);
    WindShapefileGrid cardinal(config(path));
    const auto uv=cardinal.sample(100005,5000010);
    near(uv[0],-10*std::sin(angle*std::acos(-1.)/180),"cardinal u");
    near(uv[1],-10*std::cos(angle*std::acos(-1.)/180),"cardinal v");
  }
  for (auto &point : p) point[2]=0;
  shape(dir/"calm.shp",p);
  WindShapefileGrid calm(config(dir/"calm.shp"));
  near(calm.sample(100005,5000010)[0],0,"calm u");
  near(calm.sample(100005,5000010)[1],0,"calm v");
}

void checkOutput(const std::string &path)
{
  netCDF::NcFile file(path,netCDF::NcFile::read);
  check(file.getDim("t").getSize()==1,"one output timestep");
  bool horizontalFlow = false;
  for (const auto &name : {"u","v","w"}) {
    auto var = file.getVar(name); check(!var.isNull(),"output velocity variable");
    size_t n=1; for (const auto &d : var.getDims()) n*=d.getSize();
    std::vector<float> values(n); var.getVar(values.data());
    check(std::all_of(values.begin(),values.end(),[](float v){return std::isfinite(v) && std::abs(v)<1e5;}),
          "finite output "+std::string(name));
    if (std::string(name)!="w")
      horizontalFlow |= std::any_of(values.begin(),values.end(),[](float v){return std::abs(v)>1e-6;});
  }
  check(horizontalFlow,"nonzero horizontal flow");
}

void checkTurbulence(const std::string &path)
{
  netCDF::NcFile file(path,netCDF::NcFile::read);
  check(file.getDim("t").getSize()==1,"one turbulence timestep");
  for (const auto &entry : file.getVars()) {
    const auto &var=entry.second;
    if (var.getType()!=netCDF::ncFloat) continue;
    size_t n=1; for (const auto &d : var.getDims()) n*=d.getSize();
    std::vector<float> values(n); var.getVar(values.data());
    check(std::all_of(values.begin(),values.end(),[](float v){return std::isfinite(v) && std::abs(v)<1e30;}),
          "finite turbulence "+entry.first);
  }
}

void compare(const std::string &first, const std::string &second, bool initial = false)
{
  netCDF::NcFile a(first,netCDF::NcFile::read), b(second,netCDF::NcFile::read);
  const std::vector<std::string> names = initial ? std::vector<std::string>{"u0","v0","w0"}
                                                : std::vector<std::string>{"u","v","w"};
  bool agrees = true;
  for (const auto &name : names) {
    auto av=a.getVar(name), bv=b.getVar(name);
    size_t n=1, m=1;
    for (const auto &d : av.getDims()) n*=d.getSize();
    for (const auto &d : bv.getDims()) m*=d.getSize();
    check(n==m,"comparison dimensions");
    std::vector<float> x(n),y(n); av.getVar(x.data()); bv.getVar(y.data());
    double maxDifference=0, magnitude=0;
    size_t worst=0;
    for (size_t i=0;i<n;++i) {
      check(std::isfinite(x[i]) && std::isfinite(y[i]),"finite comparison");
      if (std::abs(x[i]-y[i]) > maxDifference) { maxDifference=std::abs(x[i]-y[i]); worst=i; }
      magnitude=std::max(magnitude,double(std::abs(x[i])));
    }
    std::cout << name << " max CPU/GPU difference: " << maxDifference << " at flat index " << worst << '\n';
    agrees &= initial ? maxDifference == 0 : maxDifference < 2e-4*(1+magnitude);
  }
  check(agrees,"CPU/GPU velocity agreement");
}

void solverRegression(bool gpu)
{
  qes::Domain d(12,9,8,3,5,2);
  near(std::get<1>(d.getDomainSize()),5,"non-square dy");
  WINDSGeneralData w(d);
  std::fill(w.e.begin(),w.e.end(),1.0/(d.dx()*d.dx())); w.f=w.e;
  std::fill(w.g.begin(),w.g.end(),1.0/(d.dy()*d.dy())); w.h=w.g;
  std::fill(w.m.begin(),w.m.end(),1.0/(d.dz()*d.dz())); w.n=w.m;
  std::fill(w.icellflag.begin(),w.icellflag.end(),1);
  // QES does not scale the prescribed top row's coefficients. Reconstruction
  // must obtain its lower-face coefficient from the last solved row below.
  for(int j=0;j<d.ny()-1;++j) for(int i=0;i<d.nx()-1;++i) {
    const auto c=d.cell(i,j,d.nz()-2);
    w.e[c]=w.f[c]=w.g[c]=w.h[c]=w.m[c]=w.n[c]=1;
  }
  const auto phi=[&](int i,int j,int k) {
    const double top=d.nz()-2.5;
    return 1000.0*(i*(d.nx()-2-i)/30.0)*(j*(d.ny()-2-j)/16.0)
      *(top*top-(k-.5)*(k-.5))/(top*top);
  };
  std::fill(w.u0.begin(),w.u0.end(),5);
  std::fill(w.v0.begin(),w.v0.end(),1);
  std::fill(w.w0.begin(),w.w0.end(),0);
  for(int k=1;k<d.nz()-1;++k) for(int j=1;j<d.ny()-1;++j) for(int i=1;i<d.nx()-1;++i) {
    const auto f=d.face(i,j,k);
    w.u0[f]-=.5*(phi(i,j,k)-phi(i-1,j,k))/d.dx();
    w.v0[f]-=.5*(phi(i,j,k)-phi(i,j-1,k))/d.dy();
    w.w0[f]-=.5*(phi(i,j,k)-phi(i,j,k-1))/d.dz();
  }
  // A gradient perturbation of a known divergence-free flow must be removed.
  std::vector<std::unique_ptr<Solver>> solvers;
  solvers.emplace_back(new Solver_CPU(d,1e-9));
  solvers.emplace_back(new Solver_CPU_RB(d,1e-9));
#ifdef HAS_CUDA
  if(gpu) {
    SolverFactory factory;
    solvers.emplace_back(factory.create(3,d,1e-9));
    solvers.emplace_back(factory.create(4,d,1e-9));
  }
#else
  check(!gpu,"GPU solver test requires CUDA build");
#endif
  for(auto &solver:solvers) {
    solver->solve(&w,10000);
    check(solver->diagnostics.converged,"manufactured solution converged");
    check(solver->diagnostics.maxDivergence<2e-5,"manufactured solution mass conservation");
    for(int k=1;k<d.nz()-1;++k) for(int j=1;j<d.ny()-1;++j) for(int i=1;i<d.nx()-1;++i) {
      const auto f=d.face(i,j,k);
      near(w.u[f],5,"manufactured u",3e-5);
      near(w.v[f],1,"manufactured v",3e-5);
      near(w.w[f],0,"manufactured w (including top correction face)",3e-5);
    }
    solver->resetLambda();
    solver->solve(&w,1);
    check(!solver->diagnostics.converged,"iteration cap must not report convergence");
  }
}

void checkMass(const std::string &path,bool verify)
{
  qes::Domain d(path);
  WINDSGeneralData w(d);
  netCDF::NcFile file(path,netCDF::NcFile::read);
  file.getVar("u").getVar(w.u.data());
  file.getVar("v").getVar(w.v.data());
  file.getVar("w").getVar(w.w.data());
  file.getVar("icellflag").getVar(w.icellflag.data());
  if(verify) {
    for(int k=1;k<d.nz()-1;++k) for(int j=0;j<d.ny()-1;++j) for(int i=0;i<d.nx()-1;++i) {
      const int flag=w.icellflag[d.cell(i,j,k)];
      if(flag!=0 && flag!=2) continue;
      const auto f=d.face(i,j,k);
      check(w.u[f]==0 && w.u[f+1]==0,"zero solid x-face velocity");
      check(w.v[f]==0 && w.v[f+d.nx()]==0,"zero solid y-face velocity");
      check(w.w[f]==0 && w.w[f+d.nx()*d.ny()]==0,"zero solid z-face velocity");
    }
  }
  double maximum=0,squares=0,netFlux=0; size_t count=0;
  for(int k=1;k<d.nz()-2;++k) for(int j=1;j<d.ny()-2;++j) for(int i=1;i<d.nx()-2;++i) {
    const int flag=w.icellflag[d.cell(i,j,k)];
    if(flag==0||flag==2||flag==7||flag==8) continue;
    const auto f=d.face(i,j,k);
    const double div=(double(w.u[f+1])-w.u[f])/d.dx()
      +(double(w.v[f+d.nx()])-w.v[f])/d.dy()
      +(double(w.w[f+d.nx()*d.ny()])-w.w[f])/d.dz_array[k];
    check(std::isfinite(div),"finite divergence");
    maximum=std::max(maximum,std::abs(div)); squares+=div*div; ++count;
    netFlux+=div*d.dx()*d.dy()*d.dz_array[k];
  }
  check(count>0,"mass-check air cells");
  std::cout << "DIVERGENCE max=" << maximum << " RMS=" << std::sqrt(squares/count)
            << " s^-1; net solved-region flux=" << netFlux << " m^3/s\n";
  if(verify) check(maximum<1e-5,"mass conservation max divergence <1e-5 s^-1");
}

int main(int argc, char **argv)
{
  try {
    GDALAllRegister();
    if(argc==2 && std::string(argv[1])=="--solver-tests") {
      solverRegression(false);
    } else if(argc==2 && std::string(argv[1])=="--solver-tests-gpu") {
      solverRegression(true);
    } else if(argc==3 && std::string(argv[1])=="--check-mass") {
      checkMass(argv[2],true);
    } else if(argc==3 && std::string(argv[1])=="--report-mass") {
      checkMass(argv[2],false);
    } else if (argc==3 && std::string(argv[1])=="--check-output") {
      checkOutput(argv[2]);
    } else if (argc==3 && std::string(argv[1])=="--check-turb") {
      checkTurbulence(argv[2]);
    } else if (argc==4 && std::string(argv[1])=="--compare") {
      compare(argv[2],argv[3]);
    } else if (argc==4 && std::string(argv[1])=="--compare-initial") {
      compare(argv[2],argv[3],true);
    } else if (argc==4 && std::string(argv[1])=="--real") {
      auto c=config(argv[2]); c.speedField="speed"; c.directionField="dir";
      WindShapefileGrid grid(c);
      check(grid.nx()==std::stoul(argv[3]) && grid.ny()==std::stoul(argv[3]),"real source dimensions");
      auto b=grid.bounds();
      grid.validateCoverage(b[0],b[1],b[2],b[3]);
      for (int j=0;j<=100;++j) for(int i=0;i<=100;++i) {
        auto uv=grid.sample(b[0]+i*(b[2]-b[0])/100,b[1]+j*(b[3]-b[1])/100);
        check(std::isfinite(uv[0]) && std::isfinite(uv[1]),"real coverage sampling");
      }
      std::cout << "Real grid " << grid.nx() << " x " << grid.ny() << " accepted\n";
    } else if (argc==1) {
      auto dir=fs::temp_directory_path()/("qes-wind-shapefile-"+
        std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
      fs::create_directory(dir);
      std::cout << "Test fixtures: " << dir << '\n';
      synthetic(dir);
    } else throw std::runtime_error("Usage: windShapefileTests [--solver-tests | --solver-tests-gpu | --real file.shp size | --check-output file.nc | --check-turb file.nc | --check-mass workspace.nc | --report-mass workspace.nc | --compare first.nc second.nc | --compare-initial first.nc second.nc]");
    std::cout << "PASS: " << checks << " checks\n";
    return 0;
  } catch (const std::exception &e) {std::cerr << "FAIL: " << e.what() << '\n'; return 1;}
}
