#include "WindShapefileGrid.h"
#include <gdal_priv.h>
#include <ogrsf_frmts.h>
#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <sstream>
#include <limits>

namespace {
constexpr double tolerance = 1e-5; // metres; preserves small spacing at large eastings
using Dataset = std::unique_ptr<GDALDataset, decltype(&GDALClose)>;
using Feature = std::unique_ptr<OGRFeature, decltype(&OGRFeature::DestroyFeature)>;
void fail(const std::string &message) { throw std::runtime_error("windShapefile: " + message); }
std::vector<double> axis(std::vector<double> values)
{
  std::sort(values.begin(), values.end());
  std::vector<double> result;
  for (double value : values)
    if (result.empty() || value - result.back() > tolerance) result.push_back(value);
  if (result.size() < 2) fail("regular grid needs at least two distinct coordinates per axis");
  const double step = (result.back() - result.front()) / (result.size() - 1);
  for (size_t i = 0; i < result.size(); ++i)
    if (std::abs(result[i] - result.front() - i * step) > tolerance)
      fail("irregular point spacing is unsupported");
  return result;
}
size_t indexOf(const std::vector<double> &a, double x)
{
  const auto i = static_cast<size_t>(std::llround((x - a.front()) / (a[1] - a[0])));
  if (i >= a.size() || std::abs(a[i] - x) > tolerance) fail("point is off the regular grid");
  return i;
}
}

WindShapefileGrid::WindShapefileGrid(const WindShapefileInput &input)
{
  input.validate();
  GDALAllRegister();
  Dataset ds(static_cast<GDALDataset *>(GDALOpenEx(input.file.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY,
                                                nullptr, nullptr, nullptr)), GDALClose);
  if (!ds || ds->GetLayerCount() != 1) fail("cannot open a single-layer point shapefile: " + input.file);
  OGRLayer *layer = ds->GetLayer(0);
  const auto *srs = layer->GetSpatialRef();
  if (!srs || !srs->IsProjected() || std::abs(srs->GetLinearUnits() - 1) > 1e-9)
    fail("a projected metre-based CRS (.prj) is required; reproject first");
  char *wkt = nullptr;
  srs->exportToWkt(&wkt);
  projection = wkt;
  CPLFree(wkt);
  auto *definition = layer->GetLayerDefn();
  const int speed = definition->GetFieldIndex(input.speedField.c_str());
  const int direction = definition->GetFieldIndex(input.directionField.c_str());
  for (int field : {speed, direction}) {
    if (field < 0) fail("missing speed/direction attribute field");
    const auto type = definition->GetFieldDefn(field)->GetType();
    if (type != OFTInteger && type != OFTInteger64 && type != OFTReal) fail("wind fields must be numeric");
  }
  std::vector<std::array<double, 4>> points;
  std::vector<double> xvalues, yvalues;
  layer->ResetReading();
  while (Feature feature{layer->GetNextFeature(), OGRFeature::DestroyFeature}) {
    const auto *geometry = feature->GetGeometryRef();
    if (!geometry || geometry->IsEmpty() || wkbFlatten(geometry->getGeometryType()) != wkbPoint)
      fail("all features must be nonempty points");
    if (!feature->IsFieldSetAndNotNull(speed) || !feature->IsFieldSetAndNotNull(direction))
      fail("missing speed/direction value");
    const auto *point = geometry->toPoint();
    const double x = point->getX(), y = point->getY();
    const double magnitude = feature->GetFieldAsDouble(speed) * input.speedMultiplier();
    const double angle = feature->GetFieldAsDouble(direction);
    const double maxAngle = input.directionUnits == "degrees" ? 360 : 2 * std::acos(-1.0);
    if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(magnitude) || magnitude < 0
        || magnitude > std::numeric_limits<float>::max()
        || !std::isfinite(angle) || angle < 0 || angle > maxAngle)
      fail("invalid coordinate, speed or direction (expected 0..360 degrees or 0..2pi radians)");
    const double radians = angle * (input.directionUnits == "degrees" ? std::acos(-1.0) / 180 : 1);
    // Meteorological FROM direction relative to the projected grid north.
    points.push_back({x, y, -magnitude * std::sin(radians), -magnitude * std::cos(radians)});
    xvalues.push_back(x);
    yvalues.push_back(y);
  }
  xs = axis(std::move(xvalues));
  ys = axis(std::move(yvalues));
  if (xs.size() > points.size() / ys.size() || xs.size() * ys.size() != points.size())
    fail("incomplete rectangular grid (missing or extra points)");
  wind.resize(points.size());
  std::vector<bool> filled(points.size(), false);
  for (const auto &p : points) {
    const size_t index = indexOf(xs, p[0]) + xs.size() * indexOf(ys, p[1]);
    if (filled[index]) fail("duplicate grid point");
    filled[index] = true;
    wind[index] = {p[2], p[3]};
  }
}

void WindShapefileGrid::validateCoverage(double xmin, double ymin, double xmax, double ymax) const
{
  if (!std::isfinite(xmin) || !std::isfinite(ymin) || !std::isfinite(xmax) || !std::isfinite(ymax)
      || xmax < xmin || ymax < ymin || xmin < xs.front() - tolerance || xmax > xs.back() + tolerance
      || ymin < ys.front() - tolerance || ymax > ys.back() + tolerance) {
    std::ostringstream message;
    message.precision(15);
    message << "domain outside wind point coverage; requested [" << xmin << "," << ymin << "]..["
            << xmax << "," << ymax << "], available [" << xs.front() << "," << ys.front()
            << "]..[" << xs.back() << "," << ys.back() << "]";
    fail(message.str());
  }
}

std::array<double, 2> WindShapefileGrid::sample(double x, double y) const
{
  validateCoverage(x, y, x, y);
  x = std::clamp(x, xs.front(), xs.back());
  y = std::clamp(y, ys.front(), ys.back());
  const size_t i = std::min(static_cast<size_t>(std::max<std::ptrdiff_t>(0,
      std::upper_bound(xs.begin(), xs.end(), x) - xs.begin() - 1)), xs.size() - 2);
  const size_t j = std::min(static_cast<size_t>(std::max<std::ptrdiff_t>(0,
      std::upper_bound(ys.begin(), ys.end(), y) - ys.begin() - 1)), ys.size() - 2);
  const double tx = (x - xs[i]) / (xs[i+1] - xs[i]);
  const double ty = (y - ys[j]) / (ys[j+1] - ys[j]);
  std::array<double, 2> result{};
  for (int c = 0; c < 2; ++c)
    result[c] = (1-ty) * ((1-tx)*wind[i + j*nx()][c] + tx*wind[i+1 + j*nx()][c])
                  + ty * ((1-tx)*wind[i + (j+1)*nx()][c] + tx*wind[i+1 + (j+1)*nx()][c]);
  return result;
}

void WindShapefileGrid::validateVectorCRS(const std::string &file) const
{
  Dataset ds(static_cast<GDALDataset *>(GDALOpenEx(file.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY,
                                                nullptr, nullptr, nullptr)), GDALClose);
  if (!ds) fail("cannot open building shapefile for coordinate validation");
  OGRSpatialReference source(projection.c_str());
  for (int i = 0; i < ds->GetLayerCount(); ++i) {
    const auto *srs = ds->GetLayer(i)->GetSpatialRef();
    if (!srs || !source.IsSame(srs))
      fail("building and wind shapefile CRS differ; reproject first");
  }
}

void WindShapefileGrid::validateDEM(const std::string &file, double ox, double oy,
                                   int originFlag, double offsetX, double offsetY,
                                   double simOriginX, double simOriginY) const
{
  Dataset ds(static_cast<GDALDataset *>(GDALOpen(file.c_str(), GA_ReadOnly)), GDALClose);
  if (!ds) fail("cannot open DEM for coordinate validation");
  OGRSpatialReference source(projection.c_str());
  const auto *demSrs = ds->GetSpatialRef();
  if (!demSrs || !source.IsSame(demSrs)) fail("DEM and wind shapefile CRS differ; reproject first");
  double gt[6];
  if (ds->GetGeoTransform(gt) != CE_None || gt[1] <= 0 || gt[5] >= 0 || gt[2] != 0 || gt[4] != 0)
    fail("DEM must be north-up with a valid geotransform");
  const double demX = gt[0], demY = gt[3] + ds->GetRasterYSize() * gt[5];
  const double dx = originFlag == 0 ? offsetX : simOriginX - demX;
  const double dy = originFlag == 0 ? offsetY : simOriginY - demY;
  // QES's DEM reader crops at integer raster indices. Reject fractional-pixel
  // origins instead of allowing the wind and terrain to shift relative to each other.
  if (std::abs(ox - (demX + dx)) > tolerance || std::abs(oy - (demY + dy)) > tolerance
      || std::abs(dx / gt[1] - std::round(dx / gt[1])) > tolerance
      || std::abs(dy / -gt[5] - std::round(dy / -gt[5])) > tolerance)
    fail("originX/originY must match the DEM crop origin on DEM pixel boundaries");
}
