#include "WindShapefileInput.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace {
double profileTerm(double z, double z0, double inverseL)
{
  const double pi = std::acos(-1.0);
  double psi;
  if (z * inverseL >= 0) {
    psi = 4.7 * z * inverseL;
  } else {
    const double x = std::pow(1 - 15 * z * inverseL, 0.25);
    psi = -2 * std::log((1 + x) / 2) - std::log((1 + x*x) / 2)
          + 2 * std::atan(x) - pi / 2;
  }
  return std::log1p(z / z0) + psi;
}
}

void WindShapefileInput::parseValues()
{
  // Reject typos, repeated options, and unsupported multi-time configurations.
  const std::vector<std::string> allowed = {"file", "speedField", "directionField",
    "speedUnits", "originX", "originY", "height", "siteZ0", "reciprocal",
    "boundaryLayerFlag", "timeStamp", "directionUnits"};
  for (const auto &entry : tree) {
    if (entry.first == "<xmlcomment>") continue;
    if (std::find(allowed.begin(), allowed.end(), entry.first) == allowed.end()
        || tree.count(entry.first) != 1)
      throw ParseException("windShapefile: unknown or repeated option " + entry.first);
  }
  parsePrimitive(true, file, "file");
  parsePrimitive(false, speedField, "speedField");
  parsePrimitive(false, directionField, "directionField");
  parsePrimitive(true, speedUnits, "speedUnits");
  parsePrimitive(false, directionUnits, "directionUnits");
  parsePrimitive(true, originX, "originX");
  parsePrimitive(true, originY, "originY");
  parsePrimitive(true, height, "height");
  parsePrimitive(true, siteZ0, "siteZ0");
  parsePrimitive(true, reciprocal, "reciprocal");
  parsePrimitive(false, boundaryLayerFlag, "boundaryLayerFlag");
  parsePrimitive(true, timeStamp, "timeStamp");
  validate();
  QEStime checkedTime(timeStamp);
  file = QESfs::get_absolute_path(file);
}

double WindShapefileInput::speedMultiplier() const
{
  if (speedUnits == "mps" || speedUnits == "m/s") return 1;
  if (speedUnits == "mph") return 0.44704;
  if (speedUnits == "kph" || speedUnits == "km/h") return 1.0 / 3.6;
  if (speedUnits == "kts" || speedUnits == "knots") return 1852.0 / 3600;
  throw std::runtime_error("windShapefile: speedUnits must be mps, m/s, mph, kph, km/h, kts or knots");
}

void WindShapefileInput::validate() const
{
  if (file.empty() || speedField.empty() || directionField.empty() || timeStamp.empty())
    throw std::runtime_error("windShapefile: file, fields and timeStamp must not be empty");
  speedMultiplier();
  if (directionUnits != "degrees" && directionUnits != "radians")
    throw std::runtime_error("windShapefile: directionUnits must be degrees or radians");
  if (!std::isfinite(originX) || !std::isfinite(originY)
      || !std::isfinite(height) || height <= 0 || !std::isfinite(siteZ0) || siteZ0 <= 0
      || !std::isfinite(reciprocal))
    throw std::runtime_error("windShapefile: finite origin/stability and positive height/siteZ0 required");
  if (boundaryLayerFlag != 0 && boundaryLayerFlag != 1)
    throw std::runtime_error("windShapefile: only uniform (0) and log (1) profiles supported");
  const double reference = profileTerm(height, siteZ0, reciprocal);
  if (boundaryLayerFlag == 1 && (!std::isfinite(reference) || reference <= 0))
    throw std::runtime_error("windShapefile: invalid log-profile reference height/stability");
}

double WindShapefileInput::profileScale(double z) const
{
  if (!std::isfinite(z)) throw std::runtime_error("windShapefile: non-finite terrain height");
  if (z <= 0) return 0;
  if (boundaryLayerFlag == 0) return 1;
  const double reference = profileTerm(height, siteZ0, reciprocal);
  // Same free-convection limit as WindProfilerSensorType, independent of speed.
  if (0.4 / reference >= 2) return 1;
  const double scale = profileTerm(z, siteZ0, reciprocal) / reference;
  if (!std::isfinite(scale) || scale < 0)
    throw std::runtime_error("windShapefile: invalid vertical profile; check height, siteZ0 and reciprocal");
  return scale;
}
