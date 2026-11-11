#pragma once
#include "util/ParseInterface.h"
#include <memory>
#include <string>

class WindVTKGrid;
class WINDSInputData;

// Coordinates and velocities are in the declared projected CRS axes, in metres
// and metres/second. originZ is in the same vertical datum as the VTK and TIFF.
class WindVTKInput : public ParseInterface
{
public:
  std::string file, terrainFile, velocityUnits, timeStamp;
  int epsg = 0;
  double originX = 0, originY = 0, originZ = 0;
  double siteZ0 = 0.1, reciprocal = 0, terrainTolerance = 0.1;
  std::shared_ptr<WindVTKGrid> field;
  void parseValues() override;
  void prepare(WINDSInputData *);
  double localGround(double x, double y) const;
};
