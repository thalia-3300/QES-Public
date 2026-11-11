#pragma once

#include "util/ParseInterface.h"
#include "util/QESFileSystemHandler.h"
#include "util/QEStime.h"
#include <string>

// One projected regular point grid, at one explicitly specified height/time.
class WindShapefileInput : public ParseInterface
{
public:
  std::string file, speedField = "speed", directionField = "dir", speedUnits;
  std::string timeStamp;
  std::string directionUnits = "degrees";
  double originX = 0, originY = 0; // projected coordinates of QES local (0,0)
  double height = 0, siteZ0 = 0.1, reciprocal = 0;
  int boundaryLayerFlag = 1; // 0: uniform; 1: QES stability-corrected log profile

  void parseValues() override;
  void validate() const;
  double speedMultiplier() const;
  double profileScale(double heightAboveGround) const;
};
