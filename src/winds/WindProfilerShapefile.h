#pragma once
#include "WindProfilerType.h"
#include "WindShapefileGrid.h"
#include "qes/Domain.h"

class WindProfilerShapefile : public WindProfilerType
{
public:
  WindProfilerShapefile(const WINDSInputData *, const qes::Domain &);
  void interpolateWindProfile(const WINDSInputData *, WINDSGeneralData *) override;
private:
  WindShapefileGrid grid;
};
