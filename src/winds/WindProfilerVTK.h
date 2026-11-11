#pragma once
#include "WindProfilerType.h"

class WindProfilerVTK : public WindProfilerType
{
public:
  void interpolateWindProfile(const WINDSInputData *, WINDSGeneralData *) override;
  static void enforceMasks(WINDSGeneralData *);
  static void finishVelocity(WINDSGeneralData *);
  static void reportChange(const WINDSGeneralData *);
  static void reportFinal(const WINDSGeneralData *);
};
