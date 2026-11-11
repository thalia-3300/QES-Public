#pragma once
#include "WindShapefileInput.h"
#include <array>
#include <vector>

// No sensor-by-domain weight arrays: storage is O(number of source points).
class WindShapefileGrid
{
public:
  explicit WindShapefileGrid(const WindShapefileInput &input);
  std::array<double, 2> sample(double projectedX, double projectedY) const;
  void validateVectorCRS(const std::string &file) const;
  void validateCoverage(double xmin, double ymin, double xmax, double ymax) const;
  void validateDEM(const std::string &file, double expectedOriginX, double expectedOriginY,
                   int originFlag, double offsetX, double offsetY,
                   double simOriginX, double simOriginY) const;
  size_t nx() const { return xs.size(); }
  size_t ny() const { return ys.size(); }
  std::array<double, 4> bounds() const { return {xs.front(), ys.front(), xs.back(), ys.back()}; }

private:
  std::vector<double> xs, ys;
  std::vector<std::array<double, 2>> wind;
  std::string projection;
};
