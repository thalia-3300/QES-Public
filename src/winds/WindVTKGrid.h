#pragma once
#include <array>
#include <string>
#include <vector>

// Legacy binary STRUCTURED_GRID, separable increasing x/y, vertical columns.
// Trilinear isoparametric interpolation uses the physical height of each of
// the two bounding bilinear surfaces. No AGL profile or extrapolation.
class WindVTKGrid
{
public:
  explicit WindVTKGrid(const std::string &file);
  std::array<double, 3> sample(double x, double y, double z) const;
  double ground(double x, double y) const;
  double groundMaximum(double xmin, double ymin, double xmax, double ymax) const;
  double top(double x, double y) const;
  void validateTerrain(const std::string &file, int epsg, double tolerance,
                       double xmin, double ymin, double xmax, double ymax) const;
  void validateVectorCRS(const std::string &file, int epsg) const;
  int nx() const { return xs.size(); }
  int ny() const { return ys.size(); }
  int nz() const { return levels; }
  const std::vector<double> &x() const { return xs; }
  const std::vector<double> &y() const { return ys; }
private:
  struct Column { size_t ids[4]; double weights[4]; };
  Column column(double x, double y) const;
  double height(const Column &, int k) const;
  std::vector<double> xs, ys, zs;
  std::vector<std::array<double, 3>> velocity;
  int levels = 0;
};
