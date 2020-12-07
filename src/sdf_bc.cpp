#include "sdf_bc.h"
#include <iostream>
std::function<std::pair<double,double>(const Eigen::Vector2d)> sdf_bc(
  const std::function<double(const Eigen::Vector2d)> sdf,
  const std::function<double(const Eigen::Vector2d)> bc)
{
  double fdeps = std::sqrt(1.0e-16);
  Eigen::Vector2d epsvec1(0, fdeps);
  Eigen::Vector2d epsvec2(fdeps, 0);

  auto sdfgradx = [&](Eigen::Vector2d p) {
    return (sdf(p+epsvec1)-sdf(p))/fdeps;
  };
  auto sdfgrady = [&](Eigen::Vector2d p) {
    return (sdf(p+epsvec2)-sdf(p))/fdeps;
  };
  auto sdfgrad2 = [&](Eigen::Vector2d p) {
    return pow(sdfgradx(p),2) + pow(sdfgrady(p),2);
  };
  auto closest = [&](Eigen::Vector2d p) {
    return p - sdf(p) * Eigen::Vector2d(sdfgradx(p), sdfgrady(p)) / sdfgrad2(p);
  };

  // Eigen::Vector2d x(3,0);
  // std::cout << closest(x) << std::endl;

  return [&](Eigen::Vector2d p) {
    return std::make_pair(sdf(p), bc(closest(p)));
  };
}