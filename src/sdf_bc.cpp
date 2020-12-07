#include "sdf_bc.h"
#include <iostream>
std::function<std::pair<double, double>(const Eigen::Vector2d)>
sdf_bc(const std::function<double(const Eigen::Vector2d)> sdf, const std::function<double(const Eigen::Vector2d)> bc) {
  double fdeps = std::sqrt(1.0e-16);
  Eigen::Vector2d epsvec1(fdeps, 0);
  Eigen::Vector2d epsvec2(0, fdeps);

  auto sdfgradx = [=](Eigen::Vector2d p) { return (sdf(p + epsvec1) - sdf(p)) / fdeps; };
  auto sdfgrady = [=](Eigen::Vector2d p) { return (sdf(p + epsvec2) - sdf(p)) / fdeps; };
  auto sdfgrad2 = [=](Eigen::Vector2d p) { return pow(sdfgradx(p), 2) + pow(sdfgrady(p), 2); };
  auto closest = [=](Eigen::Vector2d p) -> Eigen::Vector2d {
    return p - sdf(p) * Eigen::Vector2d(sdfgradx(p), sdfgrady(p)) / sdfgrad2(p);
  };

  Eigen::Vector2d x(2, 2);
  Eigen::Vector2d y(sqrt(2) / 2, sqrt(2) / 2);
  std::cout << closest(x) << std::endl;
  std::cout << "bc(closest(x))=" << bc(closest(x)) << std::endl;
  std::cout << "bc(y)=" << bc(y) << std::endl;

  return [=](Eigen::Vector2d p) {
    // std::cout << p << std::endl;
    // std::cout << "bc(closest(p))=" << bc(closest(p)) << std::endl;
    return std::make_pair(sdf(p), bc(closest(p)));
  };
}

std::function<std::pair<double, double>(const Eigen::Vector3d)>
sdf_bc3d(const std::function<double(const Eigen::Vector3d)> sdf, const std::function<double(const Eigen::Vector3d)> bc) {
  double fdeps = std::sqrt(1.0e-16);
  Eigen::Vector3d epsvec1(fdeps, 0, 0);
  Eigen::Vector3d epsvec2(0, fdeps, 0);
  Eigen::Vector3d epsvec3(0, 0, fdeps);

  auto sdfgradx = [=](Eigen::Vector3d p) { return (sdf(p + epsvec1) - sdf(p)) / fdeps; };
  auto sdfgrady = [=](Eigen::Vector3d p) { return (sdf(p + epsvec2) - sdf(p)) / fdeps; };
  auto sdfgradz = [=](Eigen::Vector3d p) { return (sdf(p + epsvec3) - sdf(p)) / fdeps; };
  auto sdfgrad2 = [=](Eigen::Vector3d p) { return p.squaredNorm(); };
  auto closest = [=](Eigen::Vector3d p) -> Eigen::Vector3d {
    auto grad = Eigen::Vector3d(sdfgradx(p), sdfgrady(p),sdfgradz(p));
    return p - sdf(p) * grad.normalized();
  };

  return [=](Eigen::Vector3d p) {
    return std::make_pair(sdf(p), bc(closest(p)));
  };
}