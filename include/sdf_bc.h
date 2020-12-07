#pragma once

#include <functional>
#include <Eigen/Core>

std::function<std::pair<double,double>(const Eigen::Vector2d)> sdf_bc(
  const std::function<double(const Eigen::Vector2d)> sdf,
  const std::function<double(const Eigen::Vector2d)> bc);

std::function<std::pair<double,double>(const Eigen::Vector3d)> sdf_bc3d(
  const std::function<double(const Eigen::Vector3d)> sdf,
  const std::function<double(const Eigen::Vector3d)> bc);