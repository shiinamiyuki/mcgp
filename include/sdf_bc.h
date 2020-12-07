#pragma once

#include <functional>
#include <Eigen/Core>

std::function<std::pair<double,double>(const Eigen::Vector2d)> sdf_bc(
  const std::function<double(const Eigen::Vector2d)> sdf,
  const std::function<double(const Eigen::Vector2d)> bc);