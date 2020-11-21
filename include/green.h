#pragma once
#include <Eigen/Core>

double lapg3d(
  Eigen::Vector3d & x, 
  Eigen::Vector3d & y, 
  double R);

double lapg3d(
  double r,
  double R);