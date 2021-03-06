#pragma once
#include <Eigen/Core>
#include <igl/PI.h>
double lapg3d(
  Eigen::Vector3d & x, 
  Eigen::Vector3d & y, 
  double R);

double lapg3d(
  double r,
  double R);


inline double harmonic_green(const Eigen::Vector3d  & x, const Eigen::Vector3d  & y, double R){
  double r =(x-y).norm();
  return 1/(2*igl::PI)*std::log(R / r);
}

Eigen::Vector3d lapdg3d(
  Eigen::Vector3d & x, 
  Eigen::Vector3d & y, 
  double R);

double lapg2d(
  double r,
  double R);

Eigen::Vector2d lapdg2d(
  Eigen::Vector2d & x, 
  Eigen::Vector2d & y, 
  double R);