#pragma once

#include <Eigen/Core>

/* Decomposite vector field X into X = grad u + curl A + Y
*/
void helm3d(
  const std::function<double(const Eigen::Vector3d)> sdf,
  const Eigen::Vector3d center,
  const double Rmax,
  const std::function<Eigen::Vector3d(Eigen::Vector3d)> X,
  const std::function<double(Eigen::Vector3d)> divX,
  const std::function<Eigen::Vector3d(Eigen::Vector3d)> curlX,
  const Eigen::MatrixXd &P,
  const int nWalks,
  Eigen::MatrixXd &gradu,
  Eigen::MatrixXd &curlA,
  Eigen::MatrixXd &Y) ;