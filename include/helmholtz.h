#pragma once

#include <Eigen/Core>

void helm3d(
  const std::function<double(const Eigen::Vector3d)> sdf,
  const std::function<Eigen::Vector3d(Eigen::Vector3d)> X,
  const Eigen::MatrixXd &P,
  Eigen::MatrixXd gradu,
  Eigen::MatrixXd curlA,
  Eigen::MatrixXd Y
);