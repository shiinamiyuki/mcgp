#pragma once
#include <Eigen/Core>


template<int Degree=2>
double moving_least_squares(const Eigen::VectorXd & f, const Eigen::MatrixXd & xs, const Eigen::Vector3d &p);