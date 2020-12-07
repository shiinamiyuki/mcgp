#pragma once
#include <Eigen/Core>



// second degree mls
double moving_least_squares(const Eigen::VectorXd & f, const Eigen::MatrixXd & xs, const Eigen::Vector3d &p);