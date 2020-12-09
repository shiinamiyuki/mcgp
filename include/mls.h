#pragma once
#include <Eigen/Core>


/*
1, 2, 3 degree moving least squares
Inputs:
  xs        #xs by 3 list of sample points
  f         #xs list of scalar function value evaluated at xs
  p         the location where MLS interpolates f
*/
template<int Degree=2>
double moving_least_squares(const Eigen::VectorXd & f, const Eigen::MatrixXd & xs, const Eigen::Vector3d &p);