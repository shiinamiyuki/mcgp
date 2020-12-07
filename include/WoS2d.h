#pragma once

#include <Eigen/Core>

void walk_on_spheres2d(const std::function<std::pair<double,double>(const Eigen::Vector2d)> &sdf_bc,
                     const std::function<double(const Eigen::Vector2d)> &f,
                     const Eigen::MatrixXd &P,
                     int num_walks,
                     Eigen::VectorXd &U,
                     Eigen::MatrixXd &U_grad);