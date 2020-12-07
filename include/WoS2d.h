#pragma once

#include <Eigen/Core>

void walk_on_spheres(const std::function<std::pair<double,double>(const Eigen::Vector2d)> &sdf_bc,
                     const Eigen::MatrixXd &P,
                     const std::function<double(const Eigen::Vector2d)> &f,
                     int num_walks,
                     Eigen::VectorXd &U,
                     Eigen::MatrixXd &U_grad);