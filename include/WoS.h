#pragma once

#include <Eigen/Core>
#include <functional>
#include <igl/embree/EmbreeIntersector.h>
// Solve ∆u = -f over space at given poinst P subject to B on the given boundary
// mesh (V,F)
//
// Inputs:
//   V  #V by 3 list of surface mesh vertex positions
//   F  #F by 3 list of triangles 
//   B  #V by 1 list of Dirichlet boundary conditions
//   P  #P by 3 list of query positions
//   F  R^3 -> R function
// Outputs:
//   U  #P by 1 list of values at query positions
//   U_grad #P by 3 matrix of gradient of u

void walk_on_spheres(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                     const Eigen::VectorXd &B, const Eigen::MatrixXd &P,
                     const std::function<double(const Eigen::Vector3d)> &f,
                     int num_walks,
                     Eigen::VectorXd &U,
                     Eigen::MatrixXd &U_grad);


// Solve ∆u = -f over space at given poinst P subject to B on the boundary given by implicit surface
//
// Inputs:
//   sdf_bc R^3 -> pair<signed distance, boundary condition
//   P  #P by 3 list of query positions
//   F  R^3 -> R function
// Outputs:
//   U  #P by 1 list of values at query positions
//   U_grad #P by 3 matrix of gradient of u

void walk_on_spheres3d(const std::function<std::pair<double,double>(const Eigen::Vector3d)> &sdf_bc,
                     const Eigen::MatrixXd &P,
                     const std::function<double(const Eigen::Vector3d)> &f,
                     int num_walks,
                     Eigen::VectorXd &U,
                     Eigen::MatrixXd &U_grad);