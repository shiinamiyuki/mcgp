#pragma once

#include <functional>
#include <Eigen/Core>

/* Construct a pair<sdf, bc> function by given individual sdf, bc function
*/
std::function<std::pair<double,double>(const Eigen::Vector2d)> sdf_bc(
  const std::function<double(const Eigen::Vector2d)> sdf,
  const std::function<double(const Eigen::Vector2d)> bc);

std::function<std::pair<double,double>(const Eigen::Vector3d)> sdf_bc3d(
  const std::function<double(const Eigen::Vector3d)> sdf,
  const std::function<double(const Eigen::Vector3d)> bc);

// Construct a pair<sdf, bc> function from a mesh(V, F)
/*
Inputs:
//   V  #V by 3 list of surface mesh vertex positions
//   F  #F by 3 list of triangles
//   B  #V by 1 list of Dirichlet boundary conditions
*/
std::function<std::pair<double,double>(const Eigen::Vector3d)> sdf_bc_mesh(
  const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::VectorXd &B);