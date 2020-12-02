#pragma once
// Solve ∆u = 0 over space at given poinst P subject to B on the given boundary
// mesh (V,F)
//
// Inputs:
//   V  #V by 3 list of surface mesh vertex positions
//   F  #F by 3 list of triangles 
//   B  #V by 1 list of Dirichlet boundary conditions
//   P  #P by 3 list of query positions
// Outputs:
//   U  #P by 1 list of values at query positions
#include <Eigen/Core>
#include <functional>
#include <igl/embree/EmbreeIntersector.h>
void walk_on_spheres(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXd & B,
  const Eigen::MatrixXd & P,
  Eigen::VectorXd & U);


// Solve ∆u = f over space at given poinst P subject to B on the given boundary
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

void walk_on_spheres(
  const Eigen::MatrixXd &V,
  const Eigen::MatrixXi &F,
  const Eigen::VectorXd &B,
  const Eigen::MatrixXd &P,
  const std::function<double(const Eigen::Vector3d)> &f,
  Eigen::VectorXd &U,
  Eigen::MatrixXd &U_grad);


inline void walk_on_spheres(
  const Eigen::MatrixXd &V,
  const Eigen::MatrixXi &F,
  const Eigen::VectorXd &B,
  const Eigen::MatrixXd &P,
  const std::function<double(const Eigen::Vector3d)> &f,
  Eigen::VectorXd &U){
    Eigen::MatrixXd U_grad;
  walk_on_spheres(V,F,B,P,f, U, U_grad);
}

void walk_on_spheres(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                     const Eigen::VectorXd &B, const Eigen::MatrixXd &P,
                     const std::function<double(const Eigen::Vector3d)> &f,
                     Eigen::VectorXd &U,
                     Eigen::MatrixXd &U_grad);