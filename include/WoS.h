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

void walk_on_spheres_mesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::VectorXd &B,
                          const Eigen::MatrixXd &P, const std::function<double(const Eigen::Vector3d)> &f,
                          int num_walks, Eigen::VectorXd &U, Eigen::MatrixXd &U_grad);

// Solve ∆u = -f over space at given poinst P subject to B on the boundary given by implicit surface
//
// Inputs:
//   sdf_bc         R^3 -> pair<signed distance, boundary condition
//   P              #P by 3 list of query positions
//   F              R^3 -> R function
//   center         the center of the domain
//   Rmax           anything walks beyond Rmax from the center will be terminated
//                  (ensures the exterior problem can be solved in reasonable time)
// Outputs:
//   U  #P by 1 list of values at query positions
//   U_grad #P by 3 matrix of gradient of u

void walk_on_spheres3d(const std::function<std::pair<double, double>(const Eigen::Vector3d)> &sdf_bc,
                       const Eigen::MatrixXd &P, const std::function<double(const Eigen::Vector3d)> &f,
                       const Eigen::Vector3d &center, double Rmax, int num_walks, Eigen::VectorXd &U,
                       Eigen::MatrixXd &U_grad);

// //  Solve ∆u = -f over space quad defined by p = center + s * u + t * v, u, v in [-1, 1]
// void walk_on_spheres3d_region(const std::function<std::pair<double, double>(const Eigen::Vector3d)> &sdf_bc,
//                               const std::function<double(const Eigen::Vector3d)> &f, const Eigen::Vector3d &center,
//                               const Eigen::Vector3d &s, const Eigen::Vector3d &t, const Eigen::Vector2i &resolution,
//                               int num_walks, Eigen::VectorXd &U, Eigen::MatrixXd &U_grad);

struct WoSPointCloud {
  Eigen::MatrixXd P;
  Eigen::VectorXd U;
  Eigen::MatrixXd U_grad;
  void resize(int N) {
    P.conservativeResize(N, 3);
    U.conservativeResize(N);
    U_grad.conservativeResize(N, 3);
  }
  int n_points() const { return P.rows(); }
};

//  Solve ∆u = -f over space parameterized by region
// region     [0,1)^3 -> R^3
void walk_on_spheres3d_region(const std::function<std::pair<double, double>(const Eigen::Vector3d)> &sdf_bc,
                              const std::function<double(const Eigen::Vector3d)> &f,
                              const std::function<Eigen::Vector3d(Eigen::Vector3d)> &region,
                              const Eigen::Vector3d &center, double Rmax, size_t walks_per_point, size_t total_walks,
                              WoSPointCloud &point_cloud);

void wos_point_cloud_interpolate(const WoSPointCloud &point_cloud, const Eigen::MatrixXd &P, Eigen::VectorXd &U,
                                 Eigen::MatrixXd &U_grad);
