#pragma once

#include <Eigen/Core>
#include <igl/centroid.h>

// find the centroid and the radius of the bounding sphere of the mesh(V, F)
inline void mesh_centroid_bound(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,Eigen::Vector3d &centroid, double &radius) {
  igl::centroid(V, F, centroid);
  double Rmax = 0.0;
  for (int i = 0; i < V.rows(); i++) {
    Rmax = std::max((V.row(i).transpose() - centroid).norm(), Rmax);
  }
  radius = Rmax;
}