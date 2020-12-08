#include <iostream>
#include "WoS.h"
#include "WoS2d.h"
#include "sdf_bc.h"
#include <igl/Timer.h>
#include <igl/png/writePNG.h>
#include <igl/opengl/glfw/Viewer.h>
int main(int argc, char **argv) {
  int wpp = 16;
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::VectorXd B;

  igl::read_triangle_mesh((argc > 1 ? argv[1] : "../data/knot.obj"), V, F);
  auto solf = [](Eigen::Vector3d v) -> double { return 1 / (Eigen::Vector3d(1.0, 1.0, 1.0) - v).norm(); };
  auto solf_grad = [](Eigen::Vector3d v) -> Eigen::Vector3d {
    return (Eigen::Vector3d(1.0, 1.0, 1.0) - v) / pow((Eigen::Vector3d(1.0, 1.0, 1.0) - v).norm(), 3);
  };
  B.resize(V.rows());
  for (int i = 0; i < V.rows(); i++) {
    B[i] = solf(V.row(i));
  }
  int nquery = 1;
  Eigen::MatrixXd P(nquery, 3), U_grad(nquery, 3), sol_grad(nquery, 3);
  Eigen::VectorXd U(nquery), sol(nquery);
  P << 0.4, 0.2, 0.1;
  igl::Timer timer;
  timer.start();
  walk_on_spheres_mesh(
    V, F, B, P, [](Eigen::Vector3d p) { return 0.0; }, 128000, U, U_grad);
  timer.stop();
  std::cout << "took " << timer.getElapsedTimeInSec() << "s" << std::endl;
  std::cout << "U:" << U << std::endl;
  std::cout << "U_grad:" << U_grad << std::endl;
  std::cout << "reference: U:" << solf(P.row(0)) << "\nU_grad: " << solf_grad(P.row(0)) << std::endl;
}