#include <iostream>
#include "WoS.h"
#include "WoS2d.h"
#include "sdf_bc.h"
#include <igl/PI.h>

int main() {
  // reference solution
  auto solf = [](Eigen::Vector3d v) { return cos(2.0 * igl::PI * v[0]) * sin(2.0 * igl::PI * v[1]); };
  auto solf_grad = [](Eigen::Vector3d v) -> Eigen::Vector3d {
    return Eigen::Vector3d(-2.0 * igl::PI * sin(2.0 * igl::PI * v[0]) * sin(2.0 * igl::PI * v[1]),
                           2.0 * igl::PI * cos(2.0 * igl::PI * v[0]) * cos(2.0 * igl::PI * v[1]), 0.0);
  };
  auto solf_lap = [](Eigen::Vector3d v) -> double {
    return 8.0 * igl::PI * igl::PI * cos(2.0 * igl::PI * v[0]) * sin(2.0 * igl::PI * v[1]);
  };
  auto sdf = [](Eigen::Vector3d p) { return p.norm() - 1.0; };
  auto bc = [&](Eigen::Vector3d p) { return solf(p); };
  auto sdfbc = sdf_bc3d(sdf, bc);
  int nquery = 1;
  Eigen::MatrixXd P(nquery, 3), U_grad(nquery, 3), sol_grad(nquery, 3);
  Eigen::VectorXd U(nquery), sol(nquery);
  P << 0.4, 0.2, 0.1;
  std::cout << "solving ..." << std::endl;
  walk_on_spheres3d(sdfbc, P, solf_lap, Eigen::Vector3d::Zero(), 2.0, 1024000, U, U_grad);
  std::cout << "U:" << U << std::endl;
  std::cout << "U_grad:" << U_grad << std::endl;
  std::cout << "reference: U:" << solf(P.row(0)) << "\nU_grad: " << solf_grad(P.row(0)) << std::endl;
}