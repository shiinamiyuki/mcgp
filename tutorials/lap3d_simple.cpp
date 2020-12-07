#include <iostream>
#include "WoS.h"
#include "WoS2d.h"
#include "sdf_bc.h"

int main() {

  auto solf = [](Eigen::Vector3d v) { return 1 / (Eigen::Vector3d(1.0, 1.0, 1.0) - v).norm(); };
  auto solf_grad = [](Eigen::Vector3d v) {
    return (Eigen::Vector3d(1.0, 1.0, 1.0) - v) / pow((Eigen::Vector3d(1.0, 1.0, 1.0) - v).norm(), 3);
  };

  auto sdf = [](Eigen::Vector3d p) { return p.norm() - 1.0; };
  auto bc = [&](Eigen::Vector3d p) { return solf(p); };
  auto sdfbc = sdf_bc3d(sdf, bc);
  int nquery = 1;
  Eigen::MatrixXd P(nquery, 3), U_grad(nquery, 3), sol_grad(nquery, 3);
  Eigen::VectorXd U(nquery), sol(nquery);
  P << 0.4, 0.2, 0.1;
  walk_on_spheres3d(
    sdfbc, P, [](Eigen::Vector3d p) { return 0.0; }, Eigen::Vector3d::Zero(), 2.0, 1024000, U, U_grad);
  std::cout << "U:" << U << std::endl;
  std::cout << "U_grad:" << U_grad << std::endl;
  std::cout << "reference: U:" << solf(P.row(0)) << "\nU_grad: " << solf_grad(P.row(0)) << std::endl;
}