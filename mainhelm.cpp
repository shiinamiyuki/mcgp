#include <iostream>
#include "WoS.h"
#include "helmholtz.h"

int main() {


  auto sdf = [](Eigen::Vector3d p) { return p.norm() - 1.0; };
  Eigen::Vector3d center(0,0,0);
  double Rmax = 2;
  // auto X = [](Eigen::Vector3d p) -> Eigen::Vector3d { return Eigen::Vector3d(p[0],-2*p[1],p[2]); };
  auto X = [](Eigen::Vector3d p) -> Eigen::Vector3d { return Eigen::Vector3d(2*p[0]*p[1],-p[1]*p[1],exp(p[0]*p[1])); };
  auto divX = [](Eigen::Vector3d p) -> double { return 0.0; };
  // auto curlX = [](Eigen::Vector3d p) -> Eigen::Vector3d { return Eigen::Vector3d::Zero(); };
  auto curlX = [](Eigen::Vector3d p) -> Eigen::Vector3d { return Eigen::Vector3d(exp(p[0]*p[1])*p[0],-exp(p[0]*p[1])*p[1],-2*p[0]); };

  Eigen::MatrixXd P(1,3), gradu, curlA, Y;
  P << 0.1, 0.2, 0.3;
  helm3d(sdf,center,Rmax,X,divX,curlX,P,100000,gradu,curlA,Y);

  std::cout << gradu << std::endl;
  std::cout << curlA << std::endl;
  std::cout << Y << std::endl;
  std::cout << P << std::endl;
  return 0;
}