#include "helmholtz.h"
#include "WoS.h"
#include "sdf_bc.h"
void helm3d(
  const std::function<double(const Eigen::Vector3d)> sdf,
  const Eigen::Vector3d center,
  const double Rmax,
  const std::function<Eigen::Vector3d(Eigen::Vector3d)> X,
  const std::function<double(Eigen::Vector3d)> divX,
  const std::function<Eigen::Vector3d(Eigen::Vector3d)> curlX,
  const Eigen::MatrixXd &P,
  const int nWalks,
  Eigen::MatrixXd &gradu,
  Eigen::MatrixXd &curlA,
  Eigen::MatrixXd &Y) 
{
  
  auto bc = [](Eigen::Vector3d p) -> double { return 0.0; };

  auto sdfbc = sdf_bc3d(sdf, bc);

  Eigen::VectorXd u, A0, A1, A2;
  Eigen::MatrixXd gradA0, gradA1, gradA2;

  gradu.resize(P.rows(), 3);
  walk_on_spheres3d(sdfbc,P,[=](Eigen::Vector3d p){ return -divX(p); },center,Rmax,nWalks,u,gradu);
  // std::cout << gradu << std::endl;

  auto curlX0 = [=](Eigen::Vector3d p) -> double { return -curlX(p)[0]; };
  auto curlX1 = [=](Eigen::Vector3d p) -> double { return -curlX(p)[1]; };
  auto curlX2 = [=](Eigen::Vector3d p) -> double { return -curlX(p)[2]; };


  walk_on_spheres3d(sdfbc,P,curlX0,center,Rmax,nWalks,A0,gradA0);
  walk_on_spheres3d(sdfbc,P,curlX1,center,Rmax,nWalks,A1,gradA1);
  walk_on_spheres3d(sdfbc,P,curlX2,center,Rmax,nWalks,A2,gradA2);

  curlA.resize(P.rows(), 3);
  Eigen::Vector3d dx(1, 0, 0);
  Eigen::Vector3d dy(0, 1, 0);
  Eigen::Vector3d dz(0, 0, 1);

  for (int i = 0; i < P.rows(); i++) {
    double A0dz = gradA0.row(i).dot(dz);
    double A0dy = gradA0.row(i).dot(dy);
    double A1dz = gradA1.row(i).dot(dz);
    double A1dx = gradA1.row(i).dot(dx);
    double A2dx = gradA2.row(i).dot(dx);
    double A2dy = gradA2.row(i).dot(dy);
    curlA.row(i) = Eigen::Vector3d(A2dy-A1dz, A0dz-A2dx, A1dx-A0dy);
  }


  Y.resize(P.rows(), 3);

  for (int i = 0; i < P.rows(); i++) {
    Y.row(i) = X(P.row(i)) - (Eigen::Vector3d)(gradu.row(i)) - (Eigen::Vector3d)(curlA.row(i));
  }



}