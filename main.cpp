#include <Eigen/Core>
#include <igl/read_triangle_mesh.h>
#include <igl/parula.h>
#include <iostream>
#include "WoS.h"
#include <igl/opengl/glfw/Viewer.h>

int main(int argc, char *argv[]){
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(
    (argc>1?argv[1]:"../data/icosphere.obj"),V,F);

  igl::opengl::glfw::Viewer viewer;
  const auto & update = [&]()
  {
    viewer.data().set_mesh(V,F);
    viewer.data().compute_normals();

  };

  int nV = V.rows();
  std::cout << nV << std::endl;
  std::cout << V.mean() << std::endl;
  Eigen::VectorXd B(nV);
  auto solf = [](Eigen::Vector3d v) { return (v[0]*v[0] - v[1]*v[1] + v[2]); }; // x^2 - y^2 + z

  for (int i = 0; i < nV; i++) {
    B[i] = solf(V.row(i));
  }
  std::cout << B.mean() << std::endl;
  int nquery = 3;
  Eigen::MatrixXd P(nquery, 3);
  Eigen::VectorXd U(nquery), sol(nquery);
  P << 0., 0., 0.,
       0.05,0.06,0.07,
       -0.1,-0.02,0.05;
  walk_on_spheres(V, F, B, P, U);

  for (int i = 0; i < nquery; i++) {
    sol[i] = solf(P.row(i));
  }


  std::cout << "approx: " << std::endl << U << std::endl;
  std::cout << "true sol: " << std::endl << sol << std::endl;


  update();
  viewer.launch();
  return 0;
}