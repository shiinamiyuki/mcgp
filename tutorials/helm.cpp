#include <iostream>
#include "WoS.h"
#include "WoS2d.h"
#include "sdf_bc.h"
#include <igl/png/writePNG.h>
#include <igl/opengl/glfw/Viewer.h>
#include <helmholtz.h>
int main(int argc, char **argv) {
  int wpp = 16;
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh((argc > 1 ? argv[1] : "../data/icosphere.obj"), V, F);
  wpp = argc > 2 ? std::stoi(argv[2]) : wpp;

  auto X = [](Eigen::Vector3d p) -> Eigen::Vector3d {
    return Eigen::Vector3d(2 * p[0] * p[1], -p[1] * p[1], exp(p[0] * p[1]));
  };
  auto divX = [](Eigen::Vector3d p) -> double { return 0.0; };
  // auto curlX = [](Eigen::Vector3d p) -> Eigen::Vector3d { return Eigen::Vector3d::Zero(); };
  auto curlX = [](Eigen::Vector3d p) -> Eigen::Vector3d {
    return Eigen::Vector3d(exp(p[0] * p[1]) * p[0], -exp(p[0] * p[1]) * p[1], -2 * p[0]);
  };
  auto sdf = [](Eigen::Vector3d p) { return p.norm() - 1.0; };
  Eigen::MatrixXd gradu, curlA, Y;
  Eigen::MatrixXd P1, P2;

  int cnt = 0;
  P1.resize(1000, 3);
  for (int x = 0; x < 10; x++) {
    for (int y = 0; y < 10; y++) {
      for (int z = 0; z < 10; z++) {
        Eigen::Vector3d p(x, y, z);
        p /= 10.0;
        p = 2.0 * p - Eigen::Vector3d::Ones();
        if (sdf(p) > 0.0) {
          continue;
        }
        P1.row(cnt) = p;
        cnt++;
      }
    }
  }
  P1.conservativeResize(cnt, Eigen::NoChange);
  P2.resizeLike(P1);

  helm3d(sdf, Eigen::Vector3d::Zero(), 2.0, X, divX, curlX, P1, 4096, gradu, curlA, Y);
  std::cout << curlA << std::endl;
  for (int i = 0; i < cnt; i++) {
    Eigen::Vector3d p = P1.row(i);
    P2.row(i) = p + Eigen::Vector3d(curlA.row(i));
  }
  Eigen::MatrixXd P;
  Eigen::MatrixXi PI;
  P.resize(P1.rows() + P2.rows(), 3);
  PI.resize(P1.rows(), 2);
  for (int i = 0; i < P1.rows(); i++) {
    P.row(i) = P1.row(i);
    P.row(i + P1.rows()) = P2.row(i);
    PI.row(i) = Eigen::Vector2i(i, i + P1.rows());
  }

  igl::opengl::glfw::Viewer viewer;
  viewer.data().clear();
  viewer.data().set_mesh(V, F);
  viewer.data().compute_normals();

  viewer.data().point_size = 4;
  viewer.data().set_edges(P, PI, Eigen::RowVector3d(1, 0, 0));
  viewer.launch();
}