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

  // auto X = [](Eigen::Vector3d p) -> Eigen::Vector3d {
  //   return Eigen::Vector3d(2 * p[0] * p[1], -p[1] * p[1], exp(p[0] * p[1]));
  // };
  auto X = [](Eigen::Vector3d p) -> Eigen::Vector3d {
    return Eigen::Vector3d(p[0]*p[0]*cos(p[1]), p[0]*p[1]*p[2], exp(p[0] * p[1]));
  };
  auto divX = [](Eigen::Vector3d p) -> double { return p[0]*p[2]+2*p[0]*cos(p[1]); };
  // auto curlX = [](Eigen::Vector3d p) -> Eigen::Vector3d { return Eigen::Vector3d::Zero(); };
  // auto curlX = [](Eigen::Vector3d p) -> Eigen::Vector3d {
  //   return Eigen::Vector3d(exp(p[0] * p[1]) * p[0], -exp(p[0] * p[1]) * p[1], -2 * p[0]);
  // };
  auto curlX = [](Eigen::Vector3d p) -> Eigen::Vector3d {
    return Eigen::Vector3d(exp(p[0]*p[1])-p[0]*p[1], -exp(p[0] * p[1]) * p[1], p[0]*p[0]*sin(p[1]));
  };
  auto sdf = [](Eigen::Vector3d p) { return p.norm() - 1.0; };
  Eigen::MatrixXd gradu, curlA, Y;
  Eigen::MatrixXd P1, P2, P3, P4;

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
  std::cout << "solving ..." << std::endl;
  P3.resizeLike(P1);
  // P4.resizeLike(P1);

  helm3d(sdf, Eigen::Vector3d::Zero(), 2.0, X, divX, curlX, P1, 4096, gradu, curlA, Y);
  std::cout << curlA << std::endl;
  for (int i = 0; i < cnt; i++) {
    Eigen::Vector3d p = P1.row(i);
    P2.row(i) = p + Eigen::Vector3d(gradu.row(i));
    P3.row(i) = p + Eigen::Vector3d(curlA.row(i));
    // P4.row(i) = p + X(p);
  }
  Eigen::MatrixXd P, PP, PPP; // P for gradu, PP for curlA, PPP for X
  Eigen::MatrixXi PI;
  P.resize(P1.rows() + P2.rows(), 3);
  PP.resize(P1.rows() + P3.rows(), 3);
  PI.resize(P1.rows(), 2);
  for (int i = 0; i < P1.rows(); i++) {
    P.row(i) = P1.row(i);
    P.row(i + P1.rows()) = P2.row(i);
    PP.row(i) = P1.row(i);
    PP.row(i + P1.rows()) = P3.row(i);
    // PPP.row(i) = P1.row(i);
    // PPP.row(i + P1.rows()) = P4.row(i);
    PI.row(i) = Eigen::Vector2i(i, i + P1.rows());
  }

  igl::opengl::glfw::Viewer viewer;
  viewer.data().clear();
  viewer.data().set_mesh(V, F);
  viewer.data().compute_normals();

  viewer.data().point_size = 4;
  viewer.data().set_edges(P, PI, Eigen::RowVector3d(1, 1, 0));
  viewer.data().show_faces = false;
  viewer.data().show_lines = false;

  int plot_ugrad = true;
  const auto & update = [&]()
  {
    if(plot_ugrad)
    {
      viewer.data().set_edges(P, PI, Eigen::RowVector3d(1, 1, 0));
    }else
    {
      viewer.data().set_edges(PP, PI, Eigen::RowVector3d(1, 0.47, 0.45));
    }
  };
  viewer.callback_key_pressed = 
    [&](igl::opengl::glfw::Viewer &, unsigned int key, int)
  {
    switch(key)
    {
      case ' ':
        plot_ugrad ^= 1;
        break;
      default:
        return false;
    }
    update();
    return true;
  };
  viewer.launch();
}