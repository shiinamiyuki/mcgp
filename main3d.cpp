#include <Eigen/Core>
#include <igl/read_triangle_mesh.h>
#include <igl/parula.h>
#include <igl/png/writePNG.h>
#include <igl/PI.h>
#include <iostream>
#include "WoS.h"
#include "WoS2d.h"
#include "sdf_bc.h"
#include <igl/opengl/glfw/Viewer.h>

int wpp = 16;

void test_lap3d() {
  auto solf = [](Eigen::Vector3d v) { return 1 / (Eigen::Vector3d(1.0, 1.0, 1.0) - v).norm(); };
  auto solf_grad = [](Eigen::Vector3d v) {
    return (Eigen::Vector3d(1.0, 1.0, 1.0) - v) / pow((Eigen::Vector3d(1.0, 1.0, 1.0) - v).norm(), 3);
  };
  //   auto solf = [](Eigen::Vector3d v) { return v.sum(); };
  //   auto solf_grad = [](Eigen::Vector3d v) { return Eigen::Vector3d::Ones(); };

  auto sdf = [](Eigen::Vector3d p) { return p.norm() - 1.0; };
  auto bc = [&](Eigen::Vector3d p) { return solf(p); };
  auto sdfbc = sdf_bc3d(sdf, bc);
  int nquery = 1;
  Eigen::MatrixXd P(nquery, 3), U_grad(nquery, 3), sol_grad(nquery, 3);
  Eigen::VectorXd U(nquery), sol(nquery);
  P << 0.4, 0.2, 0.1;
  walk_on_spheres3d(
    sdfbc, P, [](Eigen::Vector3d p) { return 0.0; }, wpp, U, U_grad);
  std::cout << "U:" << U << std::endl;
  std::cout << "U_grad:" << U_grad << std::endl;
  std::cout << "reference: U:" << solf(P.row(0)) << "\nU_grad: " << solf_grad(P.row(0)) << std::endl;
}

void test_poi3d() {
  // auto solf = [](Eigen::Vector3d v) { return sin(v[0]) * cos(v[1]); };
  // auto solf_grad = [](Eigen::Vector3d v) -> Eigen::Vector3d {
  //   return Eigen::Vector3d(cos(v[0]) * cos(v[1]), -sin(v[0]) * sin(v[1]), 0.0);
  // };
  // auto solf_lap = [](Eigen::Vector3d v) -> double { return -2.0 * sin(v[0]) * cos(v[1]); };
  // //   auto solf = [](Eigen::Vector3d v) { return v.sum(); };
  // //   auto solf_grad = [](Eigen::Vector3d v) { return Eigen::Vector3d::Ones(); };
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
  walk_on_spheres3d(sdfbc, P, solf_lap, wpp, U, U_grad);
  std::cout << "P:" << P << std::endl;
  std::cout << "U:" << U << std::endl;
  std::cout << "U_grad:" << U_grad << std::endl;
  std::cout << "reference: U:" << solf(P.row(0)) << "\nU_grad: " << solf_grad(P.row(0)) << std::endl;
}

int main(int argc, char *argv[]) {
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh((argc > 1 ? argv[1] : "../data/icosphere.obj"), V, F);
  wpp = argc > 2 ? std::stoi(argv[2]) : wpp;
  // test_lapint(V,F);
  // test_lapintgrad(V,F);
  // test_lap3d();
  // test_poi3d();
  igl::opengl::glfw::Viewer viewer;
  const auto &update = [&]() {
    viewer.data().set_mesh(V, F);
    viewer.data().compute_normals();
  };

  int nV = V.rows();
  std::cout << nV << std::endl;
  std::cout << V.mean() << std::endl;
  Eigen::VectorXd B(nV);
  // auto solf = [](Eigen::Vector3d v) { return (v[0]*v[0] - v[1]*v[1] + v[2] + sin(v[2])); }; // x^2 - y^2 + z
  // auto solf = [](Eigen::Vector3d v) { return sin(igl::PI * v[1]) * cos(igl::PI * v[0]);};
  // auto solf = [](Eigen::Vector3d v) { return 1 / (Eigen::Vector3d(1.0, 1.0, 1.0) - v).norm(); };
  // auto solf_grad = [](Eigen::Vector3d v) {
  //   return (Eigen::Vector3d(1.0, 1.0, 1.0) - v) / pow((Eigen::Vector3d(1.0, 1.0, 1.0) - v).norm(), 3);
  // };
  auto solf = [](Eigen::Vector3d v) { return cos(2.0 * igl::PI * v[0]) * sin(2.0 * igl::PI * v[1]); };
  auto solf_grad = [](Eigen::Vector3d v) -> Eigen::Vector3d {
    return Eigen::Vector3d(-2.0 * igl::PI * sin(2.0 * igl::PI * v[0]) * sin(2.0 * igl::PI * v[1]),
                           2.0 * igl::PI * cos(2.0 * igl::PI * v[0]) * cos(2.0 * igl::PI * v[1]), 0.0);
  };
  auto solf_lap = [](Eigen::Vector3d v) -> double {
    return 8.0 * igl::PI * igl::PI * cos(2.0 * igl::PI * v[0]) * sin(2.0 * igl::PI * v[1]);
  };
  for (int i = 0; i < nV; i++) {
    B[i] = solf(V.row(i));
  }
  std::cout << B.mean() << std::endl;

  int w = 64, h = 64;
  int nquery = w * h;
  Eigen::MatrixXd P(nquery, 3), U_grad(nquery, 3), sol_grad(nquery, 3);
  Eigen::VectorXd U(nquery), sol(nquery);
  // P << 0., 0., 0.;
  //     //  0.05,0.06,0.07,
  //     //  -0.1,-0.02,0.05;
  for (int x = 0; x < w; x++) {
    for (int y = 0; y < h; y++) {
      int i = x + y * w;
      P.row(i) = Eigen::Vector3d(2 * (double(x) / w) - 1, 2 * (double(y) / h) - 1, 0.0);
    }
  }
  // walk_on_spheres(V, F, B, P, [](Eigen::Vector3d v)->double {
  //   return -2 * igl::PI * igl::PI * sin(igl::PI * v[1]) * cos(igl::PI * v[0]);
  //  }, U);
  // walk_on_spheres(V, F, B, P, [](Eigen::Vector3d v)->double { return -sin(v[2]); }, U, U_grad);
  auto sdf = [](Eigen::Vector3d p) { return p.norm() - 1.0; };
  auto bc = [&](Eigen::Vector3d p) { return solf(p); };
  auto sdfbc = sdf_bc3d(sdf, bc);
  walk_on_spheres3d(sdfbc, P, solf_lap, wpp, U, U_grad);

  for (int i = 0; i < nquery; i++) {
    sol[i] = solf(P.row(i));
    Eigen::Vector3d tmpg = solf_grad(P.row(i));
    if (i == 0)
      std::cout << P.row(i) << ", " << tmpg << std::endl;
    if (!(std::isfinite(tmpg[0]) && std::isfinite(tmpg[1]) && std::isfinite(tmpg[2]))) {
      sol_grad.row(i) = U_grad.row(i);
    } else {
      sol_grad.row(i) = tmpg;
    }
  }
  std::cout << "sol grad:" << sol_grad.block(0, 0, 20, 3) << std::endl;
  std::cout << "Ugrad:" << U_grad.block(0, 0, 20, 3) << std::endl;

  std::cout << "grad error: " << (sol_grad - U_grad).mean() << std::endl;

  auto write_solution = [=](const Eigen::VectorXd &u, const char *path) {
    Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> R(w, h);
    Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> G(w, h);
    Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> B(w, h);
    Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> A(w, h);
    for (int x = 0; x < w; x++) {
      for (int y = 0; y < h; y++) {
        int i = x + y * w;
        auto val = fmin(1.0, fmax(u[i], 0.0));
        R(x, y) = (unsigned char)(val * 255);
        G(x, y) = (unsigned char)(val * 255);
        B(x, y) = (unsigned char)(val * 255);
        A(x, y) = 255;
      }
    }
    igl::png::writePNG(R, G, B, A, path);
  };
  // std::cout << "approx: " << std::endl << U << std::endl;
  // std::cout << "true sol: " << std::endl << sol << std::endl;
  std::cout << "MSE: " << ((U - sol).array() * (U - sol).array()).mean() << std::endl;
  write_solution(U, "../images/mcgp.png");
  write_solution(sol, "../images/groundtruth.png");

  update();
  viewer.launch();
  return 0;
}