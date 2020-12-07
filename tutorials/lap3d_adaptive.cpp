#include <iostream>
#include "WoS.h"
#include "WoS2d.h"
#include "sdf_bc.h"
#include <igl/png/writePNG.h>
int main() {
  auto sdf = [](Eigen::Vector3d p) { return p.norm() - 1.0; };
  auto bc = [&](Eigen::Vector3d p) { return abs(p.x()) > 0.5 ? 1.0 : 0.0; };
  auto sdfbc = sdf_bc3d(sdf, bc);
  auto region = [](Eigen::Vector3d p) -> Eigen::Vector3d {
    return Eigen::Vector3d(p.x() * 2.0 - 1.0, p.y() * 2.0 - 1.0, 0.0);
  };
  WoSPointCloud point_cloud;
  constexpr int w = 256, h = 256;
  constexpr int wpp = 160;
  constexpr size_t total = w * h * wpp;
  walk_on_spheres3d_region(
    sdfbc, [](Eigen::Vector3d v) -> double { return 0.0; }, region, Eigen::Vector3d::Zero(), 2.0, total, point_cloud);

  int nquery = w * h;
  Eigen::MatrixXd P(nquery, 3), U_grad(nquery, 3), sol_grad(nquery, 3);
  Eigen::VectorXd U(nquery), sol(nquery);

  for (int x = 0; x < w; x++) {
    for (int y = 0; y < h; y++) {
      int i = x + y * w;
      P.row(i) = Eigen::Vector3d(2 * (double(x) / w) - 1, 2 * (double(y) / h) - 1, 0.0);
    }
  }
  wos_point_cloud_interpolate(point_cloud, P, U, U_grad);

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
  write_solution(U, "../images/adaptive.png");
  walk_on_spheres3d(
    sdfbc, P, [](Eigen::Vector3d p) { return 0.0; }, Eigen::Vector3d::Zero(), 2.0, wpp, U, U_grad);
  write_solution(U, "../images/uniform.png");
}