#include <iostream>
#include "WoS.h"
#include "WoS2d.h"
#include "sdf_bc.h"
#include <igl/png/writePNG.h>
#include <igl/opengl/glfw/Viewer.h>
int main(int argc, char **argv) {
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh((argc > 1 ? argv[1] : "../data/icosphere.obj"), V, F);
  auto sdf = [](Eigen::Vector3d p) { return p.norm() - 1.0; };
  auto bc = [&](Eigen::Vector3d p) {
    auto theta = atan2(p.x(), p.y());
    if (theta < 0.0) {
      theta += 2.0 * igl::PI;
    }
    int i = (theta / (2 * igl::PI) * 36);
    return i % 2 == 0.0 ? 1.0 : 0.0;
  };
  auto sdfbc = sdf_bc3d(sdf, bc);
  auto region = [](Eigen::Vector3d p) -> Eigen::Vector3d {
    return Eigen::Vector3d(p.x() * 2.0 - 1.0, p.y() * 2.0 - 1.0, 0.0);
  };
  WoSPointCloud point_cloud;
  constexpr int w = 256, h = 256;
  constexpr int wpp = 1024;
  constexpr size_t total = w * h * wpp;
  std::cout << "computing point cloud..." << std::endl;
  walk_on_spheres3d_region(
    sdfbc, [](Eigen::Vector3d v) -> double { return 0.0; }, region, Eigen::Vector3d::Zero(), 2.0, 1024, 4096 * 4, total,
    point_cloud);

  int nquery = w * h;
  Eigen::MatrixXd P(nquery, 3), U_grad(nquery, 3), sol_grad(nquery, 3);
  Eigen::VectorXd U(nquery), sol(nquery);

  for (int x = 0; x < w; x++) {
    for (int y = 0; y < h; y++) {
      int i = x + y * w;
      P.row(i) = Eigen::Vector3d(2 * (double(x) / w) - 1, 2 * (double(y) / h) - 1, 0.0);
    }
  }
  std::cout << "interpolating point cloud..." << std::endl;
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
    std::cout << "write to " << path << std::endl;
  };
  write_solution(U, "./adaptive.png");

  std::cout << "uniform sampling ..." << std::endl;

  WoSPointCloud cloud_uniform;
  walk_on_spheres3d_region(
    sdfbc, [](Eigen::Vector3d v) -> double { return 0.0; }, region, Eigen::Vector3d::Zero(), 2.0, total / (4096 * 4),
    4096 * 4, total, cloud_uniform);
  wos_point_cloud_interpolate(cloud_uniform, P, U, U_grad);
  write_solution(U, "./uniform.png");
  igl::opengl::glfw::Viewer viewer;
  viewer.data().clear();
  viewer.data().set_mesh(V, F);
  viewer.data().compute_normals();
  viewer.data().point_size = 4;
  viewer.data().set_points(point_cloud.P, Eigen::RowVector3d(0, 0, 0));
  viewer.callback_key_pressed = [&](igl::opengl::glfw::Viewer &, unsigned int key, int mod) {
    switch (key) {
    case 'U':
    case 'u': {
      viewer.data().set_points(cloud_uniform.P, Eigen::RowVector3d(0, 0, 0));
    } break;
    case 'A':
    case 'a': {
      viewer.data().set_points(point_cloud.P, Eigen::RowVector3d(0, 0, 0));
    } break;
    default:
      return false;
    }
    return true;
  };
  viewer.launch();
  return 0;
}