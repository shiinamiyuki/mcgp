#include <iostream>
#include "WoS.h"
#include "WoS2d.h"
#include "sdf_bc.h"
#include <igl/Timer.h>
#include <igl/png/writePNG.h>
#include <igl/opengl/glfw/Viewer.h>
int main(int argc, char **argv) {
  int wpp = 16;
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::VectorXd B;
  wpp = argc > 2 ? std::stoi(argv[2]) : wpp;

  int equation = 0;
  equation = argc > 3 ? std::stoi(argv[3]) : equation;
  igl::read_triangle_mesh((argc > 1 ? argv[1] : "../data/knot.obj"), V, F);
  for (int i = 0; i < V.rows(); i++) {
    double x = V(i, 0);
    double z = V(i, 2);
    V(i, 0) = z;
    V(i, 2) = -x;
  }
  Eigen::MatrixXd VP;
  Eigen::MatrixXi FP;
  // igl::read_triangle_mesh("../data/plane.obj", VP, FP);
  auto solf = [](Eigen::Vector3d v) -> double { return 1 / (Eigen::Vector3d(1, 1, 1) - v).norm(); };
  B.resize(V.rows());
  for (int i = 0; i < V.rows(); i++) {
    if (equation == 1) {
      auto x = V.row(i).x();
      double fract = std::fmod((x + 1) * 2, 1.0);
      B[i] = (fract > 0.5) ? 1.0 : 0.0;
    } else {
      B[i] = solf(V.row(i));
    }
  }
  int w = 256, h = 256;
  int nquery = w * h;
  Eigen::MatrixXd P(nquery, 3), U_grad(nquery, 3), sol_grad(nquery, 3);
  Eigen::VectorXd U(nquery), sol(nquery);

  for (int x = 0; x < w; x++) {
    for (int y = 0; y < h; y++) {
      int i = x + y * w;
      P.row(i) = Eigen::Vector3d(2 * (double(x) / w) - 1, 2 * (double(y) / h) - 1, 0.0);
    }
  }

  std::cout << "solving ..." << std::endl;
  igl::Timer timer;
  timer.start();
  walk_on_spheres_mesh(
    V, F, B, P, [](Eigen::Vector3d p) { return 0.0; }, wpp, U, U_grad);
  timer.stop();
  std::cout << "took " << timer.getElapsedTimeInSec() << "s" << std::endl;
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
  write_solution(U, "mesh_fancy.png");
  Eigen::VectorXd ref;
  ref.resize(P.rows());
  for (int i = 0; i < P.rows(); i++) {
    ref[i] = solf(P.row(i));
  }
  write_solution(ref, "mesh_fancy_ref.png");
  igl::opengl::glfw::Viewer viewer;
  const int xid = viewer.selected_data_index;
  viewer.append_mesh();
  const int yid = viewer.selected_data_index;
  viewer.data_list[yid].set_mesh(V, F);
  // viewer.data_list[yid].set_colors(Eigen::Vector3d::Ones());
  {
    Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> R, G, B, A;
    R.resize(w, h);
    G.resize(w, h);
    B.resize(w, h);
    A.resize(w, h);
    {
      for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y++) {
          int i = x + y * w;
          auto val = fmin(1.0, fmax(U[i], 0.0));
          R(x, y) = (unsigned char)(val * 255);
          G(x, y) = (unsigned char)(val * 255);
          B(x, y) = (unsigned char)(val * 255);
          A(x, y) = 255;
        }
      }
    }
    VP.resize(4, 3);
    FP.resize(2, 3);
    VP << -1, 1, 0, -1, -1, 0, 1, -1, 0, 1, 1, 0;
    FP << 0, 1, 2, 2, 3, 0;
    Eigen::MatrixXd uv;
    uv.resize(4, 2);
    uv << 0, 1, 0, 0, 1, 0, 1, 1;
    // uv << 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0;
    // uv = Eigen::MatrixXd::Ones(4, 2) - uv;
    viewer.data_list[xid].set_mesh(VP, FP);
    viewer.data_list[xid].set_uv(uv);
    viewer.data_list[xid].set_colors(Eigen::Vector3d::Ones());
    viewer.data_list[xid].set_texture(R, G, B);
    viewer.data_list[xid].show_texture = 1;
    viewer.data_list[xid].double_sided = true;

    //  viewer.data_list[xid].
    // viewer.data_list[xid].updateGL( viewer.data_list[xid], false, viewer.data_list[xid].meshgl);
    // std::cout << VP << std::endl;
  }
  // viewer.data().show_texture = true;
  viewer.launch();
  return 0;
}