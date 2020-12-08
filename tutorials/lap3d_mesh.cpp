#include <iostream>
#include "WoS.h"
#include "WoS2d.h"
#include "sdf_bc.h"
#include <igl/png/writePNG.h>
#include <igl/opengl/glfw/Viewer.h>
int main() {
  int wpp = 16;
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh((argc > 1 ? argv[1] : "../data/icosphere.obj"), V, F);
  wpp = argc > 2 ? std::stoi(argv[2]) : wpp;
}