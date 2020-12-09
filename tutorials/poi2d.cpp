#include "WoS2d.h"
#include "sdf_bc.h"
#include "igl/PI.h"
#include <igl/png/writePNG.h>
#include <iostream>

int main() {
  auto sdf = [](Eigen::Vector2d p) { return sqrt(pow(p[0],2) + pow(p[1],2)) - 1; };
  auto solf = [](Eigen::Vector2d p) { return cos(2*igl::PI*p[0])*sin(2*igl::PI*p[1]); };
  auto solgrad = [](Eigen::Vector2d p)->Eigen::Vector2d { 
    return Eigen::Vector2d(-2*igl::PI*sin(2*igl::PI*p[0])*sin(2*igl::PI*p[1]),
            2*igl::PI*cos(2*igl::PI*p[0])*cos(2*igl::PI*p[1])); 
  };
  auto f = [](Eigen::Vector2d p) { return 8*pow(igl::PI,2)*cos(2*igl::PI*p[0])*sin(2*igl::PI*p[1]); };
  auto sdfbc = sdf_bc(sdf, solf);
  // double sd, bcval;
  // Eigen::Vector2d p(2,2);
  // std::tie(sd,bcval) = sdfbc(p);
  // std::cout << sd << " " << bcval << std::endl;

  int w = 128, h = 128;
  int nquery = w * h;
  Eigen::MatrixXd P(nquery, 2), Ugrad(nquery, 2), sol_grad(nquery,2);
  Eigen::VectorXd U(nquery), sol(nquery);

  for(int x = 0; x < w; x++){
    for(int y = 0; y < h; y++){
      int i = x + y * w;
      double xx = 2*(double(x)/w)-1;
      double yy = 2*(double(y)/h)-1;
      Eigen::Vector2d xy(xx,yy);
      if (sdf(xy) < 0) {    
        P.row(i) = xy;
      } else {
        P.row(i) = Eigen::Vector2d::Zero();
      }
    }
  }


  walk_on_spheres2d(sdfbc,f,P,4096,U,Ugrad);


  for(int x = 0; x < w; x++){
    for(int y = 0; y < h; y++){
      int i = x + y * w;
      double xx = 2*(double(x)/w)-1;
      double yy = 2*(double(y)/h)-1;
      Eigen::Vector2d xy(xx,yy);
      if (sdf(xy) >= 0) {    
        U[i] = 0;
        sol[i] = 0;
      } else {
        sol[i] = solf(xy);
      }
    }
  }

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
  write_solution(U, "./poi2d.png");
  write_solution(sol, "./poi2d_groundtruth.png");

  return 0;
}
