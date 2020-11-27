#include <Eigen/Core>
#include <igl/read_triangle_mesh.h>

int main(int argc, char *argv[]){
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(
    (argc>1?argv[1]:"../data/sphere.obj"),V,F);
  
  int nV = V.rows();
  Eigen::VectorXd B(nV);
  auto solf = [](Eigen::Vector3d v) { return v[0]*v[0] - v[1]*v[1] + v[2]; }; // x^2 - y^2 + z

  for (int i = 0; i < nV; i++) {
    B[i] = solf(V.row(i));
  }
}