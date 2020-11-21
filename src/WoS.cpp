#include <igl/AABB.h>
#include <igl/PI.h>
#include <igl/barycentric_coordinates.h>
#include <igl/barycentric_interpolation.h>
#include <igl/point_mesh_squared_distance.h>
#include <limits>
#define _USE_MATH_DEFINES
#include <cmath>
#include <green.h>
// returns a random value in the range [rMin,rMax]
// Copied from
// http://www.cs.cmu.edu/~kmcrane/Projects/MonteCarloGeometryProcessing/WoSLaplace2D.cpp.html
double random(double rMin, double rMax) {
  const double rRandMax = 1. / (double)RAND_MAX;
  double u = rRandMax * (double)rand();
  return u * (rMax - rMin) + rMin;
}

Eigen::VectorXd random(int N, double rMin, double rMax) {
  Eigen::VectorXd result(N);
  for (int i = 0; i < N; i++) {
    result[i] = random(rMin, rMax);
  }
  return result;
}

Eigen::Vector3d uniform_sphere_sampling() {
  double theta = random(0, igl::PI);
  double phi = random(0, 2 * igl::PI);

  Eigen::Vector3d new_direction;
  new_direction[0] = sin(theta) * cos(phi);
  new_direction[1] = sin(theta) * sin(phi);
  new_direction[2] = cos(theta);
  return new_direction;
}
Eigen::Vector3d uniform_ball_sampling() {
  return uniform_sphere_sampling() * sqrt(random(0, 1));
}
double sphere_volume(double R) { return 4.0 / 3.0 * igl::PI * R * R * R; }
// single point estimator for ∆u = f
double walk_on_spheres_single_point(
    const igl::AABB<Eigen::MatrixXd, 3> &aabb, const Eigen::MatrixXd &V,
    const Eigen::MatrixXi &F, const Eigen::VectorXd &B,
    const std::function<double(const Eigen::Vector3d)> &f,
    const Eigen::Vector3d &P) {
  const double eps = 0.01;
  const int nWalks = 128;
  const int maxSteps = 32;

  double sum = 0;
  int closest_face;
  for (int j = 0; j < nWalks; j++) { // j is a dummy var
    Eigen::Vector3d x(P);            // tmp
    int steps = 0;
    double R = 10000000.;
    double u = 0;
    do {

      Eigen::RowVector3d _closest;
      double sqrR = aabb.squared_distance(V, F, x, closest_face, _closest);

      Eigen::Vector3d new_direction = uniform_sphere_sampling();
      R = sqrt(sqrR);
      Eigen::Vector3d x_k1 = x + new_direction * R;
      Eigen::Vector3d y = uniform_ball_sampling() * R;
      u += R * f(x) * lapg3d(x, y, R) * sphere_volume(R);
      x = x_k1;
      steps++;
    } while (R > eps && steps < maxSteps);
    if (R < eps) { // on boundary
      std::array<Eigen::RowVector3d, 3> triangle;
      triangle[0] = V.row(F(closest_face, 0));
      triangle[1] = V.row(F(closest_face, 1));
      triangle[2] = V.row(F(closest_face, 2));
      Eigen::RowVector3d bc;
      igl::barycentric_coordinates(Eigen::RowVector3d(x), triangle[0], triangle[1], triangle[2],
                                   bc);
      Eigen::Vector3d gx;
      gx[0] = B[F(closest_face, 0)];
      gx[1] = B[F(closest_face, 1)];
      gx[2] = B[F(closest_face, 2)];
      u += bc * gx;
    }
    sum += u;
    // dunno how author handles it, but this version will work as well.
  }
  return sum / nWalks;
}
void walk_on_spheres(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                     const Eigen::VectorXd &B, const Eigen::MatrixXd &P,
                     const std::function<double(const Eigen::Vector3d)> &f,
                     Eigen::VectorXd &U) {
  igl::AABB<Eigen::MatrixXd, 3> aabb;
  aabb.init(V, F);
  U.resize(P.rows());

  for (int i = 0; i < P.rows(); i++) {
    U[i] = walk_on_spheres_single_point(aabb, V, F, B, f, P.row(i));
  }
}
void walk_on_spheres(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                     const Eigen::VectorXd &B, const Eigen::MatrixXd &P,
                     Eigen::VectorXd &U) {
  return walk_on_spheres(
      V, F, B, P, [](const Eigen::Vector3d &x) { return 0.0; }, U);
}