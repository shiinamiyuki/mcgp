#include <igl/AABB.h>
#include <igl/PI.h>
#include <igl/barycentric_coordinates.h>
#include <igl/barycentric_interpolation.h>
#include <igl/point_mesh_squared_distance.h>
#include <limits>
#include <igl/parallel_for.h>
#include <random>
#define _USE_MATH_DEFINES
#include <cmath>
#include <green.h>
#include <WoS.h>
// returns a random value in the range [rMin,rMax]
// Copied from
// http://www.cs.cmu.edu/~kmcrane/Projects/MonteCarloGeometryProcessing/WoSLaplace2D.cpp.html
static thread_local std::random_device rd;
double random(double rMin, double rMax) {
  // const double rRandMax = 1. / (double)RAND_MAX;
  // double u = rRandMax * (double)rand();
  // return u * (rMax - rMin) + rMin;
  std::uniform_real_distribution<double> dist(rMin, rMax);
  return dist(rd);
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
    const igl::AABB<Eigen::MatrixXd, 3> & aabb, 
    igl::embree::EmbreeIntersector * ei, const Eigen::MatrixXd &V,
    const Eigen::MatrixXi &F, const Eigen::VectorXd &B,
    const std::function<double(const Eigen::Vector3d)> &f,
    const Eigen::Vector3d &P) {
  const double eps = 0.01;
  const int nWalks = 32;
  const int maxSteps = 16;

  double sum = 0;
 
  for (int j = 0; j < nWalks; j++) { // j is a dummy var
    Eigen::Vector3d x(P);            // tmp
     int closest_face;
    int steps = 0;
    double R = 10000000.;
    double u = 0;
    double R_last;
    double k = 1;
    for(int steps =0 ;; steps++) {

      // Eigen::RowVector3d _closest;
      // int closest_face;
      // double sqrR = aabb.squared_distance(V, F, x, closest_face, _closest);
      // R = sqrt(sqrR);
      // double distance;
      float distance;
      ei->distance(x.cast<float>(), distance, closest_face);
      R = distance;
      if(R < eps)break;
      if(steps < maxSteps){
        R_last = R;
      }
      Eigen::Vector3d new_direction = uniform_sphere_sampling();
      Eigen::Vector3d x_k1 = x + new_direction * R;
      Eigen::Vector3d y = uniform_ball_sampling() * R;
      u += k * R * f(x) * lapg3d(x, y, R) * sphere_volume(R);
      x = x_k1;
      if(steps >= maxSteps){
        auto continue_prob = std::fmin(1.0, R_last / R) * 0.95;
        if(random(0, 1) < continue_prob){
          k /= continue_prob;
        }else{
          // printf("exit at %d steps\n", steps);
          break;
        }
        R_last = std::min(R_last, R);
      }
    }
    if (R <= eps) { // on boundary
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
      u += k * bc * gx;
    }
    sum += u;
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
  Eigen::MatrixXf Vf = V.cast<float>();
  igl::embree::EmbreeIntersector ei;
  ei.init(Vf, F);
  igl::parallel_for(P.rows(), [&](int i){
    U[i] = walk_on_spheres_single_point(aabb, &ei, V, F, B, f, P.row(i));
  });
}
void walk_on_spheres(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                     const Eigen::VectorXd &B, const Eigen::MatrixXd &P,
                     Eigen::VectorXd &U) {
  return walk_on_spheres(
      V, F, B, P, [](const Eigen::Vector3d &x) { return 0.0; }, U);
}