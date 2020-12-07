#include <igl/AABB.h>
#include <igl/PI.h>
#include <igl/barycentric_coordinates.h>
#include <igl/barycentric_interpolation.h>
#include <igl/point_mesh_squared_distance.h>
#include <limits>
#include <igl/parallel_for.h>
#include <random>
#include <tuple>
#define _USE_MATH_DEFINES
#include <cmath>
#include <green.h>
#include <WoS.h>

// returns a random value in the range [rMin,rMax]
static thread_local std::random_device rd;
double random(double rMin, double rMax) {
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
  auto u = random(-1, 1);
  auto phi = random(0, 2 * igl::PI);
  return Eigen::Vector3d(cos(phi) * sqrt(1 - u * u), sin(phi) * sqrt(1 - u * u), u);
}
// Eigen::Vector3d uniform_ball_sampling() {
//   auto u = random(-1, 1);
//   auto phi = random(0, 2 * igl::PI);
//   auto r = pow(random(0, 1), 1.0 / 3.0);
//   return Eigen::Vector3d(r * cos(phi) * sqrt(1 - u * u), r * sin(phi) * sqrt(1 - u * u), r * u);
// }
double sphere_volume(double R) { return 4.0 / 3.0 * igl::PI * R * R * R; }

std::pair<Eigen::Vector3d, double> importance_sample_green3d(double R) {
  const double intG = R * R / 6.0;
  Eigen::Vector3d y = uniform_sphere_sampling();
  // double theta = acos(y.z());
  // sample r prop. to r * r* sin(theta)
  // r' = r / R
  //  r * r* sin(theta) = R^2 * r' * sin(theta)
  // cdf(r') = integrate R^2 * r' * sin(theta) [0, r'] / intCDF
  //         = R^2 *  sin(theta) * r'^3 / 3 / intCDF
  // R^2 *  sin(theta) * r'^3 / 3 / intCDF= u
  // cdf^{-1}(u) = (intCDF * 3 * u / sin(theta) /  R^2)^(-1/3)
  // const auto intCDF = R * R / 3 * sin(theta);
  // auto u = random(0, 1);
  // double r;
  // if (theta == 0) {
  //   r = u;
  // } else {
  //   r = pow(intCDF * 3 * u / sin(theta) / (R * R), -1.0 / 3.0);
  // }
  // r *= R;
  // int r^2*sin(theta)*G = sin(theta) / (4*pi * r * R) * (r*r*R/2.0 - r*r*r/3.0)
  // const auto F = [=](double r) { return sin(theta) / (4 * igl::PI * r * R) * (r * r * R / 2.0 - r * r * r / 3.0); };
  // const auto CDF = [=](double r) { return F(r) / F(R); };
  // const auto invCDF = [=](double u){
  // u = sin(theta) / (4 * igl::PI * r * R) * (r * r * R / 2.0 - r * r * r / 3.0) / intCDF
  // u* intCDF / sin(theta) = 1.0 / (4 * igl::PI * r * R) * (r * r * R / 2.0 - r * r * r / 3.0)
  // };
  const auto r = pow(random(0, 1), 1.0 / 3.0) * R;
  return std::pair<Eigen::Vector3d, double>(y * r, r);
}
double pdf_green3d(double r, double R) {
  const double intG = R * R / 6.0;
  return lapg3d(r, R) / intG;
}

// single point estimator for ∆u = f
std::pair<double, Eigen::Vector3d>
walk_on_sphere_single_point3d(const std::function<std::pair<double, double>(const Eigen::Vector3d)> &sdf_bc,
                              const std::function<double(const Eigen::Vector3d)> &f, const Eigen::Vector3d &P,
                              int num_walks) {
  const double eps = 0.001;
  const int nWalks = num_walks;
  const int maxSteps = 32;

  double sum = 0;
  Eigen::Vector3d sumgrad = Eigen::Vector3d::Zero();

  for (int j = 0; j < nWalks; j++) { // j is a dummy var
    Eigen::Vector3d x(P);            // tmp
    int steps = 0;
    double R = 10000000.;
    double u = 0;
    double R_last;
    double k = 1;

    Eigen::Vector3d grad = Eigen::Vector3d::Zero();

    Eigen::Vector3d first_direction;
    double first_R;
    for (int steps = 0;; steps++) {
      double sd, bc;
      std::tie(sd, bc) = sdf_bc(x);
      double R = std::abs(sd);
      if (R < eps) {
        u += bc;
        if (j != 0) {
          double oldsum = sum;
          Eigen::Vector3d oldgrad = sumgrad / j;
          sum += u - oldgrad.dot(first_R * first_direction);
          sumgrad += (u - oldsum / j) * first_direction * 3 / first_R + grad;
        } else {
          sum += u;
          sumgrad += u * first_direction * 3 / first_R + grad;
        }
        break;
      }
      if (steps < maxSteps) {
        R_last = R;
      }

      Eigen::Vector3d new_direction = uniform_sphere_sampling();

      Eigen::Vector3d x_k1 = x + new_direction * R;

      {
        double r = pow(random(0, 1), 1.0 / 3.0) * R;
        // Eigen::Vector3d y;
        Eigen::Vector3d y = x + uniform_sphere_sampling() * r;
        u += k * f(y) * lapg3d(r, R) * sphere_volume(R);
        // std::tie(y, r) = importance_sample_green3d(R);
        // y = x + y;
        // u += k * f(y) * lapg3d(r, R) / pdf_green3d(r, R);
        // printf("%lf\n", lapg3d(r, R) / pdf_green3d(r, R));
      }

      if (steps == 0) {
        first_direction =
          new_direction
            .normalized(); // the normalize is unnecessary since the direction is always normalized.; yes, it is.
        first_R = R;
        double r = pow(random(0, 1), 1.0 / 3.0) * R;
        Eigen::Vector3d y = x + uniform_sphere_sampling() * r;
        grad = sphere_volume(R) * f(y) * lapdg3d(x, y, R);
      }

      x = x_k1;

      if (steps >= maxSteps) {
        auto continue_prob = std::fmin(1.0, R_last / R) * 0.95;
        if (std::isfinite(continue_prob) && random(0, 1) < continue_prob) {
          k /= continue_prob;
        } else {
          break;
        }
        R_last = std::min(R_last, R);
      }
    }
  }
  return std::make_pair(sum / nWalks, sumgrad / nWalks);
}

void walk_on_spheres3d(const std::function<std::pair<double, double>(const Eigen::Vector3d)> &sdf_bc,
                       const Eigen::MatrixXd &P, const std::function<double(const Eigen::Vector3d)> &f, int num_walks,
                       Eigen::VectorXd &U, Eigen::MatrixXd &U_grad) {
  U.resize(P.rows());
  U_grad.resize(P.rows(), 3);

  igl::parallel_for(P.rows(), [&](int i) {
    double val;
    Eigen::Vector3d grad;
    std::tie(val, grad) = walk_on_sphere_single_point3d(sdf_bc, f, P.row(i), num_walks);
    U[i] = val;
    U_grad.row(i) = grad;
  });
}

void walk_on_spheres(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::VectorXd &B,
                     const Eigen::MatrixXd &P, const std::function<double(const Eigen::Vector3d)> &f, int num_walks,
                     Eigen::VectorXd &U, Eigen::MatrixXd &U_grad) {
  igl::embree::EmbreeIntersector ei;
  ei.init(V.cast<float>(), F);
  auto sdf_bc = [&](const Eigen::Vector3d &p) -> std::pair<double, double> {
    int closest_face;
    float distance;
    ei.distance(p.cast<float>(), distance, closest_face);
    double bc;
    {
      std::array<Eigen::RowVector3d, 3> triangle;
      triangle[0] = V.row(F(closest_face, 0));
      triangle[1] = V.row(F(closest_face, 1));
      triangle[2] = V.row(F(closest_face, 2));
      Eigen::RowVector3d uv;
      igl::barycentric_coordinates(Eigen::RowVector3d(p), triangle[0], triangle[1], triangle[2], uv);
      Eigen::Vector3d gx;
      gx[0] = B[F(closest_face, 0)];
      gx[1] = B[F(closest_face, 1)];
      gx[2] = B[F(closest_face, 2)];
      bc = uv.dot(gx);
    }
    return std::make_pair(double(distance), bc);
  };
  walk_on_spheres3d(sdf_bc, P, f, num_walks, U, U_grad);
}