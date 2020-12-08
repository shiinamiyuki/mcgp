#include <igl/AABB.h>
#include <igl/PI.h>
#include <igl/barycentric_coordinates.h>
#include <igl/barycentric_interpolation.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/centroid.h>
#include <igl/octree.h>
#include <igl/knn.h>
#include <limits>
#include <igl/parallel_for.h>
#include <random>
#include <tuple>
#define _USE_MATH_DEFINES
#include <cmath>
#include <green.h>
#include <WoS.h>
#include <mls.h>

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

template <class T>
struct VarianceTracker {
  double mean, m2;
  int count = 0;
  void update(T value) {
    if (count == 0) {
      mean = value;
      m2 = T(0.0);
    } else {
      auto delta = value - mean;
      mean += delta / T(count + 1);
      m2 += delta * (value - mean);
    }
    count++;
  }
  T variance() const {
    if (count < 2) {
      return T(-1.0);
    }
    return m2;
  }
  T mean_variance() const {
    if (count < 2) {
      return T(-1.0);
    }
    return m2 / T(count * count);
  }
};

struct WoSEstimator {
  double u_var;
  // VarianceTracker<Eigen::Vector3d> u_grad_var;
  double u = 0.0;
  Eigen::Vector3d grad;
};

// single point estimator for ∆u = f
WoSEstimator
walk_on_sphere_single_point3d(const std::function<std::pair<double, double>(const Eigen::Vector3d)> &sdf_bc,
                              const std::function<double(const Eigen::Vector3d)> &f, const Eigen::Vector3d &P,
                              const Eigen::Vector3d &center, double Rmax, int num_walks) {
  const double eps = 0.001;
  const int nWalks = num_walks;
  const int maxSteps = 32;
  VarianceTracker<double> u_var, u_grad_var;
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
          u = u - oldgrad.dot(first_R * first_direction);
          grad = (u - oldsum / j) * first_direction * 3 / first_R + grad;
        } else {
          grad = (u * first_direction * 3 / first_R + grad);
        }

        if (std::isfinite(u)) {
          u_var.update(u);
          sum += u;
        }
        if (std::isfinite(grad.sum())) {
          sumgrad += grad;
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
      if ((x - center).squaredNorm() > Rmax * Rmax) {
        break;
      }
      if (steps >= maxSteps) {
        break;
        // auto continue_prob = std::fmin(1.0, R_last / R) * 0.95;
        // if (std::isfinite(continue_prob) && random(0, 1) < continue_prob) {
        //   k /= continue_prob;
        // } else {
        //   break;
        // }
        // R_last = std::min(R_last, R);
      }
    }
  }
  WoSEstimator es;
  es.grad = 2.0 * sumgrad / nWalks;
  es.u = sum / nWalks;
  es.u_var = u_var.variance();
  return es;
}

void walk_on_spheres3d(const std::function<std::pair<double, double>(const Eigen::Vector3d)> &sdf_bc,
                       const Eigen::MatrixXd &P, const std::function<double(const Eigen::Vector3d)> &f,
                       const Eigen::Vector3d &center, double Rmax, int num_walks, Eigen::VectorXd &U,
                       Eigen::MatrixXd &U_grad) {
  U.resize(P.rows());
  U_grad.resize(P.rows(), 3);

  igl::parallel_for(P.rows(), [&](int i) {
    double val;
    Eigen::Vector3d grad;
    auto es = walk_on_sphere_single_point3d(sdf_bc, f, P.row(i), center, Rmax, num_walks);
    val = es.u;
    grad = es.grad;
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
  Eigen::Vector3d center;
  igl::centroid(V, F, center);
  double Rmax = 0.0;
  for (int i = 0; i < V.rows(); i++) {
    Rmax = std::max((V.row(i).transpose() - center).norm(), Rmax);
  }
  walk_on_spheres3d(sdf_bc, P, f, center, Rmax * 2.0, num_walks, U, U_grad);
}

void walk_on_spheres3d_region(const std::function<std::pair<double, double>(const Eigen::Vector3d)> &sdf_bc,
                              const std::function<double(const Eigen::Vector3d)> &f,
                              const std::function<Eigen::Vector3d(Eigen::Vector3d)> &region,
                              const Eigen::Vector3d &center, double Rmax, size_t walks_per_point, size_t total_walks,
                              WoSPointCloud &point_cloud) {
  point_cloud.resize(0);
  const auto walks_per_candidate = walks_per_point;
  size_t accumluate_walks = 0;
  int consecutive_zeros = 0;
  int max_iter = 8 * (total_walks / walks_per_candidate / 1024) + 1;
  int iter = 0;
  while (accumluate_walks < total_walks && iter < max_iter) {
    iter++;
    if (consecutive_zeros >= 16)
      break;
    std::vector<std::vector<int>> O_PI;
    Eigen::MatrixXi O_CH;
    Eigen::MatrixXd O_CN;
    Eigen::VectorXd O_W;
    size_t n_max_candidate = std::min<size_t>(1024, (total_walks - accumluate_walks) / walks_per_candidate);
    Eigen::MatrixXd P(n_max_candidate, 3);
    std::vector<Eigen::Vector3d> candidates;
    for (size_t i = 0; i < n_max_candidate; i++) {
      Eigen::Vector3d candidate = region(Eigen::Vector3d(random(0, 1), random(0, 1), random(0, 1)));
      P.row(i) = candidate;
    }
    if (point_cloud.n_points() > 16) {
      igl::octree(point_cloud.P, O_PI, O_CH, O_CN, O_W);
      Eigen::MatrixXi I;
      Eigen::VectorXi neighbors;
      std::cout << "computing knn for " << P.rows() << " points" << std::endl;
      igl::knn(P, point_cloud.P, 8, O_PI, O_CH, O_CN, O_W, I);
      for (int i = 0; i < P.rows(); i++) {
        Eigen::Vector3d x = P.row(i);
        neighbors = I.row(i);
        auto approx = [&](int idx) {
          return point_cloud.U[idx] + point_cloud.U_grad.row(idx).dot(x - Eigen::Vector3d(point_cloud.P.row(idx)));
        };
        VarianceTracker<double> u_var;
        double u_mean = 0.0;
        for (int ni = 0; ni < neighbors.size(); ni++) {
          int idx = neighbors[ni];
          auto f = approx(idx);
          u_var.update(f);
          u_mean += f / neighbors.size();
        }
        // std::cout << u_var.variance() << " " << u_mean << std::endl;
        if (u_mean == 0.0 && u_var.variance() > 0.01) {
          candidates.push_back(x);
        } else if (u_var.variance() / std::abs(u_mean) > 0.01) {
          candidates.push_back(x);
        }
      }
    } else {
      for (int i = 0; i < P.rows(); i++) {
        Eigen::Vector3d x = P.row(i);
        candidates.push_back(x);
      }
    }
    std::cout << "candidates: " << candidates.size() << std::endl;
    if (candidates.empty()) {
      consecutive_zeros++;
      continue;
    } else {
      consecutive_zeros = 0;
    }
    P.resize(candidates.size(), 3);
    // std::cout << P << std::endl;
    Eigen::MatrixXd U_grad;
    Eigen::VectorXd U;
    for (size_t i = 0; i < candidates.size(); i++) {
      P.row(i) = candidates[i];
    }
    walk_on_spheres3d(sdf_bc, P, f, center, Rmax, walks_per_candidate, U, U_grad);
    auto n_prev = point_cloud.n_points();
    point_cloud.resize(point_cloud.n_points() + candidates.size());
    for (size_t i = 0; i < candidates.size(); i++) {
      point_cloud.P.row(n_prev + i) = P.row(i);
      point_cloud.U_grad.row(n_prev + i) = U_grad.row(i);
      point_cloud.U[n_prev + i] = U[i];
    }
    accumluate_walks += walks_per_candidate * candidates.size();
    std::cout << "point cloud size: " << point_cloud.n_points() << std::endl;
    // std::cout << point_cloud.U << std::endl;
  }
}

void wos_point_cloud_interpolate(const WoSPointCloud &point_cloud, const Eigen::MatrixXd &P, Eigen::VectorXd &U,
                                 Eigen::MatrixXd &U_grad) {
  U.resize(P.rows());
  U_grad.resize(P.rows(), 3);
  Eigen::VectorXd grad_x = point_cloud.U_grad.col(0);
  Eigen::VectorXd grad_y = point_cloud.U_grad.col(1);
  Eigen::VectorXd grad_z = point_cloud.U_grad.col(2);
  std::vector<std::vector<int>> O_PI;
  Eigen::MatrixXi O_CH;
  Eigen::MatrixXd O_CN;
  Eigen::VectorXd O_W;
  igl::octree(point_cloud.P, O_PI, O_CH, O_CN, O_W);
  Eigen::MatrixXi I;
  igl::knn(P, point_cloud.P, 16, O_PI, O_CH, O_CN, O_W, I);
  igl::parallel_for(P.rows(), [&](int i) {
    Eigen::Matrix<double, 16, 1> knn_gx, knn_gy, knn_gz, knn_u;
    Eigen::Matrix<double, 16, 3> knn_p;
    for (int j = 0; j < 16; j++) {
      knn_gx[j] = grad_x[I(i, j)];
      knn_gy[j] = grad_y[I(i, j)];
      knn_gz[j] = grad_z[I(i, j)];
      knn_u[j] = point_cloud.U[I(i, j)];
      knn_p.row(j) = point_cloud.P.row(I(i, j));
    }
    U[i] = moving_least_squares(knn_u, knn_p, P.row(i));
    U_grad(i, 0) = moving_least_squares(knn_gx, knn_p, P.row(i));
    U_grad(i, 1) = moving_least_squares(knn_gy, knn_p, P.row(i));
    U_grad(i, 2) = moving_least_squares(knn_gz, knn_p, P.row(i));
  });
}