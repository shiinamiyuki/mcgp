#include <igl/PI.h>
#include <limits>
#include <igl/parallel_for.h>
#include <random>
#include <tuple>
#define _USE_MATH_DEFINES
#include <cmath>
#include <green.h>
#include <WoS2d.h>

// returns a random value in the range [rMin,rMax]
static thread_local std::random_device rd;
double random(double rMin, double rMax) {
  std::uniform_real_distribution<double> dist(rMin, rMax);
  return dist(rd);
}


Eigen::Vector2d uniform_circle_sampling() {
  double theta = random(0, 2 * igl::PI);

  Eigen::Vector2d new_direction;
  new_direction[0] = cos(theta);
  new_direction[1] = sin(theta);
  return new_direction;
}
// Eigen::Vector3d uniform_ball_sampling() {
//   return uniform_sphere_sampling() * sqrt(random(0, 1));
// }
// double sphere_volume(double R) { return 4.0 / 3.0 * igl::PI * R * R * R; }


// single point estimator for ∆u = f
std::pair<double, Eigen::Vector2d> walk_on_spheres2d_single_point(
    const std::function<std::pair<double,double>(const Eigen::Vector2d)> &sdf_bc,
    const std::function<double(const Eigen::Vector2d)> &f,
    const Eigen::Vector2d &P,int num_walks) {
  const double bdeps = 0.001;
  const int nWalks = num_walks;
  const int maxSteps = 250;

  double val = 0;
  Eigen::Vector2d grad = Eigen::Vector2d::Zero();
 
  for (int j = 0; j < nWalks; j++) { // j is a dummy var
  Eigen::Vector2d curp = Eigen::Vector2d(P);
  Eigen::Vector2d first_dir;
  double first_R;

    for (int k = 0; k < maxSteps; k++) {
      double sd, bc;
      std::tie(sd, bc) = sdf_bc(P);
      double R = std::abs(sd);

      Eigen::Vector2d newdir = uniform_circle_sampling() * R;

      double ykR = std::sqrt(random(0.,1.)) * R;
      Eigen::Vector2d yk = curp + uniform_circle_sampling() * ykR;
      val = val + igl::PI * R*R * f(yk) * lapg2d(ykR, R);

      if (k == 0) {
        first_dir = newdir / R;
        first_R = R;
        grad += igl::PI* R*R * f(yk) * lapdg2d(curp, yk, R);
      }

      if (R < bdeps) {
        double oldval = val;
        if (j != 0) {
          val += bc - (grad/j).dot(first_dir*first_R);
          grad += 2 * (bc-oldval/j)*first_dir/first_R;
        } else {
          val += bc;
          grad += 2*bc*first_dir/first_R;
        }
        break;
      }

      curp += newdir;
    }

  }
  return std::make_pair(val / nWalks, grad / nWalks);
}

void walk_on_spheres2d(const std::function<std::pair<double,double>(const Eigen::Vector2d)> &sdf_bc,
                     const std::function<double(const Eigen::Vector2d)> &f,
                     const Eigen::MatrixXd &P,
                     int num_walks,
                     Eigen::VectorXd &U,
                     Eigen::MatrixXd &U_grad) {

  U.resize(P.rows());
  U_grad.resize(P.rows(),2);

  igl::parallel_for(P.rows(), [&](int i){
    double val;
    Eigen::Vector2d grad;
    std::tie(val, grad) = walk_on_spheres2d_single_point(sdf_bc, f, P.row(i), num_walks);
    U[i] = val;
    U_grad.row(i) = grad;
  });
}

// void walk_on_spheres(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
//                      const Eigen::VectorXd &B, const Eigen::MatrixXd &P,
//                      Eigen::VectorXd &U) {
//   Eigen::MatrixXd U_grad;
//   return walk_on_spheres(
//       V, F, B, P, [](const Eigen::Vector3d &x) { return 0.0; }, U, U_grad);
// }

