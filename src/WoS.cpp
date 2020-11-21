#include <igl/point_mesh_squared_distance.h>
#include <limits>
#include <cmath>

// returns a random value in the range [rMin,rMax]
// Copied from http://www.cs.cmu.edu/~kmcrane/Projects/MonteCarloGeometryProcessing/WoSLaplace2D.cpp.html
double random( double rMin, double rMax ) {
   const double rRandMax = 1./(double)RAND_MAX;
   double u = rRandMax*(double)rand();
   return u*(rMax-rMin) + rMin;
}

Eigen::VectorXd random(int N, double rMin, double rMax) {
  Eigen::VectorXd result(N);
  for (int i = 0; i < N; i++) {
    result[i] = random(rMin, rMax);
  }
  return result;
}

void walk_on_spheres(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXd & B,
  const Eigen::MatrixXd & P,
  Eigen::VectorXd & U)
{
  U.resize(P.rows());

  const double eps = 0.01; 
  const int nWalks = 128;
  const int maxSteps = 32;
  // they can be set as inputs in the future

  Eigen::MatrixXd closest_points;
  Eigen::MatrixXi closest_faces;
  Eigen::VectorXd sqrRs;

  Eigen::VectorXd sum = Eigen::VectorXd::Zero(P.rows());

  for (int j = 0; j < nWalks; j++) { // j is a dummy var
    Eigen::MatrixXd PP(P); //tmp
    int steps = 0;
    double R = 10000000.;  
    do {
      igl::point_mesh_squared_distance(PP, V, F, sqrRs, closest_faces, closest_points);
      // change to AABB later
      Eigen::VectorXd theta = random(P.rows(), 0, M_PI);
      Eigen::VectorXd phi = random(P.rows(), 0, 2*M_PI);
      Eigen::VectorXd R = sqrRs.array().sqrt();

      Eigen::MatrixXd new_direction(P.rows(), 3);
      new_direction.col(0) = R.array() * theta.array().sin() * phi.array().cos();
      new_direction.col(1) = R.array() * theta.array().sin() * phi.array().sin();
      new_direction.col(2) = R.array() * theta.array().cos();

      PP += new_direction;
      steps++;
    } while (R > eps && steps < maxSteps);

    // the following loop can be improved by eigen, probably.
    for (int i = 0; i < closest_faces.rows(); i++) {
      sum[i] += (B[closest_faces(i, 0)]+B[closest_faces(i, 1)+B[closest_faces(i, 2)]])/3;
      // dunno how author handles it, but this version will work as well.
    } 
  }

  U = sum/nWalks;

}