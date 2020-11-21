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

  Eigen::Vector3d closest_point;
  Eigen::Vector3i closest_face;
  Eigen::VectorXd sqrR(1);
  for (int i = 0; i < P.rows(); i++) {
    Eigen::Vector3d x = P.row(i);
    double sum = 0;

    for (int j = 0; j < nWalks; j++) { // j is a dummy var
      int steps = 0;
      double R = 10000000.;  
      do {
        igl::point_mesh_squared_distance(x, V, F, sqrR, closest_face, closest_point);
        // change to AABB later
        double theta = random(0, M_PI);
        double phi = random(0, 2*M_PI);
        double R = std::sqrt(sqrR[0]);
        Eigen::Vector3d new_direction(R*std::sin(theta)*std::cos(phi),
                                      R*std::sin(theta)*std::sin(phi),
                                      R*std::cos(theta));

        x += new_direction;
        steps++;
      } while (R > eps && steps < maxSteps);

      sum += (B[closest_face[0]]+B[closest_face[1]+B[closest_face[2]]])/3;
      // dunno how author handles it, but this version will work as well.
    }

    U[i] = sum/nWalks;


  }
}