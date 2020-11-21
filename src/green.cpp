#include "green.h"
#include <igl/PI.h>
double lapg3d(
  Eigen::Vector3d & x, 
  Eigen::Vector3d & y, 
  double R) 
{
  double r = (x-y).norm();
  return lapg3d(r, R);
}

double lapg3d(
  double r, 
  double R)
{
  return (R-r) / (4*igl::PI*r*R);
}