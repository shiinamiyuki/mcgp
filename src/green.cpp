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
  auto g =  (R-r) / (4*igl::PI*r*R);
  if(!std::isfinite(g)){
    return 0.0;
  }
  return g;
}
Eigen::Vector3d lapdg3d(
  Eigen::Vector3d & x, 
  Eigen::Vector3d & y, 
  double R)
{
  double r = (x-y).norm();
  Eigen::Vector3d dg = (y-x)*(1/pow(r,3) - 1/pow(R,3))/(4*igl::PI);
  if(!std::isfinite(dg[0]) || !std::isfinite(dg[1]) || !std::isfinite(dg[2])){
    return Eigen::Vector3d::Zero();
  }
  return dg;
}

double lapg2d(
  double r,
  double R) 
{
  double g = std::log(R/r)/(2*igl::PI);
  if (!std::isfinite(g)) {
    g = 0;
  }
  return g;
}
Eigen::Vector2d lapdg2d(
  Eigen::Vector2d & x, 
  Eigen::Vector2d & y, 
  double R)
{
  double r = (x-y).norm();
  Eigen::Vector2d dg = (y-x)*(1/pow(r,2) - 1/pow(R,2))/(2*igl::PI);
  if(!std::isfinite(dg[0]) || !std::isfinite(dg[1])){
    return Eigen::Vector2d::Zero();
  }
  return dg;
}