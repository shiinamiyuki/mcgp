#include <mls.h>
#include <Eigen/LU>
#include <Eigen/SVD>
#include <iostream>
#include <igl/pinv.h>

template <>
double moving_least_squares<1>(const Eigen::VectorXd &f, const Eigen::MatrixXd &xs, const Eigen::Vector3d &p) {
  constexpr int nBasis = 4;
  using Basis = Eigen::Matrix<double, nBasis, 1>;
  auto basis = [](Eigen::Vector3d p) -> Basis {
    Basis b;
    b[0] = 1.0;
    b[1] = p.x();
    b[2] = p.y();
    b[3] = p.z();
    return b;
  };
  auto theta = [](double d) {
    const auto eps = 0.00001;
    return 1.0 / (d * d + eps * eps);
  };
  Eigen::Matrix<double, nBasis, nBasis> A;
  A.setZero();
  Eigen::Matrix<double, nBasis, 1> b;
  b.setZero();
  for (int i = 0; i < xs.rows(); i++) {
    Eigen::Vector3d x = Eigen::Vector3d(xs.row(i)) - p;
    A += theta(x.norm()) * basis(x) * basis(x).transpose();
    b += theta(x.norm()) * basis(x) * f[i];
  }
  Eigen::Matrix<double, nBasis, 1> c = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
  igl::pinv(A, b);
  return c[0];
}

template <>
double moving_least_squares<2>(const Eigen::VectorXd &f, const Eigen::MatrixXd &xs, const Eigen::Vector3d &p) {
  constexpr int nBasis = 10;
  using Basis = Eigen::Matrix<double, nBasis, 1>;
  auto basis = [](Eigen::Vector3d p) -> Basis {
    Basis b;
    b[0] = 1.0;
    b[1] = p.x();
    b[2] = p.y();
    b[3] = p.z();
    b[4] = p.x() * p.y();
    b[5] = p.y() * p.z();
    b[6] = p.z() * p.x();
    b[7] = p.x() * p.x();
    b[8] = p.y() * p.y();
    b[9] = p.z() * p.z();
    return b;
  };
  auto theta = [](double d) {
    const auto eps = 1e-9;
    return 1.0 / (d * d + eps * eps);
  };
  Eigen::Matrix<double, nBasis, nBasis> A;
  A.setZero();
  Eigen::Matrix<double, nBasis, 1> b;
  b.setZero();
  for (int i = 0; i < xs.rows(); i++) {
    Eigen::Vector3d x = Eigen::Vector3d(xs.row(i)) - p;
    A += theta(x.norm()) * basis(x) * basis(x).transpose();
    b += theta(x.norm()) * basis(x) * f[i];
  }
//   std::cout << "A:\n" << A << "\nb:\n" << b << std::endl;
  Eigen::Matrix<double, nBasis, 1> c = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
  igl::pinv(A, b);
  return c[0];
}

template <>
double moving_least_squares<3>(const Eigen::VectorXd &f, const Eigen::MatrixXd &xs, const Eigen::Vector3d &p) {
  constexpr int nBasis = 21;
  using Basis = Eigen::Matrix<double, nBasis, 1>;
  auto basis = [](Eigen::Vector3d p) -> Basis {
    Basis b;
    b[0] = 1.0;
    b[1] = p.x();
    b[2] = p.y();
    b[3] = p.z();
    b[4] = p.x() * p.y();
    b[5] = p.y() * p.z();
    b[6] = p.z() * p.x();
    b[7] = p.x() * p.x();
    b[8] = p.y() * p.y();
    b[9] = p.z() * p.z();
    b[10] = p.x() * p.x() * p.x();
    b[11] = p.x() * p.x() * p.y();
    b[12] = p.x() * p.x() * p.z();
    b[13] = p.x() * p.y() * p.y();
    b[15] = p.x() * p.y() * p.z();
    b[16] = p.x() * p.z() * p.z();
    b[17] = p.y() * p.y() * p.y();
    b[18] = p.y() * p.y() * p.z();
    b[19] = p.y() * p.z() * p.z();
    b[20] = p.z() * p.z() * p.z();
    return b;
  };
  auto theta = [](double d) {
    const auto eps = 0.1;
    return 1.0 / (d * d + eps * eps);
  };
  Eigen::Matrix<double, nBasis, nBasis> A;
  A.setZero();
  Eigen::Matrix<double, nBasis, 1> b;
  b.setZero();
  for (int i = 0; i < xs.rows(); i++) {
    Eigen::Vector3d x = Eigen::Vector3d(xs.row(i)) - p;
    A += theta(x.norm()) * basis(x) * basis(x).transpose();
    b += theta(x.norm()) * basis(x) * f[i];
  }
  Eigen::Matrix<double, nBasis, 1> c = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
  igl::pinv(A, b);
  return c[0];
}