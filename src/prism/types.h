#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace prism {

using Eigen::MatrixX3d;
using Eigen::Vector3d;
using Eigen::VectorXd;
using SparseMatrix = Eigen::SparseMatrix<double>;

constexpr double PRISM_EPSILON = 1e-10;
constexpr double PRISM_PI = 3.14159265358979323846;

} // namespace prism