#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace prism {

using Eigen::Matrix3d;
using Eigen::MatrixX3d;
using Eigen::Vector3d;
using Eigen::VectorXd;
using SparseMatrix = Eigen::SparseMatrix<double>;

enum class Compressibility { Compressible, Incompressible };

enum class Coord { X, Y, Z };

} // namespace prism