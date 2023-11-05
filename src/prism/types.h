#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cstddef> // std::size_t

namespace prism {

using Eigen::MatrixX3d;
using Eigen::Vector3d;
using Eigen::VectorXd;
using SparseMatrix = Eigen::SparseMatrix<double>;

enum class Compressibility { Compressible, Incompressible };

enum class Coords { X, Y, Z };

} // namespace prism