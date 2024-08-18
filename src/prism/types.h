#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <memory>

namespace prism {

using Eigen::Matrix3d;
using Eigen::MatrixX3d;
using Eigen::Vector3d;
using Eigen::VectorXd;
using SparseMatrix = Eigen::SparseMatrix<double>;

enum class Coord { X, Y, Z };

template <typename T>
using SharedPtr = std::shared_ptr<T>;

} // namespace prism