#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
#include <string>

namespace prism {

using String = std::string;
using f32 = float;
using f64 = double;
using std::size_t;

using Eigen::Matrix3d;
using Eigen::MatrixX3d;
using Eigen::Vector3d;
using Tensor3d = Eigen::Matrix3d;
using Eigen::VectorXd;

using SparseMatrix = Eigen::SparseMatrix<double>;

enum class VectorCoord : std::uint8_t { X, Y, Z };

enum Sign : std::uint8_t { Positive, Negative };

template <typename T>
using UniquePtr = std::unique_ptr<T>;

template <typename T>
using SharedPtr = std::shared_ptr<T>;

template <typename T>
using Optional = std::optional<T>;

template <typename T>
using Ref = std::reference_wrapper<T>;

inline constexpr auto NullOption = std::nullopt;
} // namespace prism
