#pragma once

#include <limits>
#include <numbers>
#include <string_view>

#include "types.h"

namespace prism {

inline constexpr std::string_view VERSION = "0.1.0";
inline constexpr f64 EPSILON = std::numeric_limits<f64>::epsilon();
inline constexpr f64 PI = std::numbers::pi;

} // namespace prism
