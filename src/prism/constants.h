#pragma once

#include <string_view>

using std::string_view_literals::operator""sv;

static constexpr auto PRISM_VERSION_STR = "0.1.0"sv;

constexpr double PRISM_EPSILON = 1e-10;
constexpr double PRISM_PI = 3.14159265358979323846;
