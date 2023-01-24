#pragma once
#include <fmt/color.h>
#include <fmt/core.h>

#include <string_view>

namespace prism {
void inline warning(std::string_view msg) {
    fmt::print(fg(fmt::color::yellow), "Warning: ");
    fmt::print("{}\n", msg);
}

void inline error(std::string_view msg) {
    fmt::print(fg(fmt::color::red), "Error: ");
    fmt::print("{}\n", msg);
}

} // namespace prism