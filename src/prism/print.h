#pragma once
#include <fmt/color.h>
#include <fmt/core.h>

#include <exception>
#include <string_view>

namespace prism {
void inline warning(std::string_view msg) {
    fmt::print(fg(fmt::color::yellow), "Warning: {}\n", msg);
}

void inline error(std::string_view msg, bool terminate) {
    fmt::print(fg(fmt::color::red), "Error: {}\n", msg);
    if (terminate) {
        std::terminate();
    }
}

} // namespace prism