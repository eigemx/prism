#pragma once
#include <fmt/color.h>
#include <fmt/core.h>

#include <string_view>

#include "version.h"

namespace prism {
void inline warning(std::string_view msg) {
    fmt::print(fg(fmt::color::yellow), "Warning: ");
    fmt::print("{}\n", msg);
}

void inline error(std::string_view msg) {
    fmt::print(fg(fmt::color::red), "Error: ");
    fmt::print("{}\n", msg);
}

void inline print_header() {
    fmt::print(fg(fmt::color::cyan), "Prism ");
    fmt::print("v{} - A finite volume CFD framework\n\n", PRISM_VERSION_STR);
}

} // namespace prism