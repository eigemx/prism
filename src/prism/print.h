#pragma once
#include <fmt/color.h>
#include <fmt/core.h>

#include <string_view>

#include "meta.h"

namespace prism {

using fmt::format;
using fmt::print;

void inline warn(std::string_view msg) {
    print(fg(fmt::color::yellow), "Warning: ");
    print("{}\n", msg);
}

void inline error(std::string_view msg) {
    print(fg(fmt::color::red), "Error: ");
    print("{}\n", msg);
}

void inline debug(std::string_view msg) {
#ifdef PRISM_DEBUG
    print(fg(fmt::color::dark_magenta), "Info: ");
    print("{}\n", msg);
#endif
}

void inline info(std::string_view msg) {
    print(fg(fmt::color::dark_cyan), "Info: ");
    print("{}\n", msg);
}

void inline print_header() {
    print(fg(fmt::color::cyan), "Prism ");
    print("v{} - A finite volume CFD framework | MIT License\n", PRISM_VERSION_STR);
}

} // namespace prism