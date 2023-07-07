#pragma once
#include <fmt/color.h>
#include <fmt/core.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>

#include <string_view>

#include "meta.h"

namespace prism {

void inline warn(std::string_view msg) {
    fmt::print(fg(fmt::color::yellow), "Warning: ");
    fmt::print("{}\n", msg);
}

void inline error(std::string_view msg) {
    fmt::print(fg(fmt::color::red), "Error: ");
    fmt::print("{}\n", msg);
}

void inline debug(std::string_view msg) {
#ifdef PRISM_DEBUG
    fmt::print(fg(fmt::color::dark_magenta), "Info: ");
    fmt::print("{}\n", msg);
#endif
}

void inline info(std::string_view msg) {
    fmt::print(fg(fmt::color::dark_cyan), "Info: ");
    fmt::print("{}\n", msg);
}

void inline print_header() {
    fmt::print(fg(fmt::color::cyan), "Prism ");
    fmt::print(
        fmt::format("v{} - A finite volume CFD framework | MIT License\n", PRISM_VERSION_STR));
}

} // namespace prism