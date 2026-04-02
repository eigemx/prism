#pragma once

#include "prism/types.h"

namespace prism::field::detail {

auto inline coordToIndex(VectorCoord coord) -> std::uint8_t {
    switch (coord) {
        case VectorCoord::X: return 0;
        case VectorCoord::Y: return 1;
        case VectorCoord::Z: return 2;
        default: break;
    }
    throw std::invalid_argument("prism::field::detail::coordToIndex(): Invalid coordinate value");
}

auto inline coordToStr(VectorCoord coord) -> std::string {
    switch (coord) {
        case prism::VectorCoord::X: {
            return "x";
        }
        case prism::VectorCoord::Y: {
            return "y";
        }
        case prism::VectorCoord::Z: {
            return "z";
        }
        default: {
            return "";
        }
    }
}

} // namespace prism::field::detail