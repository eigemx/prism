#pragma once

#include <array>
#include <unvpp/unvpp.h>

namespace prism::mesh {

constexpr std::size_t hex_vertices = 8;
constexpr std::size_t hex_faces = 6;
constexpr std::size_t quad_face_vertices = 4;

// Array of 6 hex faces, each hex face is any array of 4 vertices
template <typename T = std::array<std::array<std::size_t, quad_face_vertices>,
                                  hex_faces>>
auto inline hex_cell_faces(std::vector<std::size_t>& v) -> T {
    T array;
    array[1] = {v[0], v[4], v[7], v[3]};
    array[2] = {v[3], v[7], v[6], v[2]};
    array[3] = {v[0], v[1], v[5], v[4]};
    array[4] = {v[4], v[5], v[6], v[7]};
    array[5] = {v[0], v[3], v[2], v[1]};
    return array;
}
} // namespace prism::mesh