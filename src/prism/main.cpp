#include <iostream>
#include <Eigen/Dense>
#include <fmt/core.h>
#include <unvpp/unvpp.h>

int main() {
    auto v = unv::Vertex {
        0.0,
        0.0,
        0.0,
    };

    fmt::print("{}", v.x);
    fmt::print("Hello, {}\n", "World!");
}