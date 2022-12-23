#include <fmt/core.h>
#include <unvpp/unvpp.h>


#include <Eigen/Dense>
#include <iostream>

int main() {
    auto v = unv::Vertex{
        0.0,
        0.0,
        0.0,
    };

    auto m = v;
    std::string x;
    x.clear();

    fmt::print("{}", v.x);
    fmt::print("Hello, {}\n", "World!");
}