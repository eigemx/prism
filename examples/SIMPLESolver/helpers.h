#pragma once

#include <prism/prism.h>

#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>
#include <stdexcept>

using json = nlohmann::json;
using namespace prism;

namespace fs = std::filesystem;

struct FoamFields {
    std::vector<prism::Vector3d> velocity;
    std::vector<double> pressure;
    std::vector<double> temperature;
};

auto inline fileToJson(const fs::path& path) -> json {
    if (!std::filesystem::exists(path)) {
        throw std::runtime_error("File " + path.string() + " does not exist!");
    }
    auto file = std::ifstream(path);
    return json::parse(file);
}

auto inline readFields(const std::filesystem::path& fields_file) -> FoamFields {
    auto doc = fileToJson(fields_file);

    auto velocity = doc["U"];
    auto pressure = doc["p"].get<std::vector<double>>();

    std::vector<prism::Vector3d> velocity_vec;
    for (const auto& v : velocity) {
        velocity_vec.emplace_back(v[0], v[1], v[2]);
    }

    return {velocity_vec, pressure, {}};
}

auto inline readVelocityComponents(const std::vector<prism::Vector3d>& velocity)
    -> std::array<VectorXd, 3> {
    VectorXd u, v, w; // NOLINT
    u.resize(velocity.size());
    v.resize(velocity.size());
    w.resize(velocity.size());

    for (size_t i = 0; i < velocity.size(); i++) {
        u[i] = velocity[i].x();
        v[i] = velocity[i].y();
        w[i] = velocity[i].z();
    }

    return {u, v, w};
}

// define the format you want, you only need one instance of this...
// see https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");

// writing functions taking Eigen types as parameters,
// see https://eigen.tuxfamily.org/dox/TopicFunctionTakingEigenTypes.html
template <typename Derived>
void writeToCSVfile(const std::string& name, const Eigen::MatrixBase<Derived>& matrix) {
    std::ofstream file(name.c_str());
    file << matrix.format(CSVFormat);
    // file.close() is not necessary,
    // desctructur closes file, see https://en.cppreference.com/w/cpp/io/basic_ofstream
}