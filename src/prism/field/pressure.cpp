#include "pressure.h"

namespace prism::field {
Pressure::Pressure(std::string name, const mesh::PMesh& mesh, double value)
    : Scalar(std::move(name), mesh, value) {}

Pressure::Pressure(std::string name, const mesh::PMesh& mesh, VectorXd data)
    : Scalar(std::move(name), mesh, std::move(data)) {}

Pressure::Pressure(std::string name, const mesh::PMesh& mesh, VectorXd data, VectorXd face_data)
    : Scalar(std::move(name), mesh, std::move(data), std::move(face_data)) {}

} // namespace prism::field