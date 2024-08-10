#pragma once

#include "prism/field/scalar.h"

namespace prism::field {

class Pressure : public Scalar, public units::Pascal {
  public:
    Pressure(std::string name, const mesh::PMesh& mesh, double value);
    Pressure(std::string name, const mesh::PMesh& mesh, VectorXd data);
    Pressure(std::string name, const mesh::PMesh& mesh, VectorXd data, VectorXd face_data);
};

} // namespace prism::field