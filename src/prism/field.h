#pragma once

#include <string>
#include <variant>

#include "mesh/pmesh.h"
#include "types.h"

namespace prism {
// forward declarations
class ScalarField;
class VectorField;


class VectorField {
  public:
    VectorField(std::string name, const mesh::PMesh& mesh);
    VectorField(std::string name, const mesh::PMesh& mesh, MatrixX3d data);
    VectorField(std::string name, const mesh::PMesh& mesh, double value);

    auto inline name() const -> const std::string& { return _name; }
    auto inline data() const -> const MatrixX3d& { return _data; }
    auto inline data() -> MatrixX3d& { return _data; }

    // return the first column of data
    auto inline x() const -> ScalarField;
    auto inline y() const -> ScalarField;
    auto inline z() const -> ScalarField;

  private:
    const mesh::PMesh& _mesh;
    std::string _name;
    MatrixX3d _data;
};

class ScalarField {
  public:
    ScalarField(std::string name, const mesh::PMesh& mesh);
    ScalarField(std::string name, const mesh::PMesh& mesh, VectorXd data);
    ScalarField(std::string name, const mesh::PMesh& mesh, double value);

    auto inline name() const -> const std::string& { return _name; }
    auto inline data() const -> const VectorXd& { return _data; }
    auto inline data() -> VectorXd& { return _data; }
    void inline set_parent_vec_field(VectorField* parent) { _parent_vec_field = parent; }
    void update_parent_vec_field();

  private:
    const mesh::PMesh& _mesh;
    std::string _name;
    VectorXd _data;
    VectorField* _parent_vec_field = nullptr;
};


} // namespace prism