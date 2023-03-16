#include "field.h"

#include "print.h"

namespace prism {

// ScalarField constructors
ScalarField::ScalarField(std::string name, const mesh::PMesh& mesh)
    : _mesh(mesh), _name(std::move(name)), _data(VectorXd::Zero(mesh.cells().size())) {}

ScalarField::ScalarField(std::string name, const mesh::PMesh& mesh, VectorXd data)
    : _mesh(mesh), _name(std::move(name)), _data(std::move(data)) {}

ScalarField::ScalarField(std::string name, const mesh::PMesh& mesh, double value)
    : _mesh(mesh),
      _name(std::move(name)),
      _data(VectorXd::Constant(mesh.cells().size(), value)) {}

void ScalarField::update_parent_vec_field() {
    if (_parent_vec_field == nullptr) {
        return;
    }

    if (_name == _parent_vec_field->name() + "_x") {
        _parent_vec_field->data().col(0) = _data;
    }

    else if (_name == _parent_vec_field->name() + "_y") {
        _parent_vec_field->data().col(1) = _data;
    }

    else if (_name == _parent_vec_field->name() + "_z") {
        _parent_vec_field->data().col(2) = _data;
    }

    else {
        throw std::runtime_error(
            format("ScalarField::update_parent_vec_field: {} is not a component of {}",
                   _name,
                   _parent_vec_field->name()));
    }
}


// VectorField constructors and methods
VectorField::VectorField(std::string name, const mesh::PMesh& mesh)
    : _mesh(mesh), _name(std::move(name)), _data(MatrixX3d::Zero(mesh.cells().size(), 3)) {}

VectorField::VectorField(std::string name, const mesh::PMesh& mesh, MatrixX3d data)
    : _mesh(mesh), _name(std::move(name)), _data(std::move(data)) {}

VectorField::VectorField(std::string name, const mesh::PMesh& mesh, double value)
    : _mesh(mesh), _name(std::move(name)) {
    _data = MatrixX3d::Constant(mesh.cells().size(), 3, value);
}

auto VectorField::x() const -> ScalarField {
    return {_name + "_x", _mesh, _data.col(0)};
}

auto VectorField::y() const -> ScalarField {
    return {_name + "_y", _mesh, _data.col(1)};
}

auto VectorField::z() const -> ScalarField {
    return {_name + "_z", _mesh, _data.col(2)};
}

} // namespace prism