#include "field.h"

#include "print.h"

namespace prism {

// ScalarField constructors
ScalarField::ScalarField(std::string name, const mesh::PMesh& mesh)
    : _mesh(mesh),
      _name(std::move(name)),
      _data(std::make_shared<VectorXd>(VectorXd::Zero(mesh.cells().size()))) {}


ScalarField::ScalarField(std::string name, const mesh::PMesh& mesh, VectorXd data)
    : _mesh(mesh), _name(std::move(name)), _data(std::make_shared<VectorXd>(std::move(data))) {}


ScalarField::ScalarField(std::string name, const mesh::PMesh& mesh, double value)
    : _mesh(mesh),
      _name(std::move(name)),
      _data(std::make_shared<VectorXd>(VectorXd::Ones(mesh.cells().size()) * value)) {}


ScalarField::ScalarField(std::string name,
                         const mesh::PMesh& mesh,
                         std::shared_ptr<VectorXd> data)
    : _mesh(mesh), _name(std::move(name)), _data(std::move(data)) {}


// VectorField constructors and methods
VectorField::VectorField(std::string name, const mesh::PMesh& mesh)
    : _mesh(mesh),
      _name(std::move(name)),
      _x(std::make_shared<VectorXd>(mesh.cells().size())),
      _y(std::make_shared<VectorXd>(mesh.cells().size())),
      _z(std::make_shared<VectorXd>(mesh.cells().size())) {}


VectorField::VectorField(std::string name, const mesh::PMesh& mesh, const Vector3d& data)
    : _mesh(mesh),
      _name(std::move(name)),
      _x(std::make_shared<VectorXd>(VectorXd::Ones(mesh.cells().size()) * data(0))),
      _y(std::make_shared<VectorXd>(VectorXd::Ones(mesh.cells().size()) * data(1))),
      _z(std::make_shared<VectorXd>(VectorXd::Ones(mesh.cells().size()) * data(2))) {}

VectorField::VectorField(std::string name, const mesh::PMesh& mesh, const MatrixX3d& data)
    : _mesh(mesh),
      _name(std::move(name)),
      _x(std::make_shared<VectorXd>(data.col(0))),
      _y(std::make_shared<VectorXd>(data.col(1))),
      _z(std::make_shared<VectorXd>(data.col(2))) {}

VectorField::VectorField(std::string name, const mesh::PMesh& mesh, double value)
    : _mesh(mesh),
      _name(std::move(name)),
      _x(std::make_shared<VectorXd>(VectorXd::Ones(mesh.cells().size()) * value)),
      _y(std::make_shared<VectorXd>(VectorXd::Ones(mesh.cells().size()) * value)),
      _z(std::make_shared<VectorXd>(VectorXd::Ones(mesh.cells().size()) * value)) {}

auto VectorField::data() const -> MatrixX3d {
    MatrixX3d data(_mesh.cells().size(), 3);
    data.col(0) = *_x;
    data.col(1) = *_y;
    data.col(2) = *_z;
    return data;
}

auto VectorField::x() -> ScalarField {
    return {_name + "_x", _mesh, _x};
}

auto VectorField::y() -> ScalarField {
    return {_name + "_y", _mesh, _y};
}

auto VectorField::z() -> ScalarField {
    return {_name + "_z", _mesh, _z};
}

} // namespace prism