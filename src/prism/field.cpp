#include "field.h"

#include <stdexcept>

#include "fmt/core.h"
#include "print.h"
#include "prism/mesh/pmesh.h"

namespace prism {

void check_field_name(const std::string& name) {
    if (name.empty()) {
        throw std::runtime_error("Cannot create a Field with an empty name.");
    }
}

// ScalarField constructors
ScalarField::ScalarField(std::string name, const mesh::PMesh& mesh)
    : _mesh(&mesh),
      _name(std::move(name)),
      _data(std::make_shared<VectorXd>(VectorXd::Zero(mesh.n_cells()))) {
    check_field_name(_name);
}


ScalarField::ScalarField(std::string name, const mesh::PMesh& mesh, VectorXd data)
    : _mesh(&mesh), _name(std::move(name)), _data(std::make_shared<VectorXd>(std::move(data))) {
    check_field_name(_name);
}


ScalarField::ScalarField(std::string name, const mesh::PMesh& mesh, double value)
    : _mesh(&mesh),
      _name(std::move(name)),
      _data(std::make_shared<VectorXd>(VectorXd::Ones(mesh.n_cells()) * value)) {
    check_field_name(_name);
}


ScalarField::ScalarField(std::string name,
                         const mesh::PMesh& mesh,
                         std::shared_ptr<VectorXd> data)
    : _mesh(&mesh), _name(std::move(name)), _data(std::move(data)) {
    check_field_name(_name);
}


auto ScalarField::clone() const -> ScalarField {
    return {_name, *_mesh, *_data};
}


auto ScalarField::map(CellMapper* mapper) -> ScalarField& {
    for (std::size_t i = 0; i < _mesh->n_cells(); ++i) {
        data()[i] = mapper(_mesh->cell(i));
    }
    return *this;
}


auto ScalarField::map(CoordinatesMapper* mapper) -> ScalarField& {
    const auto n_cells = _mesh->n_cells();
    for (std::size_t i = 0; i < n_cells; ++i) {
        const auto& cell = _mesh->cell(i);
        const auto& center = cell.center();
        data()[i] = mapper(center.x(), center.y(), center.z());
    }
    return *this;
}

// VectorField constructors and methods
VectorField::VectorField(std::string name, const mesh::PMesh& mesh)
    : _mesh(&mesh),
      _name(std::move(name)),
      _x(std::make_shared<VectorXd>(mesh.n_cells())),
      _y(std::make_shared<VectorXd>(mesh.n_cells())),
      _z(std::make_shared<VectorXd>(mesh.n_cells())) {}


VectorField::VectorField(std::string name, const mesh::PMesh& mesh, const Vector3d& data)
    : _mesh(&mesh),
      _name(std::move(name)),
      _x(std::make_shared<VectorXd>(VectorXd::Ones(mesh.n_cells()) * data(0))),
      _y(std::make_shared<VectorXd>(VectorXd::Ones(mesh.n_cells()) * data(1))),
      _z(std::make_shared<VectorXd>(VectorXd::Ones(mesh.n_cells()) * data(2))) {}


VectorField::VectorField(std::string name, const mesh::PMesh& mesh, const MatrixX3d& data)
    : _mesh(&mesh),
      _name(std::move(name)),
      _x(std::make_shared<VectorXd>(data.col(0))),
      _y(std::make_shared<VectorXd>(data.col(1))),
      _z(std::make_shared<VectorXd>(data.col(2))) {}


VectorField::VectorField(std::string name, const mesh::PMesh& mesh, double value)
    : _mesh(&mesh),
      _name(std::move(name)),
      _x(std::make_shared<VectorXd>(VectorXd::Ones(mesh.n_cells()) * value)),
      _y(std::make_shared<VectorXd>(VectorXd::Ones(mesh.n_cells()) * value)),
      _z(std::make_shared<VectorXd>(VectorXd::Ones(mesh.n_cells()) * value)) {}


auto VectorField::data() const -> MatrixX3d {
    MatrixX3d data(_mesh->cells().size(), 3);
    data.col(0) = *_x;
    data.col(1) = *_y;
    data.col(2) = *_z;
    return data;
}


auto VectorField::x() -> ScalarField {
    return {_name + "_x", *_mesh, _x};
}


auto VectorField::y() -> ScalarField {
    return {_name + "_y", *_mesh, _y};
}


auto VectorField::z() -> ScalarField {
    return {_name + "_z", *_mesh, _z};
}
} // namespace prism