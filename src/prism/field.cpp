#include "field.h"

#include <fmt/core.h>

#include <stdexcept>

#include "prism/mesh/pmesh.h"

namespace prism {

void inline check_field_name(const std::string& name) {
    if (name.empty()) {
        throw std::runtime_error("Cannot create a Field with an empty name.");
    }
}

// ScalarField constructors
ScalarField::ScalarField(std::string name, const mesh::PMesh& mesh, double value)
    : _mesh(&mesh),
      _name(std::move(name)),
      _data(std::make_shared<VectorXd>(VectorXd::Ones(mesh.n_cells()) * value)) {
    check_field_name(_name);
}


ScalarField::ScalarField(std::string name, const mesh::PMesh& mesh, VectorXd data)
    : _mesh(&mesh), _name(std::move(name)), _data(std::make_shared<VectorXd>(std::move(data))) {
    check_field_name(_name);
}


ScalarField::ScalarField(std::string name,
                         const mesh::PMesh& mesh,
                         VectorXd data,
                         VectorXd face_data)
    : _mesh(&mesh),
      _name(std::move(name)),
      _data(std::make_shared<VectorXd>(std::move(data))),
      _face_data(std::make_shared<VectorXd>(std::move(face_data))) {
    check_field_name(_name);
}

auto ScalarField::clone() const -> ScalarField {
    return {_name, *_mesh, *_data, *_face_data};
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
VectorField::VectorField(std::string name, const mesh::PMesh& mesh, double value)
    : _mesh(&mesh),
      _name(std::move(name)),
      _x(_name + "_x", mesh, value),
      _y(_name + "_y", mesh, value),
      _z(_name + "_z", mesh, value) {}


VectorField::VectorField(std::string name, const mesh::PMesh& mesh, const Vector3d& data)
    : _mesh(&mesh),
      _name(std::move(name)),
      _x(_name + "_x", mesh, data[0]),
      _y(_name + "_y", mesh, data[1]),
      _z(_name + "_z", mesh, data[2]) {}


VectorField::VectorField(std::string name,
                         const mesh::PMesh& mesh,
                         const std::array<ScalarField, 3>& fields)
    : _mesh(&mesh), _name(std::move(name)), _x(fields[0]), _y(fields[1]), _z(fields[2]) {}


auto VectorField::value_at_face(std::size_t face_id) const -> Vector3d {
    if (!has_face_data()) {
        throw std::runtime_error(
            fmt::format("VectorField::face_data(): face data for vector field `{}` is not "
                        "available, however, face_data() was called.",
                        _name));
    }

    return {_x.face_data()[face_id], _y.face_data()[face_id], _z.face_data()[face_id]};
}


auto VectorField::has_face_data() const -> bool {
    return _x.face_data().size() > 0 && _y.face_data().size() > 0 && _z.face_data().size() > 0;
}

} // namespace prism