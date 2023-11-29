
#include "field.h"

#include <fmt/core.h>

#include <cassert>
#include <cstddef>
#include <stdexcept>

#include "exceptions.h"
#include "prism/mesh/boundary.h"
#include "prism/mesh/cell.h"
#include "prism/mesh/face.h"
#include "prism/mesh/pmesh.h"
#include "prism/mesh/utilities.h"

namespace prism {

void inline check_field_name(const std::string& name) {
    if (name.empty()) {
        throw std::runtime_error("Cannot create a Field with an empty name.");
    }
}

void inline check_mesh(const mesh::PMesh& mesh) {
    if (mesh.cells().empty() || mesh.faces().empty() || mesh.boundary_patches().empty()) {
        throw std::runtime_error("Cannot create a field over an empty mesh.");
    }
}

Field::Field(std::string name, const mesh::PMesh& mesh) : _name(std::move(name)), _mesh(&mesh) {
    check_field_name(_name);
    check_mesh(mesh);
}

ScalarField::ScalarField(std::string name, const mesh::PMesh& mesh, double value)
    : Field(std::move(name), mesh),
      _data(std::make_shared<VectorXd>(VectorXd::Ones(mesh.n_cells()) * value)) {}

ScalarField::ScalarField(std::string name, const mesh::PMesh& mesh, VectorXd data)
    : Field(std::move(name), mesh), _data(std::make_shared<VectorXd>(std::move(data))) {}

ScalarField::ScalarField(std::string name,
                         const mesh::PMesh& mesh,
                         VectorXd data,
                         VectorXd face_data)
    : Field(std::move(name), mesh),
      _data(std::make_shared<VectorXd>(std::move(data))),
      _face_data(std::make_shared<VectorXd>(std::move(face_data))) {}

auto ScalarField::clone() const -> ScalarField {
    if (has_face_data()) {
        return {name(), mesh(), *_data, *_face_data};
    }
    return {name(), mesh(), *_data};
}


auto ScalarField::map(CellMapper* mapper) -> ScalarField& {
    for (std::size_t i = 0; i < mesh().n_cells(); ++i) {
        data()[i] = mapper(mesh().cell(i));
    }
    return *this;
}

auto ScalarField::map(CoordinatesMapper* mapper) -> ScalarField& {
    const auto n_cells = mesh().n_cells();
    for (std::size_t i = 0; i < n_cells; ++i) {
        const auto& cell = mesh().cell(i);
        const auto& center = cell.center();
        data()[i] = mapper(center.x(), center.y(), center.z());
    }
    return *this;
}

auto ScalarField::value_at_cell(std::size_t cell_id) const -> double {
    return (*_data)[cell_id];
}

auto ScalarField::value_at_cell(const mesh::Cell& cell) const -> double {
    return value_at_cell(cell.id());
}

auto ScalarField::value_at_face(std::size_t face_id) const -> double {
    if (has_face_data()) {
        // Face data were calculataed for us before calling the constructor,
        // just return the value
        return (*_face_data)[face_id];
    }

    // We need to interpolate the value of the field at the face
    const auto& face = mesh().face(face_id);

    if (face.is_interior()) {
        return value_at_interior_face(face);
    }

    return value_at_boundary_face(face);
}

auto ScalarField::value_at_face(const mesh::Face& face) const -> double {
    return value_at_face(face.id());
}

auto ScalarField::value_at_interior_face(const mesh::Face& face) const -> double {
    const auto& owner = mesh().cell(face.owner());
    const auto& neighbor = mesh().cell(face.neighbor().value());
    const auto gc = mesh::geo_weight(owner, neighbor, face);

    double val = gc * (*_data)[owner.id()];
    val += (1 - gc) * (*_data)[neighbor.id()];

    return val;
}

auto ScalarField::value_at_boundary_face(const mesh::Face& face) const -> double {
    const auto& patch = mesh().boundary_patch(face);
    const auto& bc = patch.get_bc(name());

    switch (bc.bc_type()) {
        case mesh::BoundaryConditionType::Fixed:
        case mesh::BoundaryConditionType::Inlet: {
            return patch.get_scalar_bc(name());
        }

        // TODO: We return the field value of an empty face the same value as its owner cell.
        // This makes many schemes and gradient methods to work without a special check for an
        // empty face (like LeastSquares), check if this assumption is correct.
        case mesh::BoundaryConditionType::Empty:
        case mesh::BoundaryConditionType::Symmetry:
        case mesh::BoundaryConditionType::Outlet: {
            return (*_data)[face.owner()];
        }

        case mesh::BoundaryConditionType::FixedGradient: {
            Vector3d grad_at_boundary = patch.get_vector_bc(name());
            const auto& owner = mesh().cell(face.owner());
            Vector3d e = face.center() - owner.center();
            double d_Cf = e.norm();
            e = e / e.norm();
            grad_at_boundary = grad_at_boundary * d_Cf;

            return grad_at_boundary.dot(e) + value_at_cell(owner);
        }

        default: {
            throw error::NonImplementedBoundaryCondition(
                fmt::format("ScalarField({})::value_at_boundary_face()", name()),
                patch.name(),
                bc.bc_type_str());
        }
    }
}

VectorField::VectorField(std::string name, const mesh::PMesh& mesh, double value)
    : Field(std::move(name), mesh),
      _x(this->name() + "_x", mesh, value),
      _y(this->name() + "_y", mesh, value),
      _z(this->name() + "_z", mesh, value) {}

VectorField::VectorField(std::string name, const mesh::PMesh& mesh, const Vector3d& data)
    : Field(std::move(name), mesh),
      _x(this->name() + "_x", mesh, data[0]),
      _y(this->name() + "_y", mesh, data[1]),
      _z(this->name() + "_z", mesh, data[2]) {}

VectorField::VectorField(std::string name,
                         const mesh::PMesh& mesh,
                         const std::array<ScalarField, 3>& fields)
    : Field(std::move(name), mesh), _x(fields[0]), _y(fields[1]), _z(fields[2]) {
    // check mesh consistency
    for (const auto& field : fields) {
        if (&mesh != &field.mesh()) {
            throw std::runtime_error(
                fmt::format("VectorField constructor was given a ScalarField component with name "
                            "`{}` that is defined over a different mesh",
                            field.name()));
        }
    }

    // check sub-fields naming consistency
    if ((_x.name() != (this->name() + "_x")) || (_y.name() != (this->name() + "_y")) ||
        (_z.name() != (this->name() + "_z"))) {
        throw std::runtime_error(fmt::format(
            "All VectorField component names should end with '_x', '_y' or '_z'. VectorField "
            "constructor for `{}` vector field was given the following ScalarFields names: '{}', "
            "'{}', '{}",
            this->name(),
            _x.name(),
            _y.name(),
            _z.name()));
    }
}

auto VectorField::value_at_cell(std::size_t cell_id) const -> Vector3d {
    return operator[](cell_id);
}

auto VectorField::value_at_cell(const mesh::Cell& cell) const -> Vector3d {
    return value_at_cell(cell.id());
}

auto VectorField::value_at_face(std::size_t face_id) const -> Vector3d {
    return {_x.value_at_face(face_id), _y.value_at_face(face_id), _z.value_at_face(face_id)};
}

auto VectorField::value_at_face(const mesh::Face& face) const -> Vector3d {
    return value_at_face(face.id());
}

auto VectorField::has_face_data() const -> bool {
    return _x.has_face_data() && _y.has_face_data() && _z.has_face_data();
}

auto VectorField::operator[](std::size_t i) const -> Vector3d {
    return {_x.data()[i], _y.data()[i], _z.data()[i]};
}

TensorField::TensorField(std::string name, const mesh::PMesh& mesh, double value)
    : Field(std::move(name), mesh) {
    init_data_vec();

    const std::size_t n_cells = this->mesh().n_cells();
    for (std::size_t i = 0; i < n_cells; ++i) {
        _data.emplace_back(Matrix3d::Ones() * value);
    }
}

TensorField::TensorField(std::string name, const mesh::PMesh& mesh, Matrix3d data)
    : Field(std::move(name), mesh) {
    init_data_vec();

    const std::size_t n_cells = this->mesh().n_cells();
    for (std::size_t i = 0; i < n_cells; ++i) {
        _data.emplace_back(std::move(data));
    }
}

TensorField::TensorField(std::string name, const mesh::PMesh& mesh, std::vector<Matrix3d> data)
    : Field(std::move(name), mesh), _data(std::move(data)) {}

auto TensorField::value_at_cell(std::size_t cell_id) const -> const Matrix3d& {
    return _data[cell_id];
}

auto TensorField::value_at_cell(const mesh::Cell& cell) const -> const Matrix3d& {
    return _data[cell.id()];
}

auto TensorField::value_at_face(std::size_t face_id) const -> Matrix3d {
    const auto& mesh = this->mesh();
    const mesh::Face& face = mesh.face(face_id);
    return value_at_face(face);
}

auto TensorField::value_at_face(const mesh::Face& face) const -> Matrix3d {
    const auto& mesh = this->mesh();
    const mesh::Cell& owner = mesh.cell(face.owner());
    const mesh::Cell& neighbor = mesh.cell(face.neighbor().value());
    const double gc = mesh::geo_weight(owner, neighbor, face);

    return (gc * _data[owner.id()]) + ((1 - gc) * _data[neighbor.id()]);
}

auto TensorField::at(std::size_t i, std::size_t j, std::size_t k) -> double& {
    // i -> cell id
    // j -> row number
    // k -> column number
    return _data[i].coeffRef(j, k);
}

auto TensorField::at(std::size_t i, std::size_t j, std::size_t k) const -> double {
    return _data[i](j, k);
}

auto TensorField::operator[](std::size_t i) -> Matrix3d& {
    return _data[i];
}

auto TensorField::operator[](std::size_t i) const -> const Matrix3d& {
    return _data[i];
}

void TensorField::init_data_vec() {
    assert(mesh().n_cells() > 0 && "TensorField::init_data_vec() has a zero cells mesh");
    _data.reserve(mesh().n_cells());
}

} // namespace prism