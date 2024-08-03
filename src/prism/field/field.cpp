
#include "field.h"

#include <fmt/core.h>

#include <cassert>
#include <cstddef>
#include <memory>
#include <stdexcept>

#include "prism/exceptions.h"
#include "prism/field/boundary.h"
#include "prism/mesh/boundary.h"
#include "prism/mesh/cell.h"
#include "prism/mesh/face.h"
#include "prism/mesh/pmesh.h"
#include "prism/mesh/utilities.h"
#include "spdlog/spdlog.h"

namespace prism::field {
namespace detail {
void check_field_name(const std::string& name) {
    if (name.empty()) {
        throw std::runtime_error("Cannot create a Field with an empty name.");
    }
}

void check_mesh(const mesh::PMesh& mesh) {
    if (mesh.cells().empty() || mesh.faces().empty() || mesh.boundary_patches().empty()) {
        throw std::runtime_error("Cannot create a field over an empty mesh.");
    }
}
} // namespace detail

UniformScalar::UniformScalar(std::string name, const mesh::PMesh& mesh, double value)
    : IField(std::move(name), mesh), _value(value) {
    spdlog::debug(
        "Creating uniform scalar field: '{}' with double value = {}", this->name(), value);
}

auto UniformScalar::value_at_cell(std::size_t cell_id) const -> double { // NOLINT
    return _value;
}

auto UniformScalar::value_at_cell(const mesh::Cell& cell) const -> double { // NOLINT
    return _value;
}

auto UniformScalar::value_at_face(std::size_t face_id) const -> double { // NOLINT
    return _value;
}

auto UniformScalar::value_at_face(const mesh::Face& face) const -> double { // NOLINT
    return _value;
}

Scalar::Scalar(std::string name, const mesh::PMesh& mesh, double value, Vector* parent)
    : IField(std::move(name), mesh),
      _data(std::make_shared<VectorXd>(VectorXd::Ones(mesh.n_cells()) * value)),
      _parent(parent) {
    spdlog::debug("Creating scalar field: '{}' with double value = {}", this->name(), value);
    register_default_handlers();
}

Scalar::Scalar(std::string name, const mesh::PMesh& mesh, VectorXd data, Vector* parent)
    : IField(std::move(name), mesh),
      _data(std::make_shared<VectorXd>(std::move(data))),
      _parent(parent) {
    if (_data->size() != mesh.n_cells()) {
        throw std::runtime_error(fmt::format(
            "field::Scalar() cannot create a scalar field '{}' given a vector that has a "
            "different size than mesh's cell count.",
            this->name()));
    }

    spdlog::debug("Creating scalar field: '{}' with a vector data of size = {}",
                  this->name(),
                  _data->size());
    register_default_handlers();
}

Scalar::Scalar(std::string name,
               const mesh::PMesh& mesh,
               VectorXd data,
               VectorXd face_data,
               Vector* parent)
    : IField(std::move(name), mesh),
      _data(std::make_shared<VectorXd>(std::move(data))),
      _face_data(std::make_shared<VectorXd>(std::move(face_data))),
      _parent(parent) {
    if (_data->size() != mesh.n_cells()) {
        throw std::runtime_error(fmt::format(
            "field::Scalar() cannot create a scalar field '{}' given a vector that has a "
            "different size than mesh's cell count.",
            this->name()));
    }

    if (_face_data->size() != mesh.n_faces()) {
        throw std::runtime_error(
            fmt::format("field::Scalar() cannot create a scalar field '{}' given a face data "
                        "vector that has a different size than mesh's faces count.",
                        this->name()));
    }

    spdlog::debug(
        "Creating scalar field: '{}' with a cell data vector of size = {} and face data vector "
        "of size = {}",
        this->name(),
        _data->size(),
        _face_data->size());

    register_default_handlers();
}

void Scalar::set_face_values(VectorXd values) {
    if (values.size() != mesh().n_faces()) {
        throw std::runtime_error(fmt::format(
            "ScalarField::set_face_values(): cannot set face values for scalar field {}, to a "
            "face data vector having a different size that field's faces count.",
            name()));
    }

    if (has_face_data()) {
        spdlog::debug("Setting new face values to scalar field '{}', discarding old face values.",
                      name());
    }

    _face_data = std::make_shared<VectorXd>(std::move(values));
}

auto Scalar::value_at_cell(const mesh::Cell& cell) const -> double {
    return value_at_cell(cell.id());
}

auto Scalar::value_at_cell(std::size_t cell_id) const -> double {
    assert(_data != nullptr);           // NOLINT
    assert(cell_id < mesh().n_cells()); // NOLINT
    return (*_data)[cell_id];
}

auto Scalar::value_at_face(std::size_t face_id) const -> double {
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

auto Scalar::value_at_face(const mesh::Face& face) const -> double {
    return value_at_face(face.id());
}

auto Scalar::value_at_interior_face(const mesh::Face& face) const -> double {
    assert(face.is_interior()); // NOLINT
    const auto& owner = mesh().cell(face.owner());
    const auto& neighbor = mesh().cell(face.neighbor().value());

    const auto gc = mesh::geo_weight(owner, neighbor, face);
    double val = gc * (*_data)[owner.id()];
    val += (1 - gc) * (*_data)[neighbor.id()];

    return val;
}

auto Scalar::value_at_boundary_face(const mesh::Face& face) const -> double {
    const auto& patch = mesh().boundary_patch(face);
    const auto& bc = patch.get_bc(name());

    auto handler = _bh_manager.get_handler(bc.kind_string());

    if (handler == nullptr) {
        throw error::NonImplementedBoundaryCondition(
            fmt::format("ScalarField({})::value_at_boundary_face()", name()),
            patch.name(),
            bc.kind_string());
    }

    return handler->get(*this, face);
}

auto Scalar::parent() -> std::optional<Vector> {
    if (_parent == nullptr) {
        return std::nullopt;
    }
    return *_parent;
}

void Scalar::register_default_handlers() {
    _bh_manager.add_handler<field::boundary::Fixed>();
    _bh_manager.add_handler<field::boundary::VelocityInlet>();
    _bh_manager.add_handler<field::boundary::Empty>();
    _bh_manager.add_handler<field::boundary::Symmetry>();
    _bh_manager.add_handler<field::boundary::Outlet>();
    _bh_manager.add_handler<field::boundary::FixedGradient>();
}

Tensor::Tensor(std::string name, const mesh::PMesh& mesh, double value)
    : IField(std::move(name), mesh) {
    spdlog::debug("Creating tensor field: '{}' with double value = {}", this->name(), value);
    const std::size_t n_cells = this->mesh().n_cells();
    _data.reserve(n_cells);
    for (std::size_t i = 0; i < n_cells; ++i) {
        _data.emplace_back(Matrix3d::Ones() * value);
    }
}

Tensor::Tensor(std::string name, const mesh::PMesh& mesh, const Matrix3d& data)
    : IField(std::move(name), mesh) {
    spdlog::debug("Creating a uniform tensor field: '{}' given a Matrix3d object", this->name());

    const std::size_t n_cells = this->mesh().n_cells();
    _data.reserve(n_cells);
    for (std::size_t i = 0; i < n_cells; ++i) {
        _data.push_back(data);
    }
}

Tensor::Tensor(std::string name, const mesh::PMesh& mesh, std::vector<Matrix3d> data)
    : IField(std::move(name), mesh), _data(std::move(data)) {
    spdlog::debug("Creating a  tensor field: '{}' given a vector of Matrix3d objects",
                  this->name());

    if (_data.size() != mesh.n_cells()) {
        throw std::runtime_error(
            fmt::format("field::Tensor() cannot create a tensor field '{}' given a vector of "
                        "Matrix3d that has a different size than mesh's cell count.",
                        this->name()));
    }
}

auto Tensor::value_at_cell(std::size_t cell_id) const -> Matrix3d {
    assert(cell_id < mesh().n_cells());
    return _data[cell_id];
}

auto Tensor::value_at_cell(const mesh::Cell& cell) const -> Matrix3d {
    return _data[cell.id()];
}

auto Tensor::value_at_face(std::size_t face_id) const -> Matrix3d {
    const mesh::Face& face = this->mesh().face(face_id);
    return value_at_face(face);
}

auto Tensor::value_at_face(const mesh::Face& face) const -> Matrix3d {
    const auto& mesh = this->mesh();
    const mesh::Cell& owner = mesh.cell(face.owner());

    if (face.is_boundary()) {
        spdlog::warn(
            "field::Tensor::value_at_face() was called on a boundary face (face id = {}). "
            "Returning value of the tensor field at owner cell.",
            face.id());

        return _data[owner.id()];
    }
    const mesh::Cell& neighbor = mesh.cell(face.neighbor().value());
    const double gc = mesh::geo_weight(owner, neighbor, face);

    return (gc * _data[owner.id()]) + ((1 - gc) * _data[neighbor.id()]);
}

auto Tensor::operator[](std::size_t i) -> Matrix3d& {
    return _data[i];
}

auto Tensor::operator[](std::size_t i) const -> const Matrix3d& {
    return _data[i];
}

} // namespace prism::field