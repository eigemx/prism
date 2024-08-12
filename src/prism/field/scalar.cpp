#include "scalar.h"

#include <spdlog/spdlog.h>

#include "prism/exceptions.h"
#include "prism/mesh/utilities.h"
#include "vector.h"

auto inline coordToStr(prism::Coord coord) -> std::string {
    switch (coord) {
        case prism::Coord::X: {
            return "x";
        }
        case prism::Coord::Y: {
            return "y";
        }
        case prism::Coord::Z: {
            return "z";
        }
    }
}

namespace prism::field {

UniformScalar::UniformScalar(std::string name, const mesh::PMesh& mesh, double value)
    : IField(std::move(name), mesh), _value(value) {
    spdlog::debug(
        "Creating uniform scalar field: '{}' with double value = {}", this->name(), value);
}

auto UniformScalar::valueAtCell(std::size_t cell_id) const -> double { // NOLINT
    return _value;
}

auto UniformScalar::valueAtCell(const mesh::Cell& cell) const -> double { // NOLINT
    return _value;
}

auto UniformScalar::valueAtFace(std::size_t face_id) const -> double { // NOLINT
    return _value;
}

auto UniformScalar::valueAtFace(const mesh::Face& face) const -> double { // NOLINT
    return _value;
}

Scalar::Scalar(std::string name, const mesh::PMesh& mesh, double value, IVector* parent)
    : IField(std::move(name), mesh),
      _data(std::make_shared<VectorXd>(VectorXd::Ones(mesh.nCells()) * value)),
      _parent(parent) {
    spdlog::debug("Creating scalar field: '{}' with double value = {}", this->name(), value);
    addDefaultHandlers();
}

Scalar::Scalar(std::string name,
               const mesh::PMesh& mesh,
               double value,
               Coord coord,
               IVector* parent)
    : IField(std::move(name), mesh),
      _data(std::make_shared<VectorXd>(VectorXd::Ones(mesh.nCells()) * value)),
      _coord(coord),
      _parent(parent) {
    spdlog::debug("Creating scalar field: '{}' (as {}-coordinate) with double value = {}",
                  this->name(),
                  coordToStr(coord),
                  value);
    addDefaultHandlers();
}

Scalar::Scalar(std::string name, const mesh::PMesh& mesh, VectorXd data, IVector* parent)
    : IField(std::move(name), mesh),
      _data(std::make_shared<VectorXd>(std::move(data))),
      _parent(parent) {
    if (_data->size() != mesh.nCells()) {
        throw std::runtime_error(fmt::format(
            "field::Scalar() cannot create a scalar field '{}' given a vector that has a "
            "different size than mesh's cell count.",
            this->name()));
    }

    spdlog::debug("Creating scalar field: '{}' with a cell vector data of size = {}",
                  this->name(),
                  _data->size());
    addDefaultHandlers();
}

Scalar::Scalar(std::string name,
               const mesh::PMesh& mesh,
               VectorXd data,
               Coord coord,
               IVector* parent)
    : IField(std::move(name), mesh),
      _data(std::make_shared<VectorXd>(std::move(data))),
      _coord(coord),
      _parent(parent) {
    if (_data->size() != mesh.nCells()) {
        throw std::runtime_error(fmt::format(
            "field::Scalar() cannot create a scalar field '{}' given a vector that has a "
            "different size than mesh's cell count.",
            this->name()));
    }

    spdlog::debug(
        "Creating scalar field: '{}' (as {}-coordinate) with a cell vector data of size = {}",
        this->name(),
        coordToStr(coord),
        _data->size());
    addDefaultHandlers();
}

Scalar::Scalar(std::string name,
               const mesh::PMesh& mesh,
               VectorXd data,
               VectorXd face_data,
               IVector* parent)
    : IField(std::move(name), mesh),
      _data(std::make_shared<VectorXd>(std::move(data))),
      _face_data(std::make_shared<VectorXd>(std::move(face_data))),
      _parent(parent) {
    if (_data->size() != mesh.nCells()) {
        throw std::runtime_error(fmt::format(
            "field::Scalar() cannot create a scalar field '{}' given a vector that has a "
            "different size than mesh's cell count.",
            this->name()));
    }

    if (_face_data->size() != mesh.nFaces()) {
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

    addDefaultHandlers();
}

Scalar::Scalar(std::string name,
               const mesh::PMesh& mesh,
               VectorXd data,
               VectorXd face_data,
               Coord coord,
               IVector* parent)
    : IField(std::move(name), mesh),
      _data(std::make_shared<VectorXd>(std::move(data))),
      _face_data(std::make_shared<VectorXd>(std::move(face_data))),
      _coord(coord),
      _parent(parent) {
    if (_data->size() != mesh.nCells()) {
        throw std::runtime_error(fmt::format(
            "field::Scalar() cannot create a scalar field '{}' given a vector that has a "
            "different size than mesh's cell count.",
            this->name()));
    }

    if (_face_data->size() != mesh.nFaces()) {
        throw std::runtime_error(
            fmt::format("field::Scalar() cannot create a scalar field '{}' given a face data "
                        "vector that has a different size than mesh's faces count.",
                        this->name()));
    }

    spdlog::debug(
        "Creating scalar field: '{}' (as {}-coordinate) with a cell data vector of size = {} and "
        "face data vector "
        "of size = {}",
        this->name(),
        coordToStr(coord),
        _data->size(),
        _face_data->size());

    addDefaultHandlers();
}

void Scalar::setFaceValues(VectorXd values) {
    if (values.size() != mesh().nFaces()) {
        throw std::runtime_error(
            fmt::format("prism::field::Scalar::setFaceValues(): cannot set face values for "
                        "scalar field {}, to a "
                        "face data vector having a different size that field's faces count.",
                        name()));
    }

    if (hasFaceValues()) {
        spdlog::debug("Setting new face values to scalar field '{}', discarding old face values.",
                      name());
    }

    _face_data = std::make_shared<VectorXd>(std::move(values));
}

auto Scalar::valueAtCell(const mesh::Cell& cell) const -> double {
    return valueAtCell(cell.id());
}

auto Scalar::valueAtCell(std::size_t cell_id) const -> double {
    assert(_data != nullptr);           // NOLINT
    assert(cell_id < mesh().n_cells()); // NOLINT
    return (*_data)[cell_id];
}

auto Scalar::valueAtFace(std::size_t face_id) const -> double {
    if (hasFaceValues()) {
        // Face data were calculataed for us before calling the constructor,
        // just return the value
        return (*_face_data)[face_id];
    }

    // We need to interpolate the value of the field at the face
    const auto& face = mesh().face(face_id);

    if (face.is_interior()) {
        return valueAtInteriorFace(face);
    }

    return valueAtBoundaryFace(face);
}

auto Scalar::valueAtFace(const mesh::Face& face) const -> double {
    return valueAtFace(face.id());
}

auto Scalar::valueAtInteriorFace(const mesh::Face& face) const -> double {
    assert(face.is_interior()); // NOLINT
    const auto& owner = mesh().cell(face.owner());
    const auto& neighbor = mesh().cell(face.neighbor().value());

    const auto gc = mesh::geo_weight(owner, neighbor, face);
    double val = gc * (*_data)[owner.id()];
    val += (1 - gc) * (*_data)[neighbor.id()];

    return val;
}

auto Scalar::valueAtBoundaryFace(const mesh::Face& face) const -> double {
    const auto& patch = mesh().boundary_patch(face);
    const auto& bc = patch.getBoundaryCondition(name());

    auto handler = _bh_manager.getHandler(bc.kindString());

    if (handler == nullptr) {
        throw error::NonImplementedBoundaryCondition(
            fmt::format("prism::field::Scalar::valueAtBoundaryFace() {}", name()),
            patch.name(),
            bc.kindString());
    }

    return handler->get(*this, face);
}

auto Scalar::parent() -> IVector* {
    return _parent;
}

void Scalar::setParent(IVector* parent) {
    _parent = parent;
    // TODO: check parent and component names consistency
}

void Scalar::addDefaultHandlers() {
    _bh_manager.addHandler<field::boundary::Fixed>();
    _bh_manager.addHandler<field::boundary::Inlet>();
    _bh_manager.addHandler<field::boundary::Empty>();
    _bh_manager.addHandler<field::boundary::Symmetry>();
    _bh_manager.addHandler<field::boundary::Outlet>();
    _bh_manager.addHandler<field::boundary::FixedGradient>();
    _bh_manager.addHandler<field::boundary::NoSlip>();
}
} // namespace prism::field