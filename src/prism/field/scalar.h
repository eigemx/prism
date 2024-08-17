#pragma once

#include "boundary.h"
#include "ifield.h"
#include "prism/exceptions.h"
#include "prism/mesh/utilities.h"
#include "units.h"

namespace prism::field {

class UniformScalar : public IScalar {
  public:
    UniformScalar(std::string name, const mesh::PMesh& mesh, double value);

    auto valueAtCell(std::size_t cell_id) const -> double override;
    auto valueAtCell(const mesh::Cell& cell) const -> double override;
    auto valueAtFace(std::size_t face_id) const -> double override;
    auto valueAtFace(const mesh::Face& face) const -> double override;

  private:
    double _value {0.0};
};

template <typename Units, typename BHManagerSetter>
class GeneralScalar
    : public IScalar,
      public Units,
      public prism::boundary::BHManagersProvider<boundary::IScalarBoundaryHandler> {
  public:
    // Uniform double value constructors
    GeneralScalar(std::string name,
                  const mesh::PMesh& mesh,
                  double value,
                  IVector* parent = nullptr);
    GeneralScalar(std::string name,
                  const mesh::PMesh& mesh,
                  double value,
                  Coord coord,
                  IVector* parent = nullptr);

    // VectorXd cell values constructors
    GeneralScalar(std::string name,
                  const mesh::PMesh& mesh,
                  VectorXd data,
                  IVector* parent = nullptr);
    GeneralScalar(std::string name,
                  const mesh::PMesh& mesh,
                  VectorXd data,
                  Coord coord,
                  IVector* parent = nullptr);

    // VectorXd cell & face values constructors
    GeneralScalar(std::string name,
                  const mesh::PMesh& mesh,
                  VectorXd data,
                  VectorXd face_data,
                  IVector* parent = nullptr);
    GeneralScalar(std::string name,
                  const mesh::PMesh& mesh,
                  VectorXd data,
                  VectorXd face_data,
                  Coord coord,
                  IVector* parent = nullptr);

    GeneralScalar() = delete;
    GeneralScalar(const GeneralScalar&) = default;
    GeneralScalar(GeneralScalar&&) = default;
    auto operator=(const GeneralScalar&) -> GeneralScalar& = default;
    auto operator=(GeneralScalar&&) -> GeneralScalar& = default;
    ~GeneralScalar() override = default;

    // TODO: check that _data is not null before returning, and maybe wrap it in an optional type
    auto inline values() const -> const VectorXd& { return *_data; }
    auto inline values() -> VectorXd& { return *_data; }

    auto inline coord() const noexcept -> std::optional<Coord> override { return _coord; }
    auto inline hasFaceValues() const -> bool override { return _face_data != nullptr; }
    void setFaceValues(VectorXd values);

    auto valueAtCell(std::size_t cell_id) const -> double override;
    auto valueAtCell(const mesh::Cell& cell) const -> double override;
    auto valueAtFace(std::size_t face_id) const -> double override;
    auto valueAtFace(const mesh::Face& face) const -> double override;

    auto parent() -> IVector*;
    void setParent(IVector* parent);

    auto inline operator[](std::size_t i) const -> double { return (*_data)[i]; }
    auto inline operator[](std::size_t i) -> double& { return (*_data)[i]; }

  protected:
    auto valueAtInteriorFace(const mesh::Face& face) const -> double;
    auto valueAtBoundaryFace(const mesh::Face& face) const -> double;

  private:
    void addDefaultHandlers();

    std::shared_ptr<VectorXd> _data = nullptr;
    std::shared_ptr<VectorXd> _face_data = nullptr;
    IVector* _parent = nullptr;
    std::optional<Coord> _coord = std::nullopt;
    BHManagerSetter _setter;
};

class ScalarBHManagerSetter {
  public:
    using IScalarBHManager =
        prism::boundary::BoundaryHandlersManager<boundary::IScalarBoundaryHandler>;

    static void set(IScalarBHManager& manager);
};

using Scalar = GeneralScalar<units::Measurable, ScalarBHManagerSetter>;

template <typename Units, typename BHManagerSetter>
GeneralScalar<Units, BHManagerSetter>::GeneralScalar(std::string name,
                                                     const mesh::PMesh& mesh,
                                                     double value,
                                                     IVector* parent)
    : IScalar(std::move(name), mesh),
      _data(std::make_shared<VectorXd>(VectorXd::Ones(mesh.nCells()) * value)),
      _parent(parent) {
    spdlog::debug("Creating scalar field: '{}' with double value = {}", this->name(), value);
    addDefaultHandlers();
}

template <typename Units, typename BHManagerSetter>
GeneralScalar<Units, BHManagerSetter>::GeneralScalar(std::string name,
                                                     const mesh::PMesh& mesh,
                                                     double value,
                                                     Coord coord,
                                                     IVector* parent)
    : IScalar(std::move(name), mesh),
      _data(std::make_shared<VectorXd>(VectorXd::Ones(mesh.nCells()) * value)),
      _coord(coord),
      _parent(parent) {
    spdlog::debug("Creating scalar field: '{}' (as {}-coordinate) with double value = {}",
                  this->name(),
                  coordToStr(coord),
                  value);
    addDefaultHandlers();
}

template <typename Units, typename BHManagerSetter>
GeneralScalar<Units, BHManagerSetter>::GeneralScalar(std::string name,
                                                     const mesh::PMesh& mesh,
                                                     VectorXd data,
                                                     IVector* parent)
    : IScalar(std::move(name), mesh),
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

template <typename Units, typename BHManagerSetter>
GeneralScalar<Units, BHManagerSetter>::GeneralScalar(std::string name,
                                                     const mesh::PMesh& mesh,
                                                     VectorXd data,
                                                     Coord coord,
                                                     IVector* parent)
    : IScalar(std::move(name), mesh),
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

template <typename Units, typename BHManagerSetter>
GeneralScalar<Units, BHManagerSetter>::GeneralScalar(std::string name,
                                                     const mesh::PMesh& mesh,
                                                     VectorXd data,
                                                     VectorXd face_data,
                                                     IVector* parent)
    : IScalar(std::move(name), mesh),
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

template <typename Units, typename BHManagerSetter>
GeneralScalar<Units, BHManagerSetter>::GeneralScalar(std::string name,
                                                     const mesh::PMesh& mesh,
                                                     VectorXd data,
                                                     VectorXd face_data,
                                                     Coord coord,
                                                     IVector* parent)
    : IScalar(std::move(name), mesh),
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

template <typename Units, typename BHManagerSetter>
void GeneralScalar<Units, BHManagerSetter>::setFaceValues(VectorXd values) {
    if (values.size() != mesh().nFaces()) {
        throw std::runtime_error(
            fmt::format("prism::field::GeneralScalar<Units, BHManagerProvider, "
                        "BHManagerSetter>::setFaceValues(): cannot set face values for "
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

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::valueAtCell(const mesh::Cell& cell) const -> double {
    return valueAtCell(cell.id());
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::valueAtCell(std::size_t cell_id) const -> double {
    assert(_data != nullptr);           // NOLINT
    assert(cell_id < mesh().nCells()); // NOLINT
    return (*_data)[cell_id];
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::valueAtFace(std::size_t face_id) const -> double {
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

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::valueAtFace(const mesh::Face& face) const -> double {
    return valueAtFace(face.id());
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::valueAtInteriorFace(const mesh::Face& face) const
    -> double {
    assert(face.is_interior()); // NOLINT
    const auto& owner = mesh().cell(face.owner());
    const auto& neighbor = mesh().cell(face.neighbor().value());

    const auto gc = mesh::geo_weight(owner, neighbor, face);
    double val = gc * (*_data)[owner.id()];
    val += (1 - gc) * (*_data)[neighbor.id()];

    return val;
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::valueAtBoundaryFace(const mesh::Face& face) const
    -> double {
    const auto& patch = mesh().boundary_patch(face);
    const auto& bc = patch.getBoundaryCondition(name());

    auto handler = this->boundaryHandlersManager().getHandler(bc.kindString());

    if (handler == nullptr) {
        throw error::NonImplementedBoundaryCondition(
            fmt::format("prism::field::GeneralScalar<Units, BHManagerProvider, "
                        "BHManagerSetter>::valueAtBoundaryFace() {}",
                        name()),
            patch.name(),
            bc.kindString());
    }

    return handler->get(*this, face);
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::parent() -> IVector* {
    return _parent;
}

template <typename Units, typename BHManagerSetter>
void GeneralScalar<Units, BHManagerSetter>::setParent(IVector* parent) {
    _parent = parent;
    // TODO: check parent and component names consistency
}

template <typename Units, typename BHManagerSetter>
void GeneralScalar<Units, BHManagerSetter>::addDefaultHandlers() {
    _setter.set(this->boundaryHandlersManager());
}

void inline ScalarBHManagerSetter::set(IScalarBHManager& manager) {
    spdlog::debug(
        "prism::field::ScalarBHManagerSetter::set(): adding default boundary handlers for a "
        "scalar "
        "field instance");
    manager.addHandler<field::boundary::Fixed<Scalar>>();
    // manager.addHandler<field::boundary::VelocityInlet<Scalar>>();
    manager.addHandler<field::boundary::Empty<Scalar>>();
    manager.addHandler<field::boundary::Symmetry<Scalar>>();
    manager.addHandler<field::boundary::Outlet<Scalar>>();
    manager.addHandler<field::boundary::FixedGradient<Scalar>>();
    // this->boundaryHandlersManager().addHandler<field::boundary::NoSlip<Scalar>>();
}
} // namespace prism::field