#pragma once

#include <algorithm>

#include "ifield.h"
#include "prism/exceptions.h"
#include "prism/gradient/gradient.h"
#include "prism/log.h"
#include "prism/mesh/utilities.h"
#include "scalar_boundary.h"
#include "units.h"

/// TODO: we need to make sure that constructors are not leaving _data uninitialized, and we can
/// avoid checks for it later in member functions.

namespace prism::field {

// forward declaration for GeneralScalar
template <typename Units, typename BHManagerSetter>
class GeneralScalar;

// forward declaration for GeneralVector
template <typename Component, typename BHManagerSetter>
class GeneralVector;

class UniformScalar : public IScalar {
  public:
    UniformScalar(std::string name, const SharedPtr<mesh::PMesh>& mesh, double value);

    auto valueAtCell(std::size_t cell_id) const -> double override;
    auto valueAtCell(const mesh::Cell& cell) const -> double override;
    auto valueAtFace(std::size_t face_id) const -> double override;
    auto valueAtFace(const mesh::Face& face) const -> double override;

    auto gradAtFace(const mesh::Face& face) const -> Vector3d override;
    auto gradAtCell(const mesh::Cell& cell) const -> Vector3d override;
    auto gradAtCellStored(const mesh::Cell& cell) const -> Vector3d override;

    template <IScalarBased ScalarField>
    auto operator*(const ScalarField& other) -> ScalarField;

    template <IVectorBased VectorField>
    auto operator*(const VectorField& other) -> VectorField;

  private:
    double _value {0.0};
};

template <typename Units, typename BHManagerSetter>
class GeneralScalar
    : public IScalar,
      public Units,
      public prism::boundary::BHManagerProvider<boundary::scalar::IScalarBoundaryHandler> {
  public:
    GeneralScalar(std::string name,
                  const SharedPtr<mesh::PMesh>& mesh,
                  double value,
                  IVector* parent = nullptr);

    GeneralScalar(std::string name,
                  const SharedPtr<mesh::PMesh>& mesh,
                  double value,
                  Coord coord,
                  IVector* parent = nullptr);

    GeneralScalar(std::string name,
                  const SharedPtr<mesh::PMesh>& mesh,
                  VectorXd data,
                  IVector* parent = nullptr);

    GeneralScalar(std::string name,
                  const SharedPtr<mesh::PMesh>& mesh,
                  VectorXd data,
                  Coord coord,
                  IVector* parent = nullptr);

    GeneralScalar(std::string name,
                  const SharedPtr<mesh::PMesh>& mesh,
                  VectorXd data,
                  VectorXd face_data,
                  IVector* parent = nullptr);

    GeneralScalar(std::string name,
                  const SharedPtr<mesh::PMesh>& mesh,
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

    auto values() const -> const VectorXd&;
    auto values() -> VectorXd&;

    auto coord() const noexcept -> std::optional<Coord> override;
    auto hasFaceValues() const -> bool override;
    void setFaceValues(VectorXd values);
    void clearFaceValues();

    auto valueAtCell(std::size_t cell_id) const -> double override;
    auto valueAtCell(const mesh::Cell& cell) const -> double override;
    auto valueAtFace(std::size_t face_id) const -> double override;
    auto valueAtFace(const mesh::Face& face) const -> double override;

    auto parent() const -> const IVector*;
    void setParent(const IVector* parent);

    auto gradAtFace(const mesh::Face& face) const -> Vector3d override;
    auto gradAtCell(const mesh::Cell& cell) const -> Vector3d override;
    auto gradAtCellStored(const mesh::Cell& cell) const -> Vector3d override;

    template <typename Func>
    void updateInteriorFaces(Func func);

    template <typename Func>
    void updateFaces(Func func);

    template <typename Func>
    void updateCells(Func func);

    void setGradScheme(const SharedPtr<gradient::IGradient>& grad_scheme);
    auto clone() const -> GeneralScalar;

    auto operator[](std::size_t i) const -> double;
    auto operator[](std::size_t i) -> double&;

  protected:
    auto valueAtInteriorFace(const mesh::Face& face) const -> double;
    auto valueAtBoundaryFace(const mesh::Face& face) const -> double;

  private:
    void setGradScheme();
    void addDefaultBoundaryHandlers();

    SharedPtr<VectorXd> _data = nullptr;

    /// TODO: _face_data should not include empty faces
    SharedPtr<VectorXd> _face_data = nullptr;
    SharedPtr<gradient::IGradient> _grad_scheme = nullptr;

    const IVector* _parent = nullptr;
    std::optional<Coord> _coord = std::nullopt;
    BHManagerSetter _setter;
};

class ScalarBHManagerSetter {
  public:
    using IScalarBHManager =
        prism::boundary::BoundaryHandlersManager<boundary::scalar::IScalarBoundaryHandler>;

    static void set(IScalarBHManager& manager);
};

using Scalar = GeneralScalar<units::Measurable, ScalarBHManagerSetter>;

template <IScalarBased ScalarField>
auto UniformScalar::operator*(const ScalarField& other) -> ScalarField {
    /// TODO: this is not a good way to compare meshes, we should use a mesh equality operator.
    if (this->mesh()->cellCount() != other.mesh()->cellCount()) {
        throw std::runtime_error(
            fmt::format("UniformScalar::operator*(): cannot multiply scalar field '{}' with "
                        "scalar field '{}', because they are defined on different meshes.",
                        this->name(),
                        other.name()));
    }
    /// TODO: this ignores the units of the scalar field
    return ScalarField(fmt::format("{}{}", this->name(), other.name()),
                       this->mesh(),
                       this->_value * other.values());
}

template <IVectorBased VectorField>
auto UniformScalar::operator*(const VectorField& other) -> VectorField {
    /// TODO: this function ignores the units of the scalar field, which is not ideal, and returns
    /// a vector of components that preserves `other`'s units. We need to fix this.
    if (this->mesh()->cellCount() != other.mesh()->cellCount()) {
        throw std::runtime_error(
            fmt::format("UniformScalar::operator*(): cannot multiply scalar field '{}' with "
                        "vector field '{}', because they are defined on different meshes.",
                        this->name(),
                        other.name()));
    }

    auto output_name = fmt::format("{}{}", this->name(), other.name());
    std::vector<Vector3d> face_values(other.mesh()->faceCount(), Vector3d::Zero());
    VectorXd face_flux_values = VectorXd::Zero(other.mesh()->faceCount());

    for (const auto& face : this->mesh()->interiorFaces()) {
        const Vector3d face_value = other.valueAtFace(face) * this->_value;
        face_values[face.id()] = face_value;
        face_flux_values[face.id()] = face_value.dot(face.areaVector());
        ;
    }
    for (const auto& patch : this->mesh()->boundaryPatches()) {
        if (patch.isEmpty()) {
            continue; // skip empty patches
        }
        for (const auto& face_id : patch.facesIds()) {
            Vector3d result = other.valueAtFace(face_id) * this->_value;
            const auto& face = this->mesh()->face(face_id);
            face_values[face_id] = result;
            face_flux_values[face_id] = result.dot(face.areaVector());
        }
    }

    using Component = typename VectorField::ComponentType;
    std::array<Component, 3> components = {
        Component(output_name + "_x", this->mesh(), this->_value * other.x().values().array()),
        Component(output_name + "_y", this->mesh(), this->_value * other.y().values().array()),
        Component(output_name + "_z", this->mesh(), this->_value * other.z().values().array())};

    auto U =
        VectorField(fmt::format("{}{}", this->name(), other.name()), this->mesh(), components);
    U.setFaceValues(std::move(face_values));
    U.setFaceFluxValues(std::move(face_flux_values));
    return U;
}

template <typename Units, typename BHManagerSetter>
GeneralScalar<Units, BHManagerSetter>::GeneralScalar(std::string name,
                                                     const SharedPtr<mesh::PMesh>& mesh,
                                                     double value,
                                                     IVector* parent)
    : IScalar(std::move(name), mesh),
      _data(std::make_shared<VectorXd>(VectorXd::Ones(mesh->cellCount()) * value)),
      _parent(parent) {
    log::debug("Creating scalar field: '{}' with double value = {}", this->name(), value);
    addDefaultBoundaryHandlers();
    setGradScheme();
}

template <typename Units, typename BHManagerSetter>
GeneralScalar<Units, BHManagerSetter>::GeneralScalar(std::string name,
                                                     const SharedPtr<mesh::PMesh>& mesh,
                                                     double value,
                                                     Coord coord,
                                                     IVector* parent)
    : IScalar(std::move(name), mesh),
      _data(std::make_shared<VectorXd>(VectorXd::Ones(mesh->cellCount()) * value)),
      _coord(coord),
      _parent(parent) {
    log::debug("Creating scalar field: '{}' (as {}-coordinate) with double value = {}",
               this->name(),
               coordToStr(coord),
               value);

    addDefaultBoundaryHandlers();
    setGradScheme();
}

template <typename Units, typename BHManagerSetter>
GeneralScalar<Units, BHManagerSetter>::GeneralScalar(std::string name,
                                                     const SharedPtr<mesh::PMesh>& mesh,
                                                     VectorXd data,
                                                     IVector* parent)
    : IScalar(std::move(name), mesh),
      _data(std::make_shared<VectorXd>(std::move(data))),
      _parent(parent) {
    if (_data->size() != mesh->cellCount()) {
        throw std::runtime_error(
            fmt::format("field::GeneralScalar() cannot create a scalar field '{}' given a "
                        "vector that has a different size than mesh's cell count.",
                        this->name()));
    }

    log::debug("Creating scalar field: '{}' with a cell vector data of size = {}",
               this->name(),
               _data->size());

    addDefaultBoundaryHandlers();
    setGradScheme();
}

template <typename Units, typename BHManagerSetter>
GeneralScalar<Units, BHManagerSetter>::GeneralScalar(std::string name,
                                                     const SharedPtr<mesh::PMesh>& mesh,
                                                     VectorXd data,
                                                     Coord coord,
                                                     IVector* parent)
    : IScalar(std::move(name), mesh),
      _data(std::make_shared<VectorXd>(std::move(data))),
      _coord(coord),
      _parent(parent) {
    if (_data->size() != mesh->cellCount()) {
        throw std::runtime_error(fmt::format(
            "field::Scalar() cannot create a scalar field '{}' given a vector that has a "
            "different size than mesh's cell count.",
            this->name()));
    }

    log::debug(
        "Creating scalar field: '{}' (as {}-coordinate) with a cell vector data of size = {}",
        this->name(),
        coordToStr(coord),
        _data->size());
    addDefaultBoundaryHandlers();
    setGradScheme();
}

template <typename Units, typename BHManagerSetter>
GeneralScalar<Units, BHManagerSetter>::GeneralScalar(std::string name,
                                                     const SharedPtr<mesh::PMesh>& mesh,
                                                     VectorXd data,
                                                     VectorXd face_data,
                                                     IVector* parent)
    : IScalar(std::move(name), mesh),
      _data(std::make_shared<VectorXd>(std::move(data))),
      _face_data(std::make_shared<VectorXd>(std::move(face_data))),
      _parent(parent) {
    if (_data->size() != mesh->cellCount()) {
        throw std::runtime_error(fmt::format(
            "field::Scalar() cannot create a scalar field '{}' given a vector that has a "
            "different size than mesh's cell count.",
            this->name()));
    }

    if (_face_data->size() != mesh->faceCount()) {
        throw std::runtime_error(
            fmt::format("field::Scalar() cannot create a scalar field '{}' given a face data "
                        "vector that has a different size than mesh's faces count.",
                        this->name()));
    }

    log::debug(
        "Creating scalar field: '{}' with a cell data vector of size = {} and face data "
        "vector of size = {}",
        this->name(),
        _data->size(),
        _face_data->size());

    addDefaultBoundaryHandlers();
    setGradScheme();
}

template <typename Units, typename BHManagerSetter>
GeneralScalar<Units, BHManagerSetter>::GeneralScalar(std::string name,
                                                     const SharedPtr<mesh::PMesh>& mesh,
                                                     VectorXd data,
                                                     VectorXd face_data,
                                                     Coord coord,
                                                     IVector* parent)
    : IScalar(std::move(name), mesh),
      _data(std::make_shared<VectorXd>(std::move(data))),
      _face_data(std::make_shared<VectorXd>(std::move(face_data))),
      _coord(coord),
      _parent(parent) {
    if (_data->size() != mesh->cellCount()) {
        throw std::runtime_error(fmt::format(
            "field::Scalar() cannot create a scalar field '{}' given a vector that has a "
            "different size than mesh's cell count.",
            this->name()));
    }

    if (_face_data->size() != mesh->faceCount()) {
        throw std::runtime_error(
            fmt::format("field::Scalar() cannot create a scalar field '{}' given a face data "
                        "vector that has a different size than mesh's faces count.",
                        this->name()));
    }

    log::debug(
        "Creating scalar field: '{}' (as {}-coordinate) with a cell data vector of size = {} "
        "and "
        "face data vector of size = {}",
        this->name(),
        coordToStr(coord),
        _data->size(),
        _face_data->size());

    addDefaultBoundaryHandlers();
    setGradScheme();
}

/// TODO: delete this in favor of updateFaceValues()
template <typename Units, typename BHManagerSetter>
void GeneralScalar<Units, BHManagerSetter>::setFaceValues(VectorXd values) {
    if (values.size() != mesh()->faceCount()) {
        throw std::runtime_error(
            fmt::format("prism::field::GeneralScalar<Units, BHManagerProvider, "
                        "BHManagerSetter>::setFaceValues(): cannot set face values for scalar "
                        "field {}, to a face data vector having a different size that "
                        "field's faces count.",
                        name()));
    }

    if (hasFaceValues()) {
        log::debug(
            "GeneralScalar<Units, BHManagerSetter>::setFaceValues() was called for field "
            "'{}', which already has face values set.",
            name());
    }

    log::debug("Setting face values for field '{}'", name());
    _face_data = std::make_shared<VectorXd>(std::move(values));
}

template <typename Units, typename BHManagerSetter>
void GeneralScalar<Units, BHManagerSetter>::clearFaceValues() {
    _face_data = nullptr;
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::values() const -> const VectorXd& {
    if (_data == nullptr) {
        throw std::runtime_error(
            fmt::format("prism::field::GeneralScalar<Units, BHManagerSetter>::values() was "
                        "called for field `{}`, but the data is not initialized.",
                        name()));
    }
    return *_data;
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::values() -> VectorXd& {
    if (_data == nullptr) {
        throw std::runtime_error(
            fmt::format("prism::field::GeneralScalar<Units, BHManagerSetter>::values() was "
                        "called for field `{}`, but the data is not initialized.",
                        name()));
    }
    return *_data;
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::coord() const noexcept -> std::optional<Coord> {
    return _coord;
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::hasFaceValues() const -> bool {
    return _face_data != nullptr;
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::valueAtCell(const mesh::Cell& cell) const -> double {
    return valueAtCell(cell.id());
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::valueAtCell(std::size_t cell_id) const -> double {
    assert(_data != nullptr);              // NOLINT
    assert(cell_id < mesh()->cellCount()); // NOLINT
    return (*_data)[cell_id];
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::valueAtFace(std::size_t face_id) const -> double {
    if (hasFaceValues()) {
        // Face data were calculated already, just return the value (as in Rhie-Chow
        // correction).
        return (*_face_data)[face_id];
    }

    const auto& face = mesh()->face(face_id);

    if (face.isInterior()) {
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
    assert(face.isInterior()); // NOLINT
    const auto& owner = mesh()->cell(face.owner());
    const auto& neighbor = mesh()->cell(face.neighbor().value());

    const auto gc = mesh::geometricWeight(owner, neighbor, face);
    double val = gc * (*_data)[owner.id()];
    val += (1 - gc) * (*_data)[neighbor.id()];

    return val;
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::valueAtBoundaryFace(const mesh::Face& face) const
    -> double {
    const auto& patch = mesh()->boundaryPatch(face);
    const auto& bc = patch.getBoundaryCondition(name());

    auto handler = this->boundaryHandlersManager().getHandler(bc.kindString());

    if (handler == nullptr) {
        throw error::NonImplementedBoundaryCondition(
            fmt::format("prism::field::GeneralScalar<Units, BHManagerProvider, "
                        "BHManagerSetter>::valueAtBoundaryFace() for field `{}`",
                        name()),
            patch.name(),
            bc.kindString());
    }

    return handler->get(*this, face);
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::parent() const -> const IVector* {
    return _parent;
}

template <typename Units, typename BHManagerSetter>
void GeneralScalar<Units, BHManagerSetter>::setParent(const IVector* parent) {
    _parent = parent;
    /// TODO: check parent and component names consistency
}

template <typename Units, typename BHManagerSetter>
void GeneralScalar<Units, BHManagerSetter>::addDefaultBoundaryHandlers() {
    _setter.set(this->boundaryHandlersManager());
}

template <typename Units, typename BHManagerSetter>
void GeneralScalar<Units, BHManagerSetter>::setGradScheme(
    const SharedPtr<gradient::IGradient>& grad_scheme) {
    if (grad_scheme == nullptr) {
        throw std::runtime_error(
            "GeneralScalar::setGradScheme() was given a null gradient scheme pointer");
    }
    _grad_scheme = grad_scheme;
}

template <typename Units, typename BHManagerSetter>
void GeneralScalar<Units, BHManagerSetter>::setGradScheme() {
    // did user specify gradient scheme for the field in `fields.json`?
    auto field_infos = this->mesh()->fieldsInfo();
    auto it = std::find_if(field_infos.begin(), field_infos.end(), [this](const auto& fi) {
        return fi.name() == this->name() && fi.gradScheme().has_value();
    });

    /// TODO: this is buggy, it doesn't find the grad scheme defined in fields.json for the
    // field, also does not consider vector fields.
    if (it == field_infos.end()) {
        log::debug(
            "GeneralScalar::setGradScheme(): couldn't find a specified gradient scheme for "
            "field `{}` in `fields.json`, setting the gradient scheme to least squares.",
            this->name());

        _grad_scheme = std::make_shared<gradient::LeastSquares>(this);
        return;
    }

    auto grad_scheme_name = it->gradScheme().value();

    if (grad_scheme_name == "green-gauss" || grad_scheme_name == "greenGauss") {
        log::debug(
            "GeneralScalar::setGradScheme(): setting the gradient scheme to Green-Gauss for "
            "field `{}`",
            this->name());
        _grad_scheme = std::make_shared<gradient::GreenGauss>(this);
        return;
    }

    if (grad_scheme_name == "least-squares" || grad_scheme_name == "leastSquares") {
        log::debug(
            "GeneralScalar::setGradScheme(): setting the gradient scheme to Least-Squares "
            "for field `{}`",
            this->name());
        _grad_scheme = std::make_shared<gradient::LeastSquares>(this);
        return;
    }
    _grad_scheme = std::make_shared<gradient::LeastSquares>(this);
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::gradAtFace(const mesh::Face& face) const -> Vector3d {
    return _grad_scheme->gradAtFace(face);
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::gradAtCell(const mesh::Cell& cell) const -> Vector3d {
    return _grad_scheme->gradAtCell(cell);
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::gradAtCellStored(const mesh::Cell& cell) const
    -> Vector3d {
    return _grad_scheme->gradAtCellStored(cell);
}

template <typename Units, typename BHManagerSetter>
template <typename Func>
void GeneralScalar<Units, BHManagerSetter>::updateInteriorFaces(Func func) {
    if (!hasFaceValues()) {
        /// TODO: implement initFaceDataVector() to initialize face data vector
        VectorXd face_values(this->mesh()->faceCount());
        face_values.setZero();

        for (const auto& patch : this->mesh()->boundaryPatches()) {
            if (patch.isEmpty()) {
                continue; // skip empty patches
            }
            for (const auto& face_id : patch.facesIds()) {
                face_values[face_id] = valueAtFace(face_id);
            }
        }
        _face_data = std::make_shared<VectorXd>(std::move(face_values));
    }

    for (const auto& face : this->mesh()->interiorFaces()) {
        (*_face_data)[face.id()] = func(face);
    }
}

template <typename Units, typename BHManagerSetter>
template <typename Func>
void GeneralScalar<Units, BHManagerSetter>::updateFaces(Func func) {
    updateInteriorFaces(func);

    for (const auto& patch : this->mesh()->boundaryPatches()) {
        if (patch.isEmpty()) {
            continue; // skip empty patches
        }
        for (const auto& face_id : patch.facesIds()) {
            const auto& face = this->mesh()->face(face_id);
            (*_face_data)[face_id] = func(face);
        }
    }
}

template <typename Units, typename BHManagerSetter>
template <typename Func>
void GeneralScalar<Units, BHManagerSetter>::updateCells(Func func) {
    /// TODO: test this.
    for (const auto& cell : this->mesh()->cells()) {
        (*_data)[cell.id()] = func(cell);
    }
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::clone() const -> GeneralScalar {
    /// NOTE: cloned field is parentless.
    if (this->coord().has_value()) {
        auto clone =
            GeneralScalar(this->name(), this->mesh(), this->values(), this->coord().value());
        if (hasFaceValues()) {
            clone.setFaceValues(*_face_data);
        }
        return clone;
    }
    auto clone = GeneralScalar(this->name(), this->mesh(), this->values());
    if (hasFaceValues()) {
        clone.setFaceValues(*_face_data);
    }
    return clone;
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::operator[](std::size_t i) const -> double {
    if (_data == nullptr) {
        throw std::runtime_error(
            fmt::format("prism::field::GeneralScalar<Units, BHManagerSetter>::operator[]() was "
                        "called for field `{}`, but the data is not initialized.",
                        name()));
    }
    if (i >= _data->size()) {
        throw std::out_of_range(
            fmt::format("prism::field::GeneralScalar<Units, BHManagerSetter>::operator[]() was "
                        "called for field `{}`, but the index {} is out of range (size = {}).",
                        name(),
                        i,
                        _data->size()));
    }
    return (*_data)[i];
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::operator[](std::size_t i) -> double& {
    if (_data == nullptr) {
        throw std::runtime_error(
            fmt::format("prism::field::GeneralScalar<Units, BHManagerSetter>::operator[]() was "
                        "called for field `{}`, but the data is not initialized.",
                        name()));
    }

    if (i >= _data->size()) {
        throw std::out_of_range(
            fmt::format("prism::field::GeneralScalar<Units, BHManagerSetter>::operator[]() was "
                        "called for field `{}`, but the index {} is out of range (size = {}).",
                        name(),
                        i,
                        _data->size()));
    }
    return (*_data)[i];
}

} // namespace prism::field
