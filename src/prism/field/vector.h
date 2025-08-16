#pragma once

#include <fmt/format.h>

#include <array>

#include "ifield.h"
#include "prism/boundary.h"
#include "prism/log.h"
#include "prism/mesh/face.h"
#include "scalar.h"
#include "units.h"
#include "vector_boundary.h"

namespace prism::field {

template <typename Component, typename BHManagerSetter>
class GeneralVector
    : public IVector,
      public units::Measurable,
      public prism::boundary::BHManagerProvider<boundary::vector::IVectorBoundaryHandler> {
  public:
    GeneralVector(const std::string& name, const SharedPtr<mesh::PMesh>& mesh, double value);

    GeneralVector(const std::string& name,
                  const SharedPtr<mesh::PMesh>& mesh,
                  const Vector3d& data);

    /// TODO: replace with std::array<Component, 3> and use std::move in the constructor
    GeneralVector(std::string name,
                  const SharedPtr<mesh::PMesh>& mesh,
                  std::array<Component, 3>& components);

    /// TODO: impelement this.
    GeneralVector(std::string name,
                  const SharedPtr<mesh::PMesh>& mesh,
                  std::vector<Vector3d>& data);

    auto hasFaceValues() const -> bool override;
    auto hasFaceFluxValues() const -> bool;

    auto valueAtCell(std::size_t cell_id) const -> Vector3d override;
    auto valueAtCell(const mesh::Cell& cell) const -> Vector3d override;

    auto valueAtFace(std::size_t face_id) const -> Vector3d override;
    auto valueAtFace(const mesh::Face& face) const -> Vector3d override;

    auto fluxAtFace(std::size_t face_id) const -> double;
    auto fluxAtFace(const mesh::Face& face) const -> double;

    auto inline x() -> Component& { return _x; }
    auto inline y() -> Component& { return _y; }
    auto inline z() -> Component& { return _z; }
    auto inline x() const -> const Component& { return _x; }
    auto inline y() const -> const Component& { return _y; }
    auto inline z() const -> const Component& { return _z; }

    void setFaceValues(std::vector<Vector3d> values);
    void setFaceFluxValues(VectorXd values);

    /// TODO: updateInteriorFaces, updateFaces, updateCells and clone should all be virtual
    /// functions. Update IField.
    template <typename Func>
    void updateInteriorFaces(Func func);

    template <typename Func>
    void updateFaces(Func func);

    template <typename Func>
    void updateCells(Func func);

    auto clone() const -> GeneralVector;

    auto operator[](std::size_t i) const -> Vector3d;

    using ComponentType = Component;

  private:
    void initFaceData();
    void initFaceFluxDataVector();
    auto valueAtBoundaryFace(const mesh::Face& face) const -> Vector3d;
    auto fluxAtInteriorFace(const mesh::Face& face) const -> double;
    auto fluxAtBoundaryFace(const mesh::Face& face) const -> double;
    void addDefaultBoundaryHandlers();

    Component _x, _y, _z;
    SharedPtr<std::vector<Vector3d>> m_face_data = nullptr;
    SharedPtr<VectorXd> m_face_flux_data = nullptr;
    BHManagerSetter m_setter;
};

class VectorBHManagerSetter {
  public:
    using IVectorBHManager =
        prism::boundary::BoundaryHandlersManager<boundary::vector::IVectorBoundaryHandler>;

    static void set(IVectorBHManager& manager);
};
using Vector = GeneralVector<Scalar, VectorBHManagerSetter>;

template <typename Component, typename BHManagerSetter>
GeneralVector<Component, BHManagerSetter>::GeneralVector(const std::string& name,
                                                         const SharedPtr<mesh::PMesh>& mesh,
                                                         double value)
    : IVector(name, mesh),
      _x(this->name() + "_x", mesh, value, Coord::X),
      _y(this->name() + "_y", mesh, value, Coord::Y),
      _z(this->name() + "_z", mesh, value, Coord::Z) {
    log::debug("Creating vector field: '{}' with real value = {}", this->name(), value);
    addDefaultBoundaryHandlers();
}

template <typename Component, typename BHManagerSetter>
GeneralVector<Component, BHManagerSetter>::GeneralVector(const std::string& name,
                                                         const SharedPtr<mesh::PMesh>& mesh,
                                                         const Vector3d& data)
    : IVector(name, mesh),
      _x(this->name() + "_x", mesh, data[0], Coord::X),
      _y(this->name() + "_y", mesh, data[1], Coord::Y),
      _z(this->name() + "_z", mesh, data[2], Coord::Z) {
    log::debug("Creating vector field: '{}' with uniform vector value [{}, {}, {}]",
               this->name(),
               data.x(),
               data.y(),
               data.z());
    addDefaultBoundaryHandlers();
}

template <typename Component, typename BHManagerSetter>
GeneralVector<Component, BHManagerSetter>::GeneralVector(std::string name,
                                                         const SharedPtr<mesh::PMesh>& mesh,
                                                         std::array<Component, 3>& components)
    : IVector(std::move(name), mesh), _x(components[0]), _y(components[1]), _z(components[2]) {
    log::debug("Creating vector field: '{}'", this->name());
    // check sub-fields naming consistency
    if ((_x.name() != (this->name() + "_x")) || (_y.name() != (this->name() + "_y")) ||
        (_z.name() != (this->name() + "_z"))) {
        throw std::runtime_error(
            fmt::format("All field::Vector components names should end with '_x', '_y' or '_z'. "
                        "field::Vector constructor for `{}` vector field was given the following "
                        "field::Scalar names: '{}', '{}', '{}",
                        this->name(),
                        _x.name(),
                        _y.name(),
                        _z.name()));
    }
    addDefaultBoundaryHandlers();
}

template <typename Component, typename BHManagerSetter>
auto GeneralVector<Component, BHManagerSetter>::valueAtCell(std::size_t cell_id) const
    -> Vector3d {
    return operator[](cell_id);
}

template <typename Component, typename BHManagerSetter>
auto GeneralVector<Component, BHManagerSetter>::valueAtCell(const mesh::Cell& cell) const
    -> Vector3d {
    return valueAtCell(cell.id());
}

template <typename Component, typename BHManagerSetter>
auto GeneralVector<Component, BHManagerSetter>::valueAtFace(std::size_t face_id) const
    -> Vector3d {
    return valueAtFace(mesh()->face(face_id));
}

template <typename Component, typename BHManagerSetter>
auto GeneralVector<Component, BHManagerSetter>::valueAtFace(const mesh::Face& face) const
    -> Vector3d {
    const auto id = face.id();
    if (hasFaceValues()) {
        return (*m_face_data)[id];
    }

    if (face.isBoundary()) {
        return valueAtBoundaryFace(face);
    }
    return {_x.valueAtFace(id), _y.valueAtFace(id), _z.valueAtFace(id)};
}

template <typename Component, typename BHManagerSetter>
auto GeneralVector<Component, BHManagerSetter>::valueAtBoundaryFace(const mesh::Face& face) const
    -> Vector3d {
    const auto& patch = mesh()->boundaryPatch(face);
    const auto& bc = patch.getBoundaryCondition(name());

    auto handler = this->boundaryHandlersManager().getHandler(bc.kindString());

    if (handler == nullptr) {
        throw error::NonImplementedBoundaryCondition(
            fmt::format("prism::field::GeneralVector<Units, BHManagerProvider, "
                        "BHManagerSetter>::valueAtBoundaryFace() for field `{}`",
                        name()),
            patch.name(),
            bc.kindString());
    }

    return handler->get(*this, face);
}

template <typename Component, typename BHManagerSetter>
auto GeneralVector<Component, BHManagerSetter>::hasFaceValues() const -> bool {
    return m_face_data != nullptr;
}

template <typename Component, typename BHManagerSetter>
auto GeneralVector<Component, BHManagerSetter>::hasFaceFluxValues() const -> bool {
    return m_face_flux_data != nullptr;
}

template <typename Component, typename BHManagerSetter>
auto GeneralVector<Component, BHManagerSetter>::fluxAtFace(std::size_t face_id) const -> double {
    if (hasFaceFluxValues()) {
        return (*m_face_flux_data)[face_id];
    }

    const auto& face = this->mesh()->face(face_id);
    if (face.isInterior()) {
        return fluxAtInteriorFace(face);
    }
    return fluxAtBoundaryFace(face);
}

template <typename Component, typename BHManagerSetter>
auto GeneralVector<Component, BHManagerSetter>::fluxAtFace(const mesh::Face& face) const
    -> double {
    return fluxAtFace(face.id());
}

template <typename Component, typename BHManagerSetter>
auto GeneralVector<Component, BHManagerSetter>::fluxAtInteriorFace(const mesh::Face& face) const
    -> double {
    return valueAtFace(face).dot(face.areaVector());
}

template <typename Component, typename BHManagerSetter>
auto GeneralVector<Component, BHManagerSetter>::fluxAtBoundaryFace(const mesh::Face& face) const
    -> double {
    const auto& patch = mesh()->boundaryPatch(face);
    const auto& bc = patch.getBoundaryCondition(name());

    auto handler = this->boundaryHandlersManager().getHandler(bc.kindString());

    if (handler == nullptr) {
        throw error::NonImplementedBoundaryCondition(
            fmt::format("prism::field::GeneralVector<Units, BHManagerProvider, "
                        "BHManagerSetter>::fluxAtBoundaryFace() for field `{}`",
                        name()),
            patch.name(),
            bc.kindString());
    }

    return handler->flux(*this, face);
}

template <typename Component, typename BHManagerSetter>
auto GeneralVector<Component, BHManagerSetter>::operator[](std::size_t i) const -> Vector3d {
    return {_x.valueAtCell(i), _y.valueAtCell(i), _z.valueAtCell(i)};
}

template <typename Component, typename BHManagerSetter>
void GeneralVector<Component, BHManagerSetter>::initFaceData() {
    std::vector<Vector3d> face_values(this->mesh()->faceCount(), Vector3d::Zero());
    VectorXd face_flux_values = VectorXd::Zero(this->mesh()->faceCount());

    /// TODO: why are we looping over interior faces and boundary patches separately when we are
    /// calling valueAtFace() anyways?
    for (const auto& face : this->mesh()->interiorFaces()) {
        face_values[face.id()] = valueAtFace(face);
        face_flux_values[face.id()] = fluxAtFace(face);
    }

    for (const auto& patch : this->mesh()->boundaryPatches()) {
        if (patch.isEmpty()) {
            continue; // Skip empty patches.
        }
        for (const auto& face_id : patch.facesIds()) {
            face_values[face_id] = valueAtFace(face_id);
            face_flux_values[face_id] = fluxAtFace(face_id);
        }
    }
    m_face_data = std::make_shared<std::vector<Vector3d>>(std::move(face_values));
    m_face_flux_data = std::make_shared<VectorXd>(std::move(face_flux_values));
}

template <typename Component, typename BHManagerSetter>
template <typename Func>
void GeneralVector<Component, BHManagerSetter>::updateInteriorFaces(Func func) {
    if (!hasFaceValues()) {
        initFaceData();
    }

    for (const auto& face : this->mesh()->interiorFaces()) {
        const Vector3d updated_value = func(face);
        (*m_face_data)[face.id()] = updated_value;
        (*m_face_flux_data)[face.id()] = updated_value.dot(face.areaVector());
    }
}

template <typename Component, typename BHManagerSetter>
template <typename Func>
void GeneralVector<Component, BHManagerSetter>::updateFaces(Func func) {
    updateInteriorFaces(func);

    for (const auto& patch : this->mesh()->boundaryPatches()) {
        if (patch.isEmpty()) {
            continue; // Skip empty patches.
        }
        /// TODO: this is not efficient, because updateInteriorFaces() already iterates over
        /// boundary faces to initialize them. We need to fix.
        for (const auto& face_id : patch.facesIds()) {
            const auto& face = this->mesh()->face(face_id);
            const Vector3d updated_value = func(face);
            (*m_face_data)[face_id] = updated_value;
            (*m_face_flux_data)[face_id] = updated_value.dot(face.areaVector());
        }
    }
}

template <typename Component, typename BHManagerSetter>
auto GeneralVector<Component, BHManagerSetter>::clone() const -> GeneralVector {
    auto components = std::array<Component, 3> {
        this->x().clone(),
        this->y().clone(),
        this->z().clone(),
    };

    auto clone = GeneralVector(this->name(), this->mesh(), components);

    if (hasFaceValues()) {
        clone.setFaceValues(*m_face_data);
    }

    if (hasFaceFluxValues()) {
        clone.setFaceFluxValues(*m_face_flux_data);
    }
    return clone;
}

template <typename Component, typename BHManagerSetter>
void GeneralVector<Component, BHManagerSetter>::setFaceValues(std::vector<Vector3d> values) {
    if (values.size() != mesh()->faceCount()) {
        throw std::runtime_error(fmt::format(
            "prism::field::GeneralVector::setFaceValues(): cannot set face values for "
            "vector field {}, to a face data vector having a different size that field's "
            "faces count.",
            name()));
    }

    if (hasFaceValues()) {
        log::debug(
            "GeneralVector::setFaceValues() was called for field '{}', which already has "
            "face values set.",
            name());
    }

    log::debug("Setting face values for field '{}'", name());
    m_face_data = std::make_shared<std::vector<Vector3d>>(std::move(values));
}

template <typename Component, typename BHManagerSetter>
void GeneralVector<Component, BHManagerSetter>::setFaceFluxValues(VectorXd values) {
    if (values.size() != mesh()->faceCount()) {
        throw std::runtime_error(fmt::format(
            "prism::field::GeneralVector::setFaceFluxValues(): cannot set face flux values for "
            "vector field {}, to a face flux data vector having a different size that field's "
            "faces count.",
            name()));
    }
    if (hasFaceFluxValues()) {
        log::debug(
            "GeneralVector::setFaceFluxValues() was called for field '{}', which already has "
            "face flux values set.",
            name());
    }
    log::debug("Setting face flux values for field '{}'", name());
    m_face_flux_data = std::make_shared<VectorXd>(std::move(values));
}

template <typename Component, typename BHManagerSetter>
template <typename Func>
void GeneralVector<Component, BHManagerSetter>::updateCells(Func func) {
    for (const auto& cell : this->mesh()->cells()) {
        Vector3d update = func(cell);
        this->x()[cell.id()] = update.x();
        this->y()[cell.id()] = update.y();
        this->z()[cell.id()] = update.z();
    }
}

template <typename Component, typename BHManagerSetter>
void GeneralVector<Component, BHManagerSetter>::addDefaultBoundaryHandlers() {
    m_setter.set(this->boundaryHandlersManager());
}

} // namespace prism::field
