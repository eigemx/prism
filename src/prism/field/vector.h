#pragma once

#include <fmt/format.h>

#include <array>

#include "ifield.h"
#include "prism/log.h"
#include "scalar.h"
#include "units.h"

namespace prism::field {

template <typename Component>
class GeneralVector : public IField<Vector3d>, public IVector, public units::Measurable {
  public:
    GeneralVector(std::string name, const SharedPtr<mesh::PMesh>& mesh, double value);
    GeneralVector(std::string name, const SharedPtr<mesh::PMesh>& mesh, const Vector3d& data);

    /// TODO: replace with std::array<Component, 3> and use std::move in the constructor
    GeneralVector(std::string name,
                  const SharedPtr<mesh::PMesh>& mesh,
                  std::array<Component, 3>& fields);
    GeneralVector(std::string name,
                  const SharedPtr<mesh::PMesh>& mesh,
                  std::vector<Vector3d>& data);

    auto hasFaceValues() const -> bool override;

    auto valueAtCell(std::size_t cell_id) const -> Vector3d override;
    auto valueAtCell(const mesh::Cell& cell) const -> Vector3d override;

    auto valueAtFace(std::size_t face_id) const -> Vector3d override;
    auto valueAtFace(const mesh::Face& face) const -> Vector3d override;

    auto inline x() -> Component& { return _x; }
    auto inline y() -> Component& { return _y; }
    auto inline z() -> Component& { return _z; }
    auto inline x() const -> const Component& { return _x; }
    auto inline y() const -> const Component& { return _y; }
    auto inline z() const -> const Component& { return _z; }

    auto operator[](std::size_t i) const -> Vector3d;

    using ComponentType = Component;

  private:
    Component _x, _y, _z;
};

using Vector = GeneralVector<Scalar>;

template <typename ComponentType>
GeneralVector<ComponentType>::GeneralVector(std::string name,
                                            const SharedPtr<mesh::PMesh>& mesh,
                                            double value)
    : IField(std::move(name), mesh),
      _x(this->name() + "_x", mesh, value, Coord::X, static_cast<IVector*>(this)),
      _y(this->name() + "_y", mesh, value, Coord::Y, static_cast<IVector*>(this)),
      _z(this->name() + "_z", mesh, value, Coord::Z, static_cast<IVector*>(this)) {
    log::debug("Creating vector field: '{}' with double value = {}", this->name(), value);
}

template <typename ComponentType>
GeneralVector<ComponentType>::GeneralVector(std::string name,
                                            const SharedPtr<mesh::PMesh>& mesh,
                                            const Vector3d& data)
    : IField(std::move(name), mesh),
      _x(this->name() + "_x", mesh, data[0], Coord::X, static_cast<IVector*>(this)),
      _y(this->name() + "_y", mesh, data[1], Coord::Y, static_cast<IVector*>(this)),
      _z(this->name() + "_z", mesh, data[2], Coord::Z, static_cast<IVector*>(this)) {
    log::debug("Creating vector field: '{}' with uniform vector value", this->name());
}

template <typename ComponentType>
GeneralVector<ComponentType>::GeneralVector(std::string name,
                                            const SharedPtr<mesh::PMesh>& mesh,
                                            std::array<ComponentType, 3>& fields)
    : IField(std::move(name), mesh), _x(fields[0]), _y(fields[1]), _z(fields[2]) {
    log::debug("Creating vector field: '{}'", this->name());
    // check mesh consistency
    for (auto& field : fields) {
        if (&mesh != &field.mesh()) {
            throw std::runtime_error(
                fmt::format("field::Vector constructor was given a sub-field component with name "
                            "`{}` that is defined over a different mesh",
                            field.name()));
        }

        if (field.parent()) {
            log::warn(
                "field::Vector '{}' constructor was given a sub-field '{}' that already has a "
                "parent "
                "field::Vector",
                this->name(),
                field.name());
        }
        field.setParent(this);
    }

    // check sub-fields naming consistency
    if ((_x.name() != (this->name() + "_x")) || (_y.name() != (this->name() + "_y")) ||
        (_z.name() != (this->name() + "_z"))) {
        throw std::runtime_error(
            fmt::format("All field::Vector components names should end with '_x', '_y' or '_z'. "
                        "field::Vector constructor for `{}` vector field was given the following "
                        "field::Scalar names: '{}', "
                        "'{}', '{}",
                        this->name(),
                        _x.name(),
                        _y.name(),
                        _z.name()));
    }
}

template <typename ComponentType>
auto GeneralVector<ComponentType>::valueAtCell(std::size_t cell_id) const -> Vector3d {
    return operator[](cell_id);
}

template <typename ComponentType>
auto GeneralVector<ComponentType>::valueAtCell(const mesh::Cell& cell) const -> Vector3d {
    return valueAtCell(cell.id());
}

template <typename ComponentType>
auto GeneralVector<ComponentType>::valueAtFace(std::size_t face_id) const -> Vector3d {
    return {_x.valueAtFace(face_id), _y.valueAtFace(face_id), _z.valueAtFace(face_id)};
}

template <typename ComponentType>
auto GeneralVector<ComponentType>::valueAtFace(const mesh::Face& face) const -> Vector3d {
    return valueAtFace(face.id());
}

template <typename ComponentType>
auto GeneralVector<ComponentType>::hasFaceValues() const -> bool {
    return _x.hasFaceValues() && _y.hasFaceValues() && _z.hasFaceValues();
}

template <typename ComponentType>
auto GeneralVector<ComponentType>::operator[](std::size_t i) const -> Vector3d {
    return {_x.values()[i], _y.values()[i], _z.values()[i]};
}

} // namespace prism::field
