#pragma once

#include <concepts>

#include "prism/mesh/cell.h"
#include "prism/mesh/face.h"
#include "prism/mesh/pmesh.h"
#include "prism/types.h"

namespace prism::gradient {
class IGradient;
}
namespace prism::field {

namespace detail {
void checkFieldName(const std::string& name);
void checkMesh(const SharedPtr<mesh::PMesh>& mesh);
} // namespace detail

template <typename CellValueType>
class IField {
  public:
    IField(std::string name, const SharedPtr<mesh::PMesh>& mesh);
    IField(const IField& other) = default;
    IField(IField&& other) noexcept = default;
    auto operator=(const IField& other) -> IField& = default;
    auto operator=(IField&& other) noexcept -> IField& = default;
    virtual ~IField() = default;

    auto name() const -> const std::string&;
    auto name() -> std::string&;
    auto mesh() const -> const SharedPtr<mesh::PMesh>&;

    virtual auto hasFaceValues() const -> bool;
    virtual auto valueAtCell(std::size_t cell_id) const -> CellValueType = 0;
    virtual auto valueAtCell(const mesh::Cell& cell) const -> CellValueType = 0;
    virtual auto valueAtFace(std::size_t face_id) const -> CellValueType = 0;
    virtual auto valueAtFace(const mesh::Face& face) const -> CellValueType = 0;
    virtual auto coord() const noexcept -> Optional<Coord>;

    using ValueType = CellValueType;

  private:
    SharedPtr<mesh::PMesh> _mesh = nullptr;
    String _name;
};

template <typename T>
concept IFieldBased = std::derived_from<T, IField<typename T::ValueType>>;

class IScalar : public IField<double> {
  public:
    IScalar(std::string name, const SharedPtr<mesh::PMesh>& mesh);

    /// TODO: implement gradAtX methods with face id and cell id.
    virtual auto gradAtFace(const mesh::Face& face) const -> Vector3d = 0;
    virtual auto gradAtCell(const mesh::Cell& cell) const -> Vector3d = 0;
    virtual auto gradAtCellStored(const mesh::Cell& cell) const -> Vector3d = 0;
};

template <typename T>
concept IScalarBased = std::derived_from<T, IScalar>;

class IVector : public IField<Vector3d> {
  public:
    IVector(std::string name, const SharedPtr<mesh::PMesh>& mesh);
};

template <typename T>
concept IVectorBased = std::derived_from<T, IVector>;


template <typename CellValueType>
IField<CellValueType>::IField(std::string name, const SharedPtr<mesh::PMesh>& mesh)
    : _name(std::move(name)), _mesh(mesh) {
    detail::checkFieldName(_name);
    detail::checkMesh(mesh);
}

template <typename CellValueType>
auto IField<CellValueType>::name() const -> const std::string& {
    return _name;
}

template <typename CellValueType>
auto IField<CellValueType>::name() -> std::string& {
    return _name;
}

template <typename CellValueType>
auto IField<CellValueType>::mesh() const -> const SharedPtr<mesh::PMesh>& {
    return _mesh;
}

template <typename CellValueType>
auto IField<CellValueType>::hasFaceValues() const -> bool {
    return false;
}

template <typename CellValueType>
auto IField<CellValueType>::coord() const noexcept -> Optional<Coord> {
    return NullOption;
}

/// TODO: This function should be moved to a more appropriate place, like a utility file.
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
        default: {
            return "";
        }
    }
}

} // namespace prism::field
