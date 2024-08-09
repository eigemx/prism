#pragma once

#include <optional>
#include <string>

#include "prism/mesh/cell.h"
#include "prism/mesh/face.h"
#include "prism/mesh/pmesh.h"

namespace prism::field {

namespace detail {
void checkFieldName(const std::string& name);
void checkMesh(const mesh::PMesh& mesh);
} // namespace detail

template <typename CellValueType>
class IField {
  public:
    IField(std::string name, const mesh::PMesh& mesh);
    IField(const IField& other) = default;
    IField(IField&& other) noexcept = default;
    auto operator=(const IField& other) -> IField& = default;
    auto operator=(IField&& other) noexcept -> IField& = default;
    virtual ~IField() = default;

    auto inline name() const -> const std::string& { return _name; }
    auto inline name() -> std::string& { return _name; }

    auto inline mesh() const -> const mesh::PMesh& { return *_mesh; }

    virtual auto hasFaceValues() const -> bool { return false; }
    virtual auto valueAtCell(std::size_t cell_id) const -> CellValueType = 0;
    virtual auto valueAtCell(const mesh::Cell& cell) const -> CellValueType = 0;
    virtual auto valueAtFace(std::size_t face_id) const -> CellValueType = 0;
    virtual auto valueAtFace(const mesh::Face& face) const -> CellValueType = 0;


    using ValueType = CellValueType;

  private:
    const mesh::PMesh* _mesh = nullptr;
    std::string _name;
};

template <typename CellValueType>
IField<CellValueType>::IField(std::string name, const mesh::PMesh& mesh)
    : _name(std::move(name)), _mesh(&mesh) {
    detail::checkFieldName(_name);
    detail::checkMesh(mesh);
}

class IComponent {
  public:
    IComponent() = default;
    IComponent(const IComponent& other) = default;
    IComponent(IComponent&& other) noexcept = default;
    auto operator=(const IComponent& other) -> IComponent& = default;
    auto operator=(IComponent&& other) noexcept -> IComponent& = default;
    virtual ~IComponent() = default;

    virtual auto coord() const noexcept -> std::optional<Coord> = 0;
};

} // namespace prism::field