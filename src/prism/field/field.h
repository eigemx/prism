#pragma once

#include <array>
#include <memory>
#include <optional>
#include <string>

#include "prism/boundary.h"
#include "prism/field/boundary.h"
#include "prism/mesh/cell.h"
#include "prism/mesh/face.h"
#include "prism/mesh/pmesh.h"
#include "prism/types.h"
#include "spdlog/spdlog.h"

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


class UniformScalar : public IField<double> {
  public:
    UniformScalar(std::string name, const mesh::PMesh& mesh, double value);

    auto valueAtCell(std::size_t cell_id) const -> double override;
    auto valueAtCell(const mesh::Cell& cell) const -> double override;
    auto valueAtFace(std::size_t face_id) const -> double override;
    auto valueAtFace(const mesh::Face& face) const -> double override;

  private:
    double _value {0.0};
};

namespace detail {
template <typename ComponentType>
class Vector;
}

class Scalar;
using Vector = detail::Vector<Scalar>;

class Scalar : public IField<double> {
  public:
    Scalar(std::string name, const mesh::PMesh& mesh, double value, Vector* parent = nullptr);
    Scalar(std::string name, const mesh::PMesh& mesh, VectorXd data, Vector* parent = nullptr);
    Scalar(std::string name,
           const mesh::PMesh& mesh,
           VectorXd data,
           VectorXd face_data,
           Vector* parent = nullptr);

    // TODO: check that _data is not null before returning, and maybe wrap it in an optional type
    auto inline values() const -> const VectorXd& { return *_data; }
    auto inline values() -> VectorXd& { return *_data; }

    auto inline hasFaceValues() const -> bool override { return _face_data != nullptr; }
    void setFaceValues(VectorXd values);

    auto valueAtCell(std::size_t cell_id) const -> double override;
    auto valueAtCell(const mesh::Cell& cell) const -> double override;
    auto valueAtFace(std::size_t face_id) const -> double override;
    auto valueAtFace(const mesh::Face& face) const -> double override;

    auto parent() -> std::optional<Vector>;
    void setParent(Vector* parent);

    auto inline operator[](std::size_t i) const -> double { return (*_data)[i]; }
    auto inline operator[](std::size_t i) -> double& { return (*_data)[i]; }

    using BoundaryHandlersManager =
        prism::boundary::BoundaryHandlersManager<Scalar,
                                                 boundary::FieldBoundaryHandler<Scalar, double>>;
    auto boundaryHandlersManager() -> BoundaryHandlersManager& { return _bh_manager; }

  protected:
    auto valueAtInteriorFace(const mesh::Face& face) const -> double;
    auto valueAtBoundaryFace(const mesh::Face& face) const -> double;

  private:
    BoundaryHandlersManager _bh_manager;
    std::shared_ptr<VectorXd> _data = nullptr;
    std::shared_ptr<VectorXd> _face_data = nullptr;
    Vector* _parent = nullptr;

    void addDefaultHandlers();
};

template <typename ComponentType>
class IVector {
  public:
    IVector() = default;
    IVector(IVector&) = default;
    IVector(IVector&&) noexcept = default;
    auto operator=(const IVector&) -> IVector& = default;
    auto operator=(IVector&&) noexcept -> IVector& = default;
    virtual ~IVector() = default;

    virtual auto x() -> ComponentType& = 0;
    virtual auto y() -> ComponentType& = 0;
    virtual auto z() -> ComponentType& = 0;
};

namespace detail {
template <typename ComponentType>
class Vector : public IField<Vector3d>, public IVector<ComponentType> {
  public:
    Vector(std::string name, const mesh::PMesh& mesh, double value);
    Vector(std::string name, const mesh::PMesh& mesh, const Vector3d& data);
    Vector(std::string name, const mesh::PMesh& mesh, std::array<Scalar, 3>& fields);

    auto hasFaceValues() const -> bool override;

    auto valueAtCell(std::size_t cell_id) const -> Vector3d override;
    auto valueAtCell(const mesh::Cell& cell) const -> Vector3d override;

    auto valueAtFace(std::size_t face_id) const -> Vector3d override;
    auto valueAtFace(const mesh::Face& face) const -> Vector3d override;

    auto inline x() -> ComponentType& override { return _x; }
    auto inline y() -> ComponentType& override { return _y; }
    auto inline z() -> ComponentType& override { return _z; }

    auto operator[](std::size_t i) const -> Vector3d;

  private:
    ComponentType _x, _y, _z;
};
} // namespace detail

class Tensor : public IField<Matrix3d> {
  public:
    Tensor(std::string name, const mesh::PMesh& mesh, double value);
    Tensor(std::string name, const mesh::PMesh& mesh, const Matrix3d& data);
    Tensor(std::string name, const mesh::PMesh& mesh, std::vector<Matrix3d> data);

    auto valueAtCell(std::size_t cell_id) const -> Matrix3d override;
    auto valueAtCell(const mesh::Cell& cell) const -> Matrix3d override;

    auto valueAtFace(std::size_t face_id) const -> Matrix3d override;
    auto valueAtFace(const mesh::Face& face) const -> Matrix3d override;

    auto operator[](std::size_t i) -> Matrix3d&;
    auto operator[](std::size_t i) const -> const Matrix3d&;

  private:
    std::vector<Matrix3d> _data;
};

namespace detail {
template <typename ComponentType>
Vector<ComponentType>::Vector(std::string name, const mesh::PMesh& mesh, double value)
    : IField(std::move(name), mesh),
      _x(this->name() + "_x", mesh, value, this),
      _y(this->name() + "_y", mesh, value, this),
      _z(this->name() + "_z", mesh, value, this) {}

template <typename ComponentType>
Vector<ComponentType>::Vector(std::string name, const mesh::PMesh& mesh, const Vector3d& data)
    : IField(std::move(name), mesh),
      _x(this->name() + "_x", mesh, data[0], this),
      _y(this->name() + "_y", mesh, data[1], this),
      _z(this->name() + "_z", mesh, data[2], this) {}

template <typename ComponentType>
Vector<ComponentType>::Vector(std::string name,
                              const mesh::PMesh& mesh,
                              std::array<Scalar, 3>& fields)
    : IField(std::move(name), mesh), _x(fields[0]), _y(fields[1]), _z(fields[2]) {
    // check mesh consistency
    for (auto& field : fields) {
        if (&mesh != &field.mesh()) {
            throw std::runtime_error(fmt::format(
                "field::Vector constructor was given a field::Scalar component with name "
                "`{}` that is defined over a different mesh",
                field.name()));
        }

        if (field.parent().has_value()) {
            spdlog::warn(
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
auto Vector<ComponentType>::valueAtCell(std::size_t cell_id) const -> Vector3d {
    return operator[](cell_id);
}

template <typename ComponentType>
auto Vector<ComponentType>::valueAtCell(const mesh::Cell& cell) const -> Vector3d {
    return valueAtCell(cell.id());
}

template <typename ComponentType>
auto Vector<ComponentType>::valueAtFace(std::size_t face_id) const -> Vector3d {
    return {_x.valueAtFace(face_id), _y.valueAtFace(face_id), _z.valueAtFace(face_id)};
}

template <typename ComponentType>
auto Vector<ComponentType>::valueAtFace(const mesh::Face& face) const -> Vector3d {
    return valueAtFace(face.id());
}

template <typename ComponentType>
auto Vector<ComponentType>::hasFaceValues() const -> bool {
    return _x.hasFaceValues() && _y.hasFaceValues() && _z.hasFaceValues();
}

template <typename ComponentType>
auto Vector<ComponentType>::operator[](std::size_t i) const -> Vector3d {
    return {_x.values()[i], _y.values()[i], _z.values()[i]};
}
} // namespace detail

} // namespace prism::field