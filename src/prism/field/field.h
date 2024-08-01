#pragma once

#include <array>
#include <memory>
#include <string>

#include "prism/mesh/cell.h"
#include "prism/mesh/face.h"
#include "prism/mesh/pmesh.h"
#include "prism/types.h"

namespace prism::field {

namespace detail {
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

    virtual auto has_face_data() const -> bool { return false; }

    virtual auto value_at_cell(std::size_t cell_id) const -> CellValueType = 0;
    virtual auto value_at_cell(const mesh::Cell& cell) const -> CellValueType = 0;

    virtual auto value_at_face(std::size_t face_id) const -> CellValueType = 0;
    virtual auto value_at_face(const mesh::Face& face) const -> CellValueType = 0;

  private:
    const mesh::PMesh* _mesh = nullptr;
    std::string _name;
};

template <typename CellValueType>
IField<CellValueType>::IField(std::string name, const mesh::PMesh& mesh)
    : _name(std::move(name)), _mesh(&mesh) {
    detail::check_field_name(_name);
    detail::check_mesh(mesh);
}


class UniformScalar : public IField<double> {
  public:
    UniformScalar(std::string name, const mesh::PMesh& mesh, double value);

    auto value_at_cell(std::size_t cell_id) const -> double override;
    auto value_at_cell(const mesh::Cell& cell) const -> double override;

    auto value_at_face(std::size_t face_id) const -> double override;
    auto value_at_face(const mesh::Face& face) const -> double override;

  private:
    double _value {0.0};
};

class Scalar : public IField<double> {
  public:
    Scalar(std::string name, const mesh::PMesh& mesh, double value);
    Scalar(std::string name, const mesh::PMesh& mesh, VectorXd data);
    Scalar(std::string name, const mesh::PMesh& mesh, VectorXd data, VectorXd face_data);

    // TODO: check that _data is not null before returning, and maybe wrap it in an optional type
    auto inline data() const -> const VectorXd& { return *_data; }
    auto inline data() -> VectorXd& { return *_data; }

    auto inline has_face_data() const -> bool override { return _face_data != nullptr; }
    void set_face_values(VectorXd values);

    auto value_at_cell(std::size_t cell_id) const -> double override;
    auto value_at_cell(const mesh::Cell& cell) const -> double override;

    auto value_at_face(std::size_t face_id) const -> double override;
    auto value_at_face(const mesh::Face& face) const -> double override;

    auto inline operator[](std::size_t i) const -> double { return (*_data)[i]; }
    auto inline operator[](std::size_t i) -> double& { return (*_data)[i]; }

  protected:
    auto value_at_interior_face(const mesh::Face& face) const -> double;
    auto value_at_boundary_face(const mesh::Face& face) const -> double;

  private:
    std::shared_ptr<VectorXd> _data = nullptr;
    std::shared_ptr<VectorXd> _face_data = nullptr;
};

class Vector : public IField<Vector3d> {
  public:
    Vector(std::string name, const mesh::PMesh& mesh, double value);
    Vector(std::string name, const mesh::PMesh& mesh, const Vector3d& data);
    Vector(std::string name, const mesh::PMesh& mesh, const std::array<Scalar, 3>& fields);

    auto has_face_data() const -> bool override;

    auto value_at_cell(std::size_t cell_id) const -> Vector3d override;
    auto value_at_cell(const mesh::Cell& cell) const -> Vector3d override;

    auto value_at_face(std::size_t face_id) const -> Vector3d override;
    auto value_at_face(const mesh::Face& face) const -> Vector3d override;

    auto inline x() -> Scalar& { return _x; }
    auto inline y() -> Scalar& { return _y; }
    auto inline z() -> Scalar& { return _z; }

    auto operator[](std::size_t i) const -> Vector3d;

  private:
    Scalar _x, _y, _z;
};

class Tensor : public IField<Matrix3d> {
  public:
    Tensor(std::string name, const mesh::PMesh& mesh, double value);
    Tensor(std::string name, const mesh::PMesh& mesh, const Matrix3d& data);
    Tensor(std::string name, const mesh::PMesh& mesh, std::vector<Matrix3d> data);

    auto value_at_cell(std::size_t cell_id) const -> Matrix3d override;
    auto value_at_cell(const mesh::Cell& cell) const -> Matrix3d override;

    auto value_at_face(std::size_t face_id) const -> Matrix3d override;
    auto value_at_face(const mesh::Face& face) const -> Matrix3d override;

    auto operator[](std::size_t i) -> Matrix3d&;
    auto operator[](std::size_t i) const -> const Matrix3d&;

  private:
    std::vector<Matrix3d> _data;
};

class Pressure : public Scalar {
  public:
    Pressure(std::string name, const mesh::PMesh& mesh, double value);
    Pressure(std::string name, const mesh::PMesh& mesh, VectorXd data);
    Pressure(std::string name, const mesh::PMesh& mesh, VectorXd data, VectorXd face_data);
};


} // namespace prism::field