#pragma once

#include <array>
#include <memory>
#include <string>

#include "mesh/pmesh.h"
#include "prism/mesh/cell.h"
#include "prism/mesh/face.h"
#include "types.h"

namespace prism {

class Field {
  public:
    Field(std::string name, const mesh::PMesh& mesh);

    Field(const Field& other) = default;
    Field(Field&& other) noexcept = default;
    auto operator=(const Field& other) -> Field& = default;
    auto operator=(Field&& other) noexcept -> Field& = default;
    virtual ~Field() = default;

    auto inline name() const -> const std::string& { return _name; }
    auto inline name() -> std::string& { return _name; }

    auto inline mesh() const -> const mesh::PMesh& { return *_mesh; }

    // TODO: rename this to has_boundary_values()
    virtual auto has_face_data() const -> bool { return false; }

  private:
    const mesh::PMesh* _mesh = nullptr;
    std::string _name;
};

class ScalarField : public Field {
  public:
    ScalarField(std::string name, const mesh::PMesh& mesh, double value);
    ScalarField(std::string name, const mesh::PMesh& mesh, VectorXd data);
    ScalarField(std::string name, const mesh::PMesh& mesh, VectorXd data, VectorXd face_data);

    auto inline data() const -> const VectorXd& { return *_data; }
    auto inline data() -> VectorXd& { return *_data; }

    auto inline has_face_data() const -> bool override { return _face_data != nullptr; }
    void set_face_values(VectorXd values);

    auto value_at_cell(std::size_t cell_id) const -> double;
    auto value_at_cell(const mesh::Cell& cell) const -> double;

    auto value_at_face(std::size_t face_id) const -> double;
    auto value_at_face(const mesh::Face& face) const -> double;

    auto clone() const -> ScalarField;

    using CoordinatesMapper = double(double, double, double);
    using CellMapper = double(const mesh::Cell&);
    auto map(CellMapper* mapper) -> ScalarField&;
    auto map(CoordinatesMapper* mapper) -> ScalarField&;

    auto inline operator[](std::size_t i) const -> double { return (*_data)[i]; }
    auto inline operator[](std::size_t i) -> double& { return (*_data)[i]; }

  protected:
    auto value_at_interior_face(const mesh::Face& face) const -> double;
    auto value_at_boundary_face(const mesh::Face& face) const -> double;

  private:
    std::shared_ptr<VectorXd> _data = nullptr;
    std::shared_ptr<VectorXd> _face_data = nullptr;
};

class VectorField : public Field {
  public:
    VectorField(std::string name, const mesh::PMesh& mesh, double value);
    VectorField(std::string name, const mesh::PMesh& mesh, const Vector3d& data);
    VectorField(std::string name,
                const mesh::PMesh& mesh,
                const std::array<ScalarField, 3>& fields);

    auto has_face_data() const -> bool override;

    auto value_at_cell(std::size_t cell_id) const -> Vector3d;
    auto value_at_cell(const mesh::Cell& cell) const -> Vector3d;

    auto value_at_face(std::size_t face_id) const -> Vector3d;
    auto value_at_face(const mesh::Face& face) const -> Vector3d;

    auto inline x() -> ScalarField& { return _x; }
    auto inline y() -> ScalarField& { return _y; }
    auto inline z() -> ScalarField& { return _z; }

    auto operator[](std::size_t i) const -> Vector3d;

  private:
    ScalarField _x, _y, _z;
};

class TensorField : public Field {
  public:
    TensorField(std::string name, const mesh::PMesh& mesh, double value);
    TensorField(std::string name, const mesh::PMesh& mesh, Matrix3d data);
    TensorField(std::string name, const mesh::PMesh& mesh, std::vector<Matrix3d> data);

    auto value_at_cell(std::size_t cell_id) const -> const Matrix3d&;
    auto value_at_cell(const mesh::Cell& cell) const -> const Matrix3d&;

    auto value_at_face(std::size_t face_id) const -> Matrix3d;
    auto value_at_face(const mesh::Face& face) const -> Matrix3d;

    auto at(std::size_t i, std::size_t j, std::size_t k) -> double&;
    auto at(std::size_t i, std::size_t j, std::size_t k) const -> double;

    auto operator[](std::size_t i) -> Matrix3d&;
    auto operator[](std::size_t i) const -> const Matrix3d&;

  private:
    void init_data_vec();
    std::vector<Matrix3d> _data;
};

class PressureField : public ScalarField {
  public:
    PressureField(std::string name, const mesh::PMesh& mesh, double value);
    PressureField(std::string name, const mesh::PMesh& mesh, VectorXd data);
    PressureField(std::string name, const mesh::PMesh& mesh, VectorXd data, VectorXd face_data);
};

class VelocityField : public ScalarField {};

class VelocityVectorField : public VectorField {};


} // namespace prism