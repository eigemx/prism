#pragma once

#include <array>
#include <memory>
#include <string>

#include "mesh/pmesh.h"
#include "prism/mesh/cell.h"
#include "prism/mesh/face.h"
#include "types.h"

namespace prism {

class ScalarField {
  public:
    ScalarField(std::string name, const mesh::PMesh& mesh, double value);
    ScalarField(std::string name, const mesh::PMesh& mesh, VectorXd data);
    ScalarField(std::string name, const mesh::PMesh& mesh, VectorXd data, VectorXd face_data);

    ScalarField(const ScalarField& other) = default;
    ScalarField(ScalarField&& other) noexcept = default;
    auto operator=(const ScalarField& other) -> ScalarField& = default;
    auto operator=(ScalarField&& other) noexcept -> ScalarField& = default;
    virtual ~ScalarField() noexcept = default;

    auto inline name() const -> const std::string& { return _name; }
    auto inline name() -> std::string& { return _name; }

    virtual auto inline data() const -> const VectorXd& { return *_data; }
    auto inline data() -> VectorXd& { return *_data; }

    virtual auto inline has_face_data() const -> bool { return _face_data != nullptr; }

    virtual auto value_at_cell(std::size_t cell_id) const -> double;
    virtual auto value_at_cell(const mesh::Cell& cell) const -> double;
    virtual auto value_at_face(std::size_t face_id) const -> double;
    virtual auto value_at_face(const mesh::Face& face) const -> double;

    auto inline mesh() const -> const mesh::PMesh& { return *_mesh; }

    virtual auto clone() const -> ScalarField;

    //using CoordinatesMapper = double(double, double, double);
    //using CellMapper = double(const mesh::Cell&);
    //virtual auto map(CellMapper* mapper) -> ScalarField&;
    //virtual auto map(CoordinatesMapper* mapper) -> ScalarField&;

    virtual auto inline operator[](std::size_t i) const -> double { return (*_data)[i]; }
    auto inline operator[](std::size_t i) -> double& { return (*_data)[i]; }

  protected:
    ScalarField(std::string name, const mesh::PMesh& mesh);
    virtual auto value_at_interior_face(const mesh::Face& face) const -> double;
    virtual auto value_at_boundary_face(const mesh::Face& face) const -> double;

  private:
    const mesh::PMesh* _mesh = nullptr;
    std::string _name;
    std::shared_ptr<VectorXd> _data = nullptr;
    std::shared_ptr<VectorXd> _face_data = nullptr;
};

// Sometimes we need to define a constant field over the mesh, such as gravity
// or constant diffusion coefficient, it would be better if we can do this without
// the memory overhead of ScalarField (defining the value at each cell).
// ConstantScalarField solves this problem by overriding ScalarField getters,
// and stores just a double value `_value` instead of a vector.
class ConstScalarField : public ScalarField {
  public:
    ConstScalarField(std::string name, const mesh::PMesh& mesh, double value);

    ConstScalarField(const ConstScalarField& other) = default;
    ConstScalarField(ConstScalarField&& other) noexcept = default;
    auto operator=(const ConstScalarField& other) -> ConstScalarField& = default;
    auto operator=(ConstScalarField&& other) noexcept -> ConstScalarField& = default;
    ~ConstScalarField() noexcept = default; // NOLINT

    auto inline has_face_data() const -> bool override { return false; }

    auto value_at_cell(std::size_t cell_id) const -> double override;
    auto value_at_cell(const mesh::Cell& cell) const -> double override;

    auto clone() const -> ScalarField override;

    auto inline operator[](std::size_t i) const -> double override { return _value; } //NOLINT

  private:
    auto value_at_interior_face(const mesh::Face& face) const -> double override;

    double _value {0.0};
    VectorXd _data;
};

class VectorField {
  public:
    VectorField(std::string name, const mesh::PMesh& mesh, double value);
    VectorField(std::string name, const mesh::PMesh& mesh, const Vector3d& data);
    VectorField(std::string name,
                const mesh::PMesh& mesh,
                const std::array<ScalarField, 3>& fields);

    auto inline name() const -> const std::string& { return _name; }
    auto inline mesh() const -> const mesh::PMesh& { return *_mesh; }

    auto has_face_data() const -> bool;

    auto value_at_cell(std::size_t cell_id) const -> Vector3d;
    auto value_at_cell(const mesh::Cell& cell) const -> Vector3d;
    auto value_at_face(std::size_t face_id) const -> Vector3d;
    auto value_at_face(const mesh::Face& face) const -> Vector3d;

    auto inline x() -> ScalarField { return _x; }
    auto inline y() -> ScalarField { return _y; }
    auto inline z() -> ScalarField { return _z; }

    auto inline operator[](std::size_t i) const -> Vector3d {
        return {_x.data()[i], _y.data()[i], _z.data()[i]};
    }

  private:
    const mesh::PMesh* _mesh;
    std::string _name;
    ScalarField _x, _y, _z;
};

} // namespace prism