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
    ~ScalarField() noexcept = default;

    auto inline name() const -> const std::string& { return _name; }
    auto inline name() -> std::string& { return _name; }

    auto inline data() const -> const VectorXd& { return *_data; }
    auto inline data() -> VectorXd& { return *_data; }

    auto inline has_face_data() const -> bool { return _face_data != nullptr; }

    auto value_at_cell(std::size_t cell_id) const -> double;
    auto value_at_cell(const mesh::Cell& cell) const -> double;
    auto value_at_face(std::size_t face_id) const -> double;
    auto value_at_face(const mesh::Face& face) const -> double;

    auto inline mesh() const -> const mesh::PMesh& { return *_mesh; }

    auto clone() const -> ScalarField;

    using CoordinatesMapper = double(double, double, double);
    using CellMapper = double(const mesh::Cell&);

    auto map(CellMapper* mapper) -> ScalarField&;
    auto map(CoordinatesMapper* mapper) -> ScalarField&;

    auto inline operator[](std::size_t i) const -> const double& { return (*_data)[i]; }
    auto inline operator[](std::size_t i) -> double& { return (*_data)[i]; }

  private:
    auto value_at_interior_face(const mesh::Face& face) const -> double;
    auto value_at_boundary_face(const mesh::Face& face) const -> double;

    const mesh::PMesh* _mesh = nullptr;
    std::string _name;
    std::shared_ptr<VectorXd> _data = nullptr;
    std::shared_ptr<VectorXd> _face_data = nullptr;
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