#pragma once

#include <memory>
#include <string>

#include "mesh/pmesh.h"
#include "print.h"
#include "types.h"

namespace prism {

class ScalarField {
  public:
    ScalarField(const ScalarField& other) = default;
    ScalarField(ScalarField&& other) noexcept = default;
    auto operator=(const ScalarField& other) -> ScalarField& = default;
    auto operator=(ScalarField&& other) noexcept -> ScalarField& = default;
    ~ScalarField() noexcept = default;

    ScalarField(std::string name, const mesh::PMesh& mesh);
    ScalarField(std::string name, const mesh::PMesh& mesh, VectorXd data);
    ScalarField(std::string name, const mesh::PMesh& mesh, double value);
    ScalarField(std::string name, const mesh::PMesh& mesh, std::shared_ptr<VectorXd> data);

    auto inline name() const -> const std::string& { return _name; }
    auto inline name() -> std::string& { return _name; }

    auto inline data() const -> const VectorXd& { return *_data; }
    auto inline data() -> VectorXd& { return *_data; }

    auto inline mesh() const -> const mesh::PMesh& { return *_mesh; }

    auto clone() const -> ScalarField;

    using CoordinatesMapper = double(double, double, double);
    using CellMapper = double(const mesh::Cell&);

    auto map(CellMapper* mapper) -> ScalarField&;
    auto map(CoordinatesMapper* mapper) -> ScalarField&;

    auto inline operator[](std::size_t i) const -> const double& { return (*_data)[i]; }
    auto inline operator[](std::size_t i) -> double& { return (*_data)[i]; }
    auto inline operator()() -> VectorXd& { return *_data; }

  private:
    const mesh::PMesh* _mesh = nullptr;
    std::string _name;
    std::shared_ptr<VectorXd> _data = nullptr;
};

class VectorField {
  public:
    VectorField(std::string name, const mesh::PMesh& mesh);
    VectorField(std::string name, const mesh::PMesh& mesh, const Vector3d& data);
    VectorField(std::string name, const mesh::PMesh& mesh, const MatrixX3d& data);
    VectorField(std::string name, const mesh::PMesh& mesh, double value);

    auto inline name() const -> const std::string& { return _name; }
    auto inline mesh() const -> const mesh::PMesh& { return *_mesh; }
    auto data() const -> MatrixX3d;


    auto x() -> ScalarField;
    auto y() -> ScalarField;
    auto z() -> ScalarField;

    auto inline operator[](std::size_t i) const -> Vector3d {
        return Vector3d {(*_x)[i], (*_y)[i], (*_z)[i]};
    }

  private:
    const mesh::PMesh* _mesh;
    std::string _name;
    std::shared_ptr<VectorXd> _x = nullptr;
    std::shared_ptr<VectorXd> _y = nullptr;
    std::shared_ptr<VectorXd> _z = nullptr;
};

} // namespace prism