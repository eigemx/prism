#pragma once

#include <concepts>

#include "boundary.h"
#include "ifield.h"
#include "units.h"

namespace prism::field {

class IScalar {};

class UniformScalar : public IField<double>, public IScalar {
  public:
    UniformScalar(std::string name, const mesh::PMesh& mesh, double value);

    auto valueAtCell(std::size_t cell_id) const -> double override;
    auto valueAtCell(const mesh::Cell& cell) const -> double override;
    auto valueAtFace(std::size_t face_id) const -> double override;
    auto valueAtFace(const mesh::Face& face) const -> double override;

  private:
    double _value {0.0};
};

class IVector;

class Scalar
    : public IField<double>,
      public IScalar,
      public IComponent,
      public units::Measurable,
      public prism::boundary::BHManagersProvider<Scalar, boundary::FieldBoundaryHandler<Scalar>> {
  public:
    // Uniform double value constructors
    Scalar(std::string name, const mesh::PMesh& mesh, double value, IVector* parent = nullptr);
    Scalar(std::string name,
           const mesh::PMesh& mesh,
           double value,
           Coord coord,
           IVector* parent = nullptr);

    // VectorXd cell values constructors
    Scalar(std::string name, const mesh::PMesh& mesh, VectorXd data, IVector* parent = nullptr);
    Scalar(std::string name,
           const mesh::PMesh& mesh,
           VectorXd data,
           Coord coord,
           IVector* parent = nullptr);

    // VectorXd cell & face values constructors
    Scalar(std::string name,
           const mesh::PMesh& mesh,
           VectorXd data,
           VectorXd face_data,
           IVector* parent = nullptr);
    Scalar(std::string name,
           const mesh::PMesh& mesh,
           VectorXd data,
           VectorXd face_data,
           Coord coord,
           IVector* parent = nullptr);

    Scalar() = delete;
    Scalar(const Scalar&) = default;
    Scalar(Scalar&&) = default;
    auto operator=(const Scalar&) -> Scalar& = default;
    auto operator=(Scalar&&) -> Scalar& = default;
    ~Scalar() override = default;

    // TODO: check that _data is not null before returning, and maybe wrap it in an optional type
    auto inline values() const -> const VectorXd& { return *_data; }
    auto inline values() -> VectorXd& { return *_data; }

    auto inline coord() const noexcept -> std::optional<Coord> override { return _coord; }
    auto inline hasFaceValues() const -> bool override { return _face_data != nullptr; }
    void setFaceValues(VectorXd values);

    auto valueAtCell(std::size_t cell_id) const -> double override;
    auto valueAtCell(const mesh::Cell& cell) const -> double override;
    auto valueAtFace(std::size_t face_id) const -> double override;
    auto valueAtFace(const mesh::Face& face) const -> double override;

    auto parent() -> IVector*;
    void setParent(IVector* parent);

    auto inline operator[](std::size_t i) const -> double { return (*_data)[i]; }
    auto inline operator[](std::size_t i) -> double& { return (*_data)[i]; }

  protected:
    auto valueAtInteriorFace(const mesh::Face& face) const -> double;
    auto valueAtBoundaryFace(const mesh::Face& face) const -> double;

  private:
    void addDefaultHandlers();

    std::shared_ptr<VectorXd> _data = nullptr;
    std::shared_ptr<VectorXd> _face_data = nullptr;
    IVector* _parent = nullptr;
    std::optional<Coord> _coord = std::nullopt;
};

template <typename T>
concept ScalarBased = std::derived_from<T, IScalar>;

} // namespace prism::field