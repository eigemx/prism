#pragma once

#include "boundary.h"
#include "ifield.h"

namespace prism::field {

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
} // namespace prism::field