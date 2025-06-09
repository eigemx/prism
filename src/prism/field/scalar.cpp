#include "scalar.h"

#include "prism/log.h"

namespace prism::field {

UniformScalar::UniformScalar(std::string name, const mesh::PMesh& mesh, double value)
    : IScalar(std::move(name), mesh), _value(value) {
    log::debug("Creating uniform scalar field: '{}' with double value = {}", this->name(), value);
}

auto UniformScalar::valueAtCell(std::size_t cell_id) const -> double { // NOLINT
    return _value;
}

auto UniformScalar::valueAtCell(const mesh::Cell& cell) const -> double { // NOLINT
    return _value;
}

auto UniformScalar::valueAtFace(std::size_t face_id) const -> double { // NOLINT
    return _value;
}

auto UniformScalar::valueAtFace(const mesh::Face& face) const -> double { // NOLINT
    return _value;
}

auto UniformScalar::gradAtFace(const mesh::Face& face) const -> Vector3d { // NOLINT
    return {0.0, 0.0, 0.0};
}

auto UniformScalar::gradAtCell(const mesh::Cell& cell) const -> Vector3d { // NOLINT
    return {0.0, 0.0, 0.0};
}

auto UniformScalar::gradAtCellStored(const mesh::Cell& cell) const -> Vector3d { // NOLINT
    return {0.0, 0.0, 0.0};
}

void ScalarBHManagerSetter::set(IScalarBHManager& manager) {
    log::debug(
        "prism::field::ScalarBHManagerSetter::set(): adding default boundary handlers for a "
        "scalar field instance");
    manager.addHandler<field::boundary::Fixed<Scalar>>();
    manager.addHandler<field::boundary::Symmetry<Scalar>>();
    manager.addHandler<field::boundary::Outlet<Scalar>>();
    manager.addHandler<field::boundary::FixedGradient<Scalar>>();
    manager.addHandler<field::boundary::ZeroGradient<Scalar>>();
}

} // namespace prism::field