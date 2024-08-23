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


} // namespace prism::field