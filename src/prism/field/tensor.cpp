#include "tensor.h"

#include "prism/log.h"
#include "prism/mesh/utilities.h"

namespace prism::field {
Tensor::Tensor(std::string name, const SharedPtr<mesh::PMesh>& mesh, double value)
    : IField(std::move(name), mesh) {
    log::debug("Creating tensor field: '{}' with double value = {}", this->name(), value);
    const std::size_t n_cells = this->mesh()->cellCount();
    _data.reserve(n_cells);
    for (std::size_t i = 0; i < n_cells; ++i) {
        _data.emplace_back(Matrix3d::Identity() * value);
    }
}

Tensor::Tensor(std::string name, const SharedPtr<mesh::PMesh>& mesh, const Matrix3d& data)
    : IField(std::move(name), mesh) {
    log::debug("Creating a uniform tensor field: '{}' given a Matrix3d object", this->name());

    const std::size_t n_cells = this->mesh()->cellCount();
    _data.reserve(n_cells);
    for (std::size_t i = 0; i < n_cells; ++i) {
        _data.push_back(data);
    }
}

Tensor::Tensor(std::string name, const SharedPtr<mesh::PMesh>& mesh, std::vector<Matrix3d> data)
    : IField(std::move(name), mesh), _data(std::move(data)) {
    log::debug("Creating a tensor field: '{}' given a vector of Matrix3d objects", this->name());

    if (_data.size() != mesh->cellCount()) {
        throw std::runtime_error(
            fmt::format("field::Tensor() cannot create a tensor field '{}' given a vector of "
                        "Matrix3d that has a different size than mesh cells count.",
                        this->name()));
    }
}

auto Tensor::valueAtCell(std::size_t cell_id) const -> Matrix3d {
    assert(cell_id < mesh()->cellCount());
    assert(cell_id < _data.size());
    return _data[cell_id];
}

auto Tensor::valueAtCell(const mesh::Cell& cell) const -> Matrix3d {
    return valueAtCell(cell.id());
}

auto Tensor::valueAtFace(std::size_t face_id) const -> Matrix3d {
    const mesh::Face& face = this->mesh()->face(face_id);
    return valueAtFace(face);
}

auto Tensor::valueAtFace(const mesh::Face& face) const -> Matrix3d {
    const auto& mesh = this->mesh();
    const mesh::Cell& owner = mesh->cell(face.owner());

    if (face.isBoundary()) {
        log::warn(
            "field::Tensor::valueAtFace() was called on a boundary face (face id = {}). "
            "Returning value of the tensor field at owner cell.",
            face.id());
        return _data[owner.id()];
    }
    const mesh::Cell& neighbor = mesh->otherSharingCell(owner, face);
    const double gc = mesh::geometricWeight(owner, neighbor, face);

    /// TODO: we should replace this with a more general interpolation function
    return (gc * _data[owner.id()]) + ((1 - gc) * _data[neighbor.id()]);
}

auto Tensor::operator[](std::size_t i) -> Matrix3d& {
    return _data[i];
}

auto Tensor::operator[](std::size_t i) const -> const Matrix3d& {
    return _data[i];
}

} // namespace prism::field
