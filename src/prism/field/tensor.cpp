#include "tensor.h"

#include <spdlog/spdlog.h>

#include "prism/mesh/utilities.h"

namespace prism::field {
Tensor::Tensor(std::string name, const mesh::PMesh& mesh, double value)
    : IField(std::move(name), mesh) {
    spdlog::debug("Creating tensor field: '{}' with double value = {}", this->name(), value);
    const std::size_t n_cells = this->mesh().nCells();
    _data.reserve(n_cells);
    for (std::size_t i = 0; i < n_cells; ++i) {
        _data.emplace_back(Matrix3d::Ones() * value);
    }
}

Tensor::Tensor(std::string name, const mesh::PMesh& mesh, const Matrix3d& data)
    : IField(std::move(name), mesh) {
    spdlog::debug("Creating a uniform tensor field: '{}' given a Matrix3d object", this->name());

    const std::size_t n_cells = this->mesh().nCells();
    _data.reserve(n_cells);
    for (std::size_t i = 0; i < n_cells; ++i) {
        _data.push_back(data);
    }
}

Tensor::Tensor(std::string name, const mesh::PMesh& mesh, std::vector<Matrix3d> data)
    : IField(std::move(name), mesh), _data(std::move(data)) {
    spdlog::debug("Creating a  tensor field: '{}' given a vector of Matrix3d objects",
                  this->name());

    if (_data.size() != mesh.nCells()) {
        throw std::runtime_error(
            fmt::format("field::Tensor() cannot create a tensor field '{}' given a vector of "
                        "Matrix3d that has a different size than mesh's cell count.",
                        this->name()));
    }
}

auto Tensor::valueAtCell(std::size_t cell_id) const -> Matrix3d {
    assert(cell_id < mesh().nCells());
    return _data[cell_id];
}

auto Tensor::valueAtCell(const mesh::Cell& cell) const -> Matrix3d {
    return _data[cell.id()];
}

auto Tensor::valueAtFace(std::size_t face_id) const -> Matrix3d {
    const mesh::Face& face = this->mesh().face(face_id);
    return valueAtFace(face);
}

auto Tensor::valueAtFace(const mesh::Face& face) const -> Matrix3d {
    const auto& mesh = this->mesh();
    const mesh::Cell& owner = mesh.cell(face.owner());

    if (face.isBoundary()) {
        spdlog::warn(
            "field::Tensor::valueAtFace() was called on a boundary face (face id = {}). "
            "Returning value of the tensor field at owner cell.",
            face.id());

        return _data[owner.id()];
    }
    const mesh::Cell& neighbor = mesh.cell(face.neighbor().value());
    const double gc = mesh::geometricWeight(owner, neighbor, face);

    return (gc * _data[owner.id()]) + ((1 - gc) * _data[neighbor.id()]);
}

auto Tensor::operator[](std::size_t i) -> Matrix3d& {
    return _data[i];
}

auto Tensor::operator[](std::size_t i) const -> const Matrix3d& {
    return _data[i];
}

} // namespace prism::field