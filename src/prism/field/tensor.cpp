#include "tensor.h"

#include "prism/log.h"
#include "prism/mesh/utilities.h"

namespace prism::field {
Tensor::Tensor(std::string name, const SharedPtr<mesh::PMesh>& mesh, double value)
    : Tensor(std::move(name), mesh, Matrix3d::Identity() * value) {}

Tensor::Tensor(std::string name, const SharedPtr<mesh::PMesh>& mesh, const Matrix3d& data)
    : ITensor(std::move(name), mesh) {
    log::debug("Creating a uniform tensor field: '{}' given a Matrix3d object", this->name());

    _data.assign(mesh->cellCount(), data);
}

Tensor::Tensor(std::string name, const SharedPtr<mesh::PMesh>& mesh, std::vector<Matrix3d> data)
    : ITensor(std::move(name), mesh), _data(std::move(data)) {
    log::debug("Creating a tensor field: '{}' given a vector of Matrix3d objects", this->name());

    if (_data.size() != mesh->cellCount()) {
        throw std::runtime_error(
            fmt::format("field::Tensor() cannot create a tensor field '{}' given a vector of "
                        "Matrix3d that has a different size than mesh cells count.",
                        this->name()));
    }
}

auto Tensor::valueAtCell(std::size_t cell_id) const -> Matrix3d {
    return _data[cell_id];
}

auto Tensor::valueAtCell(const mesh::Cell& cell) const -> Matrix3d {
    return valueAtCell(cell.id());
}

auto Tensor::valueAtCellRef(std::size_t cell_id) const -> const Matrix3d& {
    return _data[cell_id];
}

auto Tensor::valueAtCellRef(const mesh::Cell& cell) const -> const Matrix3d& {
    return _data[cell.id()];
}

auto Tensor::cellValues() const -> const std::vector<Matrix3d>& {
    return _data;
}

auto Tensor::cellValues() -> std::vector<Matrix3d>& {
    return _data;
}

auto Tensor::valueAtFace(std::size_t face_id) const -> Matrix3d {
    const mesh::Face& face = this->mesh()->face(face_id);
    return valueAtFace(face);
}

auto Tensor::valueAtFace(const mesh::Face& face) const -> Matrix3d {
    if (hasFaceValues()) {
        return _face_data[face.id()];
    }

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

    return (gc * _data[owner.id()]) + ((1 - gc) * _data[neighbor.id()]);
}

auto Tensor::hasFaceValues() const -> bool {
    return !_face_data.empty();
}

void Tensor::setFaceValues(std::vector<Matrix3d> values) {
    if (values.size() != mesh()->faceCount()) {
        throw std::runtime_error(fmt::format(
            "prism::field::Tensor::setFaceValues(): cannot set face values for "
            "tensor field {}, to a face data vector having a different size that field's "
            "faces count.",
            name()));
    }

    if (hasFaceValues()) {
        log::debug(
            "Tensor::setFaceValues() was called for field '{}', which already has face values set.",
            name());
    }

    log::debug("Setting face values for field '{}'", name());
    _face_data = std::move(values);
}

void Tensor::clearFaceValues() {
    if (_face_data.empty()) {
        log::warn(
            "Tensor::clearFaceValues() was called for field '{}', but the face data is not "
            "initialized.",
            name());
        return;
    }
    _face_data.clear();
    _face_data.shrink_to_fit();
}

auto Tensor::operator[](std::size_t i) -> Matrix3d& {
    return _data[i];
}

auto Tensor::operator[](std::size_t i) const -> const Matrix3d& {
    return _data[i];
}

} // namespace prism::field
