#pragma once

#include "ifield.h"

namespace prism::field {
class Tensor : public IField<Matrix3d> {
  public:
    Tensor(std::string name, const mesh::PMesh& mesh, double value);
    Tensor(std::string name, const mesh::PMesh& mesh, const Matrix3d& data);
    Tensor(std::string name, const mesh::PMesh& mesh, std::vector<Matrix3d> data);

    auto valueAtCell(std::size_t cell_id) const -> Matrix3d override;
    auto valueAtCell(const mesh::Cell& cell) const -> Matrix3d override;

    auto valueAtFace(std::size_t face_id) const -> Matrix3d override;
    auto valueAtFace(const mesh::Face& face) const -> Matrix3d override;

    auto operator[](std::size_t i) -> Matrix3d&;
    auto operator[](std::size_t i) const -> const Matrix3d&;

  private:
    std::vector<Matrix3d> _data;
};
} // namespace prism::field