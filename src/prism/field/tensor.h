#pragma once

#include "ifield.h"
#include "units.h"

namespace prism::field {
class Tensor : public ITensor, public units::Measurable {
  public:
    Tensor(std::string name, const SharedPtr<mesh::PMesh>& mesh, double value);

    Tensor(std::string name, const SharedPtr<mesh::PMesh>& mesh, const Matrix3d& data);

    Tensor(std::string name, const SharedPtr<mesh::PMesh>& mesh, std::vector<Matrix3d> data);

    auto valueAtCell(std::size_t cell_id) const -> Matrix3d override;
    auto valueAtCell(const mesh::Cell& cell) const -> Matrix3d override;

    auto valueAtCellRef(std::size_t cell_id) const -> const Matrix3d&;
    auto valueAtCellRef(const mesh::Cell& cell) const -> const Matrix3d&;
    auto cellValues() const -> const std::vector<Matrix3d>&;
    auto cellValues() -> std::vector<Matrix3d>&;

    auto valueAtFace(std::size_t face_id) const -> Matrix3d override;
    auto valueAtFace(const mesh::Face& face) const -> Matrix3d override;

    auto hasFaceValues() const -> bool override;
    void setFaceValues(std::vector<Matrix3d> values);
    void clearFaceValues();

    auto operator[](std::size_t i) -> Matrix3d&;
    auto operator[](std::size_t i) const -> const Matrix3d&;

    using ValueType = Matrix3d;

  private:
    std::vector<Matrix3d> _data;
    std::vector<Matrix3d> _face_data;
};
} // namespace prism::field
