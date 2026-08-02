#pragma once

#include "ifield.h"
#include "units.h"

namespace prism::field {
class Tensor : public ITensor, public units::Measurable {
  public:
    Tensor(std::string name, const SharedPtr<mesh::PMesh>& mesh, double value);

    Tensor(std::string name, const SharedPtr<mesh::PMesh>& mesh, const Tensor3d& data);

    Tensor(std::string name, const SharedPtr<mesh::PMesh>& mesh, std::vector<Tensor3d> data);

    auto valueAtCell(size_t cell_id) const -> Tensor3d override;
    auto valueAtCell(const mesh::Cell& cell) const -> Tensor3d override;

    auto valueAtCellRef(size_t cell_id) const -> const Tensor3d&;
    auto valueAtCellRef(const mesh::Cell& cell) const -> const Tensor3d&;
    auto cellValues() const -> const std::vector<Tensor3d>&;
    auto cellValues() -> std::vector<Tensor3d>&;

    auto valueAtFace(size_t face_id) const -> Tensor3d override;
    auto valueAtFace(const mesh::Face& face) const -> Tensor3d override;

    auto hasFaceValues() const -> bool override;
    void setFaceValues(std::vector<Tensor3d> values);
    void clearFaceValues();

    auto operator[](size_t i) -> Tensor3d&;
    auto operator[](size_t i) const -> const Tensor3d&;

    using ValueType = Tensor3d;

  private:
    std::vector<Tensor3d> _data;
    std::vector<Tensor3d> _face_data;
};
} // namespace prism::field
