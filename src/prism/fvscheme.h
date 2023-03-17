#pragma once

#include <optional>
#include <utility>

#include "mesh/pmesh.h"
#include "types.h"

namespace prism {

class FVScheme {
  public:
    virtual inline void apply(const mesh::Cell& cell, const mesh::Face& face) {
        if (face.is_boundary()) {
            apply_boundary(cell, face);
        } else {
            apply_interior(cell, face);
        }
    }

    inline void init_linear_system(const mesh::PMesh& mesh) {
        auto n_cells = mesh.cells().size();
        A = SparseMatrix(n_cells, n_cells);
        b = VectorXd::Zero(n_cells);
    }

    virtual void finalize() = 0;

    inline auto coeff_matrix() const -> const SparseMatrix& { return A; }
    inline auto coeff_matrix() -> SparseMatrix& { return A; }

    inline auto rhs_vector() const -> const VectorXd& { return b; }
    inline auto rhs_vector() -> VectorXd& { return b; }

    virtual auto mesh() const -> const mesh::PMesh& = 0;


  private:
    virtual void apply_interior(const mesh::Cell& cell, const mesh::Face& face) = 0;
    virtual void apply_boundary(const mesh::Cell& cell, const mesh::Face& face) = 0;

    SparseMatrix A;
    VectorXd b;
};

} // namespace prism