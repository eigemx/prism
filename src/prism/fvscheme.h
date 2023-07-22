#pragma once

#include <optional>
#include <utility>

#include "field.h"
#include "mesh/pmesh.h"
#include "print.h"
#include "types.h"

namespace prism {

class FVScheme {
  public:
    FVScheme(std::size_t n_cells) : A(n_cells, n_cells), b(n_cells) { b.setZero(); }

    // apply the discretization scheme to `face` shared by `cell`
    virtual inline void apply(const mesh::Cell& cell, const mesh::Face& face) {
        if (face.is_boundary()) {
            apply_boundary(cell, face);
        } else {
            apply_interior(cell, face);
        }
    }

    // finalize() is called after all the cells are processed
    // before moving on to the next iteration
    // the default implementation does nothing
    virtual void finalize() {};

    inline auto coeff_matrix() const -> const SparseMatrix& { return A; }
    inline auto coeff_matrix() -> SparseMatrix& { return A; }

    inline auto rhs_vector() const -> const VectorXd& { return b; }
    inline auto rhs_vector() -> VectorXd& { return b; }

    // returns true if the scheme requires correction.
    // The default implementation returns ture.
    // Override this method if the scheme does not require correction.
    // If the scheme does not require correction, the solver will not zero out
    // the coefficient matrix and the RHS vector before the next iteration.
    virtual auto requires_correction() const -> bool { return true; }

    virtual auto field() -> ScalarField& = 0;

  private:
    virtual void apply_interior(const mesh::Cell& cell, const mesh::Face& face) = 0;
    virtual void apply_boundary(const mesh::Cell& cell, const mesh::Face& face) = 0;

    SparseMatrix A;
    VectorXd b;
};

} // namespace prism