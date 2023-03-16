#pragma once

#include <optional>
#include <utility>

#include "mesh/pmesh.h"
#include "types.h"

namespace prism {

class FVScheme {
  public:
    virtual void apply(const mesh::Cell& cell,
                       const mesh::Face& face,
                       const mesh::PMesh& mesh,
                       SparseMatrix& coeffs_matrix,
                       VectorXd& rhs_vec) const {
        if (face.has_neighbor()) {
            apply_interior(cell, face, mesh, coeffs_matrix, rhs_vec);
        } else {
            apply_boundary(cell, face, mesh, coeffs_matrix, rhs_vec);
        }
    }

  private:
    virtual void apply_interior(const mesh::Cell& cell,
                                const mesh::Face& face,
                                const mesh::PMesh& mesh,
                                SparseMatrix& coeffs_matrix,
                                VectorXd& rhs_vec) const = 0;

    virtual void apply_boundary(const mesh::Cell& cell,
                                const mesh::Face& face,
                                const mesh::PMesh& mesh,
                                SparseMatrix& coeffs_matrix,
                                VectorXd& rhs_vec) const = 0;
};

} // namespace prism