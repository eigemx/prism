#pragma once

#include <memory>

#include "../fvscheme.h"
#include "../mesh/pmesh.h"
#include "../types.h"

namespace prism::diffusion {

class DiffusionSchemeBase {};

class Linear : public FVScheme, public DiffusionSchemeBase {
  public:
    Linear(double kappa) : _kappa(kappa) {}


  private:
    void apply_interior(const mesh::Cell& cell,
                        const mesh::Face& face,
                        const mesh::PMesh& mesh,
                        SparseMatrix& coeffs_matrix,
                        VectorXd& rhs_vec) const override;

    void apply_boundary(const mesh::Cell& cell,
                        const mesh::Face& face,
                        const mesh::PMesh& mesh,
                        SparseMatrix& coeffs_matrix,
                        VectorXd& rhs_vec) const override;

    void apply_wall_boundary(const mesh::Cell& cell,
                             const mesh::Face& face,
                             const mesh::PMesh& mesh,
                             SparseMatrix& coeffs_matrix,
                             VectorXd& rhs_vec) const;

    double _kappa;
};

} // namespace prism::diffusion