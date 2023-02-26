#pragma once

#include "../fvscheme.h"
#include "../mesh/pmesh.h"
#include "../types.h"

namespace prism {
class GradientSchemeBase {};

class LeastSquaresGradientScheme : public FVScheme, public GradientSchemeBase {
  public:
    LeastSquaresGradientScheme(const mesh::PMesh& mesh, SparseMatrix& coeffs, VectorXd& b)
        : _mesh(mesh), _coeffs(coeffs), _b(b) {}

    auto run_once() const -> bool override { return false; }

  private:
    void apply_interior(const mesh::Cell& cell, const mesh::Face& face) const override;
    void apply_boundary(const mesh::Cell& cell, const mesh::Face& face) const override;

    const mesh::PMesh& _mesh;
    SparseMatrix& _coeffs;
    VectorXd& _b;

    bool _first_run_completed {false};
};
} // namespace prism