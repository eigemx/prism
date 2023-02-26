#pragma once

#include "../fvscheme.h"
#include "../mesh/pmesh.h"
#include "../types.h"

namespace prism {

class DiffusionSchemeBase {};

class LinearDiffusionScheme : public FVScheme, public DiffusionSchemeBase {
  public:
    LinearDiffusionScheme(double diffusion_coefficient,
                          const mesh::PMesh& mesh,
                          SparseMatrix& coeffs,
                          VectorXd& b)
        : _diffusion_coefficient(diffusion_coefficient), _mesh(mesh), _coeffs(coeffs), _b(b) {}


    auto run_once() const -> bool override { return true; }

  private:
    void apply_interior(const mesh::Cell& cell, const mesh::Face& face) const override;
    void apply_boundary(const mesh::Cell& cell, const mesh::Face& face) const override;

    double _diffusion_coefficient;
    const mesh::PMesh& _mesh;
    SparseMatrix& _coeffs;
    VectorXd& _b;
};

} // namespace prism