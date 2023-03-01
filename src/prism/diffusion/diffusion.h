#pragma once

#include "../fvscheme.h"
#include "../mesh/pmesh.h"
#include "../types.h"

namespace prism::diffusion {

class DiffusionSchemeBase {};

class Linear : public FVScheme, public DiffusionSchemeBase {
  public:
    Linear(double kappa, const mesh::PMesh& mesh, SparseMatrix& coeffs, VectorXd& b)
        : _kappa(kappa), _mesh(mesh), _coeffs(coeffs), _b(b) {}


    auto run_once() const -> bool override { return true; }

  private:
    void apply_interior(const mesh::Cell& cell, const mesh::Face& face) const override;
    void apply_boundary(const mesh::Cell& cell, const mesh::Face& face) const override;
    void apply_wall_boundary(const mesh::Cell& cell,
                             const mesh::Face& face,
                             double phi_wall) const;

    double _kappa;
    const mesh::PMesh& _mesh;
    SparseMatrix& _coeffs;
    VectorXd& _b;
};

} // namespace prism::diffusion