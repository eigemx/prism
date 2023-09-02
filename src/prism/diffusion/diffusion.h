#pragma once

#include <memory>

#include "../field.h"
#include "../fvscheme.h"
#include "../gradient/gradient.h"
#include "../mesh/pmesh.h"
#include "../types.h"

namespace prism::diffusion {

enum class NonOrthoCorrection {
    None,
    OverRelaxed,
};

template <NonOrthoCorrection Corrector = NonOrthoCorrection::None>
class Diffusion : public FVScheme {
  public:
    // Explicit gradient scheme is provided by the user
    Diffusion(double kappa,
              ScalarField& phi,
              std::shared_ptr<gradient::GradientSchemeBase> gradient_scheme)
        : _kappa(kappa),
          _phi(phi),
          _mesh(phi.mesh()),
          _gradient_scheme(std::move(gradient_scheme)),
          FVScheme(phi.mesh().n_cells()) {}

    // TODO: make each FVScheme own its own gradient scheme, and remove the
    // following clutter.
    // Each FVScheme should own its own gradient scheme, and should not share
    // it with other FVSchemes. This is because each FVScheme may run on a
    // different thread or process, and the gradient scheme is not thread-safe.
    // and to also avoid the performance overhead of locking the gradient scheme.
    Diffusion(Diffusion&& other) noexcept = default;
    auto operator=(const Diffusion& other) -> Diffusion&;
    auto operator=(Diffusion&& other) noexcept -> Diffusion& = default;
    ~Diffusion() = default;

    void apply() override;
    auto requires_correction() const -> bool override;
    auto field() -> ScalarField& override { return _phi; }

  private:
    void apply_interior(const mesh::Face& face) override;
    void apply_boundary(const mesh::Cell& cell, const mesh::Face& face) override;
    void apply_boundary_fixed(const mesh::Cell& cell, const mesh::Face& face);
    void correct_non_orhto_boundary_fixed(const mesh::Cell& cell,
                                          const mesh::Face& face,
                                          const Vector3d& T_f);

    void apply_boundary_gradient(const mesh::Cell& cell, const mesh::Face& face);

    double _kappa;
    ScalarField& _phi;
    const mesh::PMesh& _mesh;
    std::shared_ptr<gradient::GradientSchemeBase> _gradient_scheme;
};

template <NonOrthoCorrection Corrector>
auto Diffusion<Corrector>::operator=(const Diffusion& other) -> Diffusion& {
    if (this != &other) {
        _kappa = other._kappa;
        _phi = other._phi;
        _mesh = other._mesh;
        _gradient_scheme = std::make_shared<gradient::GreenGauss>(other._phi);
    }
    return *this;
}

template <NonOrthoCorrection Corrector>
void inline Diffusion<Corrector>::apply() {
    /** @brief Applies discretized diffusion equation to the mesh.
     * The discretized equation is applied per face basis, using apply_interior() and 
     * apply_boundary() functions.
     * 
     * The function will check at first if the scheme has completed the first iteration,
     * and if the scheme does not require explicit correction, it will not re-calculate
     * the scheme coefficients and will not zero out the scheme matrix and RHS vector.
     */
    for (const auto& bface : _mesh.boundary_faces()) {
        apply_boundary(_mesh.cell(bface.owner()), bface);
    }

    for (const auto& iface : _mesh.interior_faces()) {
        apply_interior(iface);
    }

    // we've inserted all the triplets, now we can collect them into the matrix
    collect();
}

template <NonOrthoCorrection Corrector>
void Diffusion<Corrector>::apply_boundary(const mesh::Cell& cell, const mesh::Face& face) {
    /**
     * @brief Applies boundary discretized diffusion equation to the cell,
     * when the current face is a boundary face. The function iteself does not
     * apply the discretized equation, rather it calls the appropriate function
     * based on the boundary type.
     *
     * @param cell The cell that owns the boundary face.
     * @param face The boundary face.
     */
    const auto& boundary_patch = _mesh.face_boundary_patch(face);
    const auto& boundary_condition = boundary_patch.get_bc(_phi.name());

    switch (boundary_condition.bc_type()) {
        // empty boundary patch, do nothing
        case mesh::BoundaryConditionType::Empty: {
            return;
        }

        // fixed boundary patch, or Dirichlet boundary condition
        case mesh::BoundaryConditionType::Fixed:
        case mesh::BoundaryConditionType::Inlet: {
            apply_boundary_fixed(cell, face);
            return;
        }

        // Symmetry boundary patch, or zero gradient boundary condition.
        // This is a special case of the general Neumann boundary condition,
        // where the gradient of the field is zero at the boundary (flux is zero),
        // and will not result in any contribution to the right hand side of the equation,
        // or the matrix coefficients. and no need for non-orthogonal correction.
        // check equation 8.41 - Chapter 8 (Moukallad et al., 2015) and the following paragraph,
        // and paragraph 8.6.8.2 - Chapter 8 in same reference.
        case mesh::BoundaryConditionType::Symmetry:
        case mesh::BoundaryConditionType::Outlet: {
            return;
        }

        // general Von Neumann boundary condition, or fixed gradient boundary condition.
        case mesh::BoundaryConditionType::FixedGradient: {
            apply_boundary_gradient(cell, face);
            return;
        }

        default:
            throw std::runtime_error(
                fmt::format("diffusion::Diffusion::apply_boundary(): "
                            "Non-implemented boundary condition type for boundary patch: '{}'",
                            boundary_patch.name()));
    }
}

template <NonOrthoCorrection Corrector>
void Diffusion<Corrector>::apply_boundary_gradient(const mesh::Cell& cell,
                                                   const mesh::Face& face) {
    /** @brief Applies boundary discretized diffusion equation to the cell,
     * when the current face is a boundary face, and the boundary condition
     * is a general Von Neumann boundary condition, or fixed gradient boundary condition.
     *
     * @param cell The cell which owns the boundary face.
     * @param face The boundary face.
     */
    // get the fixed gradient (flux) value associated with the face
    const auto& boundary_patch = _mesh.face_boundary_patch(face);
    auto flux_wall = boundary_patch.get_scalar_bc(_phi.name());

    auto cell_id = cell.id();

    // check Moukallad et al 2015 Chapter 8 equation 8.39, 8.41 and the following paragraph,
    // and paragraph 8.6.8.2
    rhs(cell_id) += -flux_wall * face.area();
}


} // namespace prism::diffusion