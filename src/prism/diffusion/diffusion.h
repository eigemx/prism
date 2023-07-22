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

// TODO: Implement a non-corrected version
template <NonOrthoCorrection Corrector = NonOrthoCorrection::None>
class Diffusion : public FVScheme {
  public:
    Diffusion() = delete;

    // Explicit gradient scheme is provided by the user
    Diffusion(double kappa,
              ScalarField& phi,
              std::shared_ptr<gradient::GradientSchemeBase> gradient_scheme)
        : _kappa(kappa),
          _phi(phi),
          _mesh(phi.mesh()),
          _gradient_scheme(std::move(gradient_scheme)),
          FVScheme(phi.mesh().n_cells()) {}

    Diffusion(const Diffusion& other)
        : _kappa(other._kappa),
          _phi(other._phi),
          _mesh(other._mesh),
          _gradient_scheme(other._gradient_scheme->clone()),
          FVScheme(other._mesh.n_cells()) {}

    Diffusion(Diffusion&& other) noexcept = default;

    auto operator=(const Diffusion& other) -> Diffusion& {
        if (this != &other) {
            _kappa = other._kappa;
            _phi = other._phi;
            _mesh = other._mesh;
            _gradient_scheme = std::make_shared<gradient::GreenGauss>(other._phi);
        }
        return *this;
    }

    auto operator=(Diffusion&& other) noexcept -> Diffusion& = default;

    ~Diffusion() = default;

    auto requires_correction() const -> bool override;

    auto field() -> ScalarField& override { return _phi; }

  private:
    void apply_interior(const mesh::Cell& cell, const mesh::Face& face) override;
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
void Diffusion<Corrector>::apply_boundary(const mesh::Cell& cell, const mesh::Face& face) {
    /**
     * @brief Applies boundary discretized diffusion equation to the cell,
     * when the current face is a boundary face. The function iteself does not
     * apply the discretized equation, rather it calls the appropriate function
     * based on the boundary type.
     *
     * @param cell The cell to which the boundary conditions are applied.
     * @param face The boundary face.
     */
    const auto& boundary_patch = _mesh.face_boundary_patch(face);
    const auto& boundary_condition = boundary_patch.get_bc(_phi.name());

    switch (boundary_condition.patch_type()) {
        // empty boundary patch, do nothing
        case mesh::BoundaryPatchType::Empty: {
            return;
        }

        // fixed boundary patch, or Dirichlet boundary condition
        case mesh::BoundaryPatchType::Fixed:
        case mesh::BoundaryPatchType::Inlet: {
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
        case mesh::BoundaryPatchType::Symmetry:
        case mesh::BoundaryPatchType::Outlet: {
            return;
        }

        // general Von Neumann boundary condition, or fixed gradient boundary condition.
        case mesh::BoundaryPatchType::FixedGradient: {
            apply_boundary_gradient(cell, face);
            return;
        }

        default:
            throw std::runtime_error(
                fmt::format("diffusion::Linear::apply_boundary(): "
                            "Non-implemented boundary type for boundary patch: '{}'",
                            boundary_patch.name()));
    }
}

template <NonOrthoCorrection Corrector>
void Diffusion<Corrector>::apply_boundary_gradient(const mesh::Cell& cell,
                                                   const mesh::Face& face) {
    // get the fixed gradient (flux) value associated with the face
    const auto& boundary_patch = _mesh.face_boundary_patch(face);
    auto flux_wall = boundary_patch.get_scalar_bc(_phi.name());

    auto cell_id = cell.id();

    // check Moukallad et al 2015 Chapter 8 equation 8.39, 8.41 and the following paragraph,
    // and paragraph 8.6.8.2
    rhs_vector()[cell_id] += -flux_wall * face.area();
}


} // namespace prism::diffusion