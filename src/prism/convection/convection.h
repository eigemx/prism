#pragma once

#include <cmath>
#include <memory>

#include "../field.h"
#include "../fvscheme.h"
#include "../gradient/gradient.h"
#include "../mesh/pmesh.h"

namespace prism::convection {

class ConvectionSchemeBase {};

// coefficients for the discretized convection equation for a face
struct CoeffsTriplet {
    double a_C {}; // cell
    double a_N {}; // neighbor
    double b {};   // source
};

// Base class for convection face interpolation schemes (upwind, central difference, ...)
class FaceInterpolationSchemeBase {
  public:
    virtual auto interpolate(double m_dot,
                             const mesh::Cell& cell,
                             const mesh::Cell& neighbor,
                             const mesh::Face& face,
                             const std::shared_ptr<gradient::GradientSchemeBase>& grad_scheme)
        -> CoeffsTriplet = 0;
};

// Central difference scheme
class CentralDifference : public FaceInterpolationSchemeBase {
  public:
    auto interpolate(double m_dot,
                     const mesh::Cell& cell,
                     const mesh::Cell& neighbor,
                     const mesh::Face& face,
                     const std::shared_ptr<gradient::GradientSchemeBase>& grad_scheme)
        -> CoeffsTriplet override;
};

// Upwind scheme
class Upwind : public FaceInterpolationSchemeBase {
  public:
    auto interpolate(double m_dot,
                     const mesh::Cell& cell,
                     const mesh::Cell& neighbor,
                     const mesh::Face& face,
                     const std::shared_ptr<gradient::GradientSchemeBase>& grad_scheme)
        -> CoeffsTriplet override;
};

// Second order upwind scheme
class SecondOrderUpwind : public FaceInterpolationSchemeBase {
  public:
    auto interpolate(double m_dot,
                     const mesh::Cell& cell,
                     const mesh::Cell& neighbor,
                     const mesh::Face& face,
                     const std::shared_ptr<gradient::GradientSchemeBase>& grad_scheme)
        -> CoeffsTriplet override;
};

// QUICK scheme
class QUICK : public FaceInterpolationSchemeBase {
  public:
    auto interpolate(double m_dot,
                     const mesh::Cell& cell,
                     const mesh::Cell& neighbor,
                     const mesh::Face& face,
                     const std::shared_ptr<gradient::GradientSchemeBase>& grad_scheme)
        -> CoeffsTriplet override;
};

// Finite volume scheme for the discretization of the convection term
template <typename FaceInterpolationScheme = Upwind>
class Convection : public FVScheme, public ConvectionSchemeBase {
  public:
    // Gradient scheme is not given, using Green-Gauss as a default explicit gradient scheme
    Convection(double rho, VectorField& U, ScalarField& phi)
        : _rho(rho),
          _U(U),
          _phi(phi),
          _mesh(phi.mesh()),
          _gradient_scheme(std::make_shared<gradient::GreenGauss>(phi)) {}

    // Gradient scheme is given
    Convection(double rho,
               VectorField& U,
               ScalarField& phi,
               std::shared_ptr<gradient::GradientSchemeBase> grad_scheme)
        : _rho(rho),
          _U(U),
          _phi(phi),
          _mesh(phi.mesh()),
          _gradient_scheme(std::move(grad_scheme)) {}

    void inline finalize() override { _main_coeffs_calculated = true; }

  private:
    void apply_interior(const mesh::Cell& cellell, const mesh::Face& face) override;
    void apply_boundary(const mesh::Cell& cellell, const mesh::Face& face) override;
    void apply_boundary_fixed(const mesh::Cell& cellell, const mesh::Face& face);
    void apply_boundary_outlet(const mesh::Cell& cellell, const mesh::Face& face);
    auto boundary_face_velocity(const mesh::Face& face) const -> Vector3d;

    double _rho;
    VectorField& _U;
    ScalarField& _phi;
    const mesh::PMesh& _mesh;
    FaceInterpolationScheme _face_interpolation_scheme;
    std::shared_ptr<gradient::GradientSchemeBase> _gradient_scheme;
    bool _main_coeffs_calculated {false};
};


auto inline face_mass_flow_rate(double rho, const Vector3d& U, const Vector3d& S) -> double {
    return rho * U.dot(S);
}

template <typename F>
void Convection<F>::apply_interior(const mesh::Cell& cell, const mesh::Face& face) {
    auto cell_id = cell.id();

    // get adjacent cell sharing `face` with `cell`
    auto is_owned = face.is_owned_by(cell_id);
    auto adjacent_cell_id = is_owned ? face.neighbor().value() : face.owner();

    const auto& adj_cell = _mesh.cell(adjacent_cell_id);

    // face area vector, always pointing outwards
    auto S_f = face.area_vector();

    if (!is_owned) {
        // face is not owned by `cell`, flip the area vector
        // so that it points outwards from `cell`
        S_f *= -1;
    }

    // interpolated velocity vector at face centroid
    const auto& U_f = 0.5 * (_U[cell_id] + _U[adjacent_cell_id]);

    auto m_dot_f = face_mass_flow_rate(_rho, U_f, S_f);

    auto [a_C, a_N, b] =
        _face_interpolation_scheme.interpolate(m_dot_f, cell, adj_cell, face, _gradient_scheme);


    if (!_main_coeffs_calculated) {
        coeff_matrix().coeffRef(cell_id, cell_id) += a_C;
        coeff_matrix().coeffRef(cell_id, adjacent_cell_id) += a_N;
    }

    rhs_vector().coeffRef(cell_id) += b;
}

template <typename F>
void Convection<F>::apply_boundary(const mesh::Cell& cell, const mesh::Face& face) {
    const auto& boundary_patch = _mesh.face_boundary_patch(face);
    const auto& boundary_condition = boundary_patch.get_bc(_phi.name());

    switch (boundary_condition.patch_type()) {
        case mesh::BoundaryPatchType::Empty: {
            return;
        }

        case mesh::BoundaryPatchType::Fixed:
        case mesh::BoundaryPatchType::Inlet: {
            apply_boundary_fixed(cell, face);
            return;
        }

        case mesh::BoundaryPatchType::Outlet:
        case mesh::BoundaryPatchType::Symmetry: {
            apply_boundary_outlet(cell, face);
            return;
        }

        default:
            throw std::runtime_error(
                format("convection::Convection::apply_boundary(): "
                       "Non-implemented boundary type for boundary patch: '{}' for field '{}'",
                       boundary_patch.name(),
                       _U.name()));
    }
}

template <typename F>
void Convection<F>::apply_boundary_fixed(const mesh::Cell& cell, const mesh::Face& face) {
    const auto& boundary_patch = _mesh.face_boundary_patch(face);
    auto phi_wall = boundary_patch.get_scalar_bc(_phi.name());

    const auto& S_f = face.area_vector();
    const auto& U_f = boundary_face_velocity(face);
    auto m_dot_f = face_mass_flow_rate(_rho, U_f, S_f);

    if (!_main_coeffs_calculated) {
        // in case owner cell is an upstream cell
        coeff_matrix().coeffRef(cell.id(), cell.id()) += std::max(m_dot_f, 0.0);
    }

    rhs_vector().coeffRef(cell.id()) += std::max(-m_dot_f * phi_wall, 0.0);
}

template <typename F>
void Convection<F>::apply_boundary_outlet(const mesh::Cell& cell, const mesh::Face& face) {
    const auto& boundary_patch = _mesh.face_boundary_patch(face);

    // The outlet boundary condition is basically a zero gradient condition,
    // so the value of the field at the outlet face centroid is the same as
    // the value of the field at the owner cell centroid.
    auto phi_outlet = _phi[face.owner()];
    auto cell_id = cell.id();

    // face area vector
    const auto& S_f = face.area_vector();

    // use owner cell velocity as the velocity at the outlet face centroid
    const auto& U_f = boundary_face_velocity(face);

    auto m_dot_f = face_mass_flow_rate(_rho, U_f, S_f);

    if (!_main_coeffs_calculated) {
        // This is an outlet, so owner cell is an upstream cell
        if (m_dot_f < 0.0) {
            warn(
                format("convection::Convection::apply_boundary_outlet(): "
                       "Reverse flow detected at outlet boundary patch '{}'. "
                       "This may cause the solution to diverge.",
                       boundary_patch.name()));
        }
        // TODO: this assumes an upwind based scheme, this is wrong for central schemes
        // and should be generalized to work for all schemes.
        coeff_matrix().coeffRef(cell_id, cell_id) += m_dot_f;
    }
}

template <typename F>
auto Convection<F>::boundary_face_velocity(const mesh::Face& face) const -> Vector3d {
    const auto& boundary_patch = _mesh.face_boundary_patch(face);
    const auto& boundary_condition = boundary_patch.get_bc(_U.name());

    switch (boundary_condition.patch_type()) {
        case mesh::BoundaryPatchType::Fixed:
        case mesh::BoundaryPatchType::Inlet: {
            return boundary_patch.get_vector_bc(_U.name());
        }

        case mesh::BoundaryPatchType::Outlet:
        case mesh::BoundaryPatchType::Symmetry: {
            return _U[face.owner()];
        }

        default:
            throw std::runtime_error(
                format("convection::Convection::boundary_face_velocity(): "
                       "Non-implemented boundary type for boundary patch: '{}' for field '{}'",
                       boundary_patch.name(),
                       _U.name()));
    }
}

} // namespace prism::convection
