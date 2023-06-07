#pragma once

#include "../field.h"
#include "../fvscheme.h"
#include "../gradient/gradient.h"
#include "../mesh/pmesh.h"

namespace prism::convection {

class ConvectionSchemeBase {};

struct InterpolationCoeffs {
    double a_C {};
    double a_D {};
    double b {};
};

class FaceInterpolationSchemeBase {
  public:
    virtual auto interpolate(double m_dot,
                             const mesh::Cell& cell,
                             const mesh::Cell& neighbor,
                             const mesh::Face& face) -> InterpolationCoeffs = 0;
};

class CentralDifference : public FaceInterpolationSchemeBase {
  public:
    auto interpolate(double m_dot,
                     const mesh::Cell& cell,
                     const mesh::Cell& neighbor,
                     const mesh::Face& face) -> InterpolationCoeffs override;
};

class Upwind : public FaceInterpolationSchemeBase {
  public:
    auto interpolate(double m_dot,
                     const mesh::Cell& cell,
                     const mesh::Cell& neighbor,
                     const mesh::Face& face) -> InterpolationCoeffs override;
};

template <typename FaceInterpolationScheme = CentralDifference,
          typename GradientScheme = gradient::GreenGauss>
class Convection : public FVScheme, public ConvectionSchemeBase {
  public:
    Convection(double rho, VectorField& U, ScalarField& phi)
        : _rho(rho), _U(U), _phi(phi), _mesh(phi.mesh()) {}

    void inline finalize() override {}

  private:
    void apply_interior(const mesh::Cell& cellell, const mesh::Face& face) override;
    void apply_boundary(const mesh::Cell& cellell, const mesh::Face& face) override;
    void apply_boundary_fixed(const mesh::Cell& cellell, const mesh::Face& face);
    void apply_boundary_outlet(const mesh::Cell& cellell, const mesh::Face& face);
    auto face_gradient(const mesh::Cell& c, const mesh::Cell& n, const mesh::Face& f) -> Vector3d;

    double _rho;
    VectorField& _U;
    ScalarField& _phi;
    mesh::PMesh _mesh;
    FaceInterpolationScheme _face_interpolation_scheme;
    GradientScheme _gradient_scheme;
};

auto inline face_mass_flow_rate(double rho, const Vector3d& U, const Vector3d& S) -> double {
    return rho * U.dot(S);
}

template <typename F, typename G>
auto Convection<F, G>::face_gradient(const mesh::Cell& c,
                                     const mesh::Cell& n,
                                     const mesh::Face& f) -> Vector3d {
    const auto& mesh = _phi.mesh();
    auto gc = mesh.cells_weighting_factor(c, n, f);

    const auto& grad_c = _grad
}

template <typename F, typename G>
void Convection<F, G>::apply_interior(const mesh::Cell& cell, const mesh::Face& face) {
    auto cell_id = cell.id();

    // get adjacent cell sharing `face` with `cell`
    auto adjacent_cell_id = face.is_owned_by(cell_id) ? face.neighbor().value() : face.owner();

    const auto& adj_cell = _mesh.cell(adjacent_cell_id);

    // face area vector
    const auto& S_f = face.area_vector();

    // interpolated velocity vector at face centroid
    const auto& U_f = 0.5 * (_U[cell_id] + _U[adjacent_cell_id]);

    auto m_dot_f = face_mass_flow_rate(_rho, U_f, S_f);

    auto [a_C, a_D, b] = _face_interpolation_scheme.interpolate(m_dot_f, cell, adj_cell, face);

    //coeff_matrix().coeffRef(cell_id, cell_id) += m_dot_f / 2;
    //coeff_matrix().coeffRef(cell_id, adjacent_cell_id) += m_dot_f / 2;
    coeff_matrix().coeffRef(cell_id, cell_id) += a_C;
    coeff_matrix().coeffRef(cell_id, adjacent_cell_id) += a_D;
    rhs_vector().coeffRef(cell_id) += b;
}

template <typename F, typename G>
void Convection<F, G>::apply_boundary(const mesh::Cell& cell, const mesh::Face& face) {
    const auto& boundary_patch = _mesh.face_boundary_patch(face);
    const auto& boundary_condition = boundary_patch.get_bc(_U.name());

    switch (boundary_condition.patch_type()) {
        case mesh::BoundaryPatchType::Empty: {
            return;
        }

        case mesh::BoundaryPatchType::Fixed:
        case mesh::BoundaryPatchType::Inlet: {
            apply_boundary_fixed(cell, face);
            return;
        }

        case mesh::BoundaryPatchType::Outlet: {
            apply_boundary_outlet(cell, face);
            return;
        }

        default:
            throw std::runtime_error(
                format("convection::Convection::apply_boundary(): "
                       "Non-implemented boundary type for boundary patch: '{}'",
                       boundary_patch.name()));
    }
}

template <typename F, typename G>
void Convection<F, G>::apply_boundary_fixed(const mesh::Cell& cell, const mesh::Face& face) {
    const auto& boundary_patch = _mesh.face_boundary_patch(face);
    auto phi_wall = boundary_patch.get_scalar_bc(_phi.name());
    auto cell_id = cell.id();

    // face area vector
    const auto& S_f = face.area_vector();

    const auto& U_f = boundary_patch.get_vector_bc(_U.name());

    auto m_dot_f = face_mass_flow_rate(_rho, U_f, S_f);

    rhs_vector().coeffRef(cell_id) += -m_dot_f * phi_wall;
}

template <typename F, typename G>
void Convection<F, G>::apply_boundary_outlet(const mesh::Cell& cell, const mesh::Face& face) {
    const auto& boundary_patch = _mesh.face_boundary_patch(face);
    auto phi_Wall = boundary_patch.get_scalar_bc(_phi.name());
    auto cell_id = cell.id();

    // face area vector
    const auto& S_f = face.area_vector();

    const auto& U_f = _U[cell_id];

    auto m_dot_f = face_mass_flow_rate(_rho, U_f, S_f);

    rhs_vector().coeffRef(cell_id) += -m_dot_f * phi_Wall;
}

auto CentralDifference::interpolate(double m_dot,
                                    const mesh::Cell& cell,
                                    const mesh::Cell& neighbor,
                                    const mesh::Face& face) -> InterpolationCoeffs {
    auto d_Cf_C = face.center() - cell.center();
    auto d_Cf_N = face.center() - neighbor.center();

    // calculate gradient of field phi at face centroid

    auto a_C = std::max(m_dot, 0.0);
    auto a_D = std::max(-m_dot, 0.0);

    return {a_C, a_D, b};
}

} // namespace prism::convection