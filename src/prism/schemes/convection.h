#pragma once

#include <spdlog/spdlog.h>

#include <cmath>

#include "convection_boundary.h"
#include "fvscheme.h"
#include "prism/field/scalar.h"
#include "prism/field/velocity.h"
#include "prism/gradient/gradient.h"
#include "prism/mesh/cell.h"
#include "prism/mesh/pmesh.h"
#include "prism/mesh/utilities.h"
#include "prism/schemes/boundary.h"
#include "prism/types.h"

namespace prism::scheme::convection {

namespace detail {
// coefficients for the discretized convection equation for a face
struct CoeffsTriplet {
    double a_C {0.0}; // cell
    double a_N {0.0}; // neighbor
    double b {0.0};   // source
};
} // namespace detail

// Finite volume scheme for the discretization of the convection term
class IConvection : public FVScheme<field::Scalar> {
  public:
    IConvection(field::Scalar rho, field::Velocity U, field::Scalar phi);

    void apply() override;

    auto inline field() -> std::optional<field::Scalar> override { return _phi; }
    auto inline U() -> const field::Velocity& { return _U; }
    auto inline rho() -> const field::Scalar& { return _rho; }

    using BoundaryHandlersManager =
        prism::boundary::BoundaryHandlersManager<IConvection,
                                                 boundary::FVSchemeBoundaryHandler<IConvection>>;
    auto boundaryHandlersManager() -> BoundaryHandlersManager& { return _bc_manager; }


  private:
    virtual auto interpolate(double m_dot,
                             const mesh::Cell& cell,
                             const mesh::Cell& neighbor,
                             const mesh::Face& face) -> detail::CoeffsTriplet = 0;

    void apply_interior(const mesh::Face& face) override;
    void apply_boundary(const mesh::Face& face) override {}
    void apply_boundary();

    field::Scalar _rho;
    field::Velocity _U;
    field::Scalar _phi;
    BoundaryHandlersManager _bc_manager;
};

// Central difference scheme
template <typename G = gradient::LeastSquares>
class CentralDifference : public IConvection, public gradient::GradientProvider<G> {
  public:
    CentralDifference(field::Scalar& rho, field::Velocity& U, field::Scalar& phi)
        : IConvection(rho, U, phi) {}

  private:
    auto interpolate(double m_dot,
                     const mesh::Cell& cell,
                     const mesh::Cell& neighbor,
                     const mesh::Face& face) -> detail::CoeffsTriplet override;
};


// Upwind scheme
template <typename G = gradient::LeastSquares>
class Upwind : public IConvection, public gradient::GradientProvider<G> {
  public:
    Upwind(field::Scalar rho, field::Velocity U, field::Scalar phi)
        : IConvection(rho, U, phi), gradient::GradientProvider<G>(phi) {}

  private:
    auto interpolate(double m_dot,
                     const mesh::Cell& cell,
                     const mesh::Cell& neighbor,
                     const mesh::Face& face) -> detail::CoeffsTriplet override;
};


// Second order upwind scheme
template <typename G = gradient::LeastSquares>
class SecondOrderUpwind : public IConvection, public gradient::GradientProvider<G> {
  public:
    SecondOrderUpwind(field::Scalar rho, field::Velocity U, field::Scalar phi)
        : IConvection(rho, U, phi), gradient::GradientProvider<G>(phi) {}

  private:
    auto interpolate(double m_dot,
                     const mesh::Cell& cell,
                     const mesh::Cell& neighbor,
                     const mesh::Face& face) -> detail::CoeffsTriplet override;
};


// QUICK scheme
template <typename G = gradient::LeastSquares>
class QUICK : public IConvection, public gradient::GradientProvider<G> {
  public:
    QUICK(field::Scalar rho, field::Velocity U, field::Scalar phi)
        : IConvection(rho, U, phi), gradient::GradientProvider<G>(phi) {}

  private:
    auto interpolate(double m_dot,
                     const mesh::Cell& cell,
                     const mesh::Cell& neighbor,
                     const mesh::Face& face) -> detail::CoeffsTriplet override;
};

template <typename G>
auto CentralDifference<G>::interpolate(double m_dot,
                                       const mesh::Cell& cell,
                                       const mesh::Cell& neighbor,
                                       const mesh::Face& face) -> detail::CoeffsTriplet {
    // in case `cell` is the upstream cell
    const Vector3d face_grad_phi = this->grad_scheme().gradient_at_face(face);
    const Vector3d d_Cf = face.center() - cell.center();
    const double a_C = std::max(m_dot, 0.0);
    double b = -std::max(m_dot, 0.0) * d_Cf.dot(face_grad_phi);

    // in case 'neighbor' is the upstream cell
    const Vector3d d_Nf = face.center() - neighbor.center();
    const double a_N = -std::max(-m_dot, 0.0);
    b += std::max(-m_dot, 0.0) * d_Nf.dot(face_grad_phi);

    return {a_C, a_N, b};
}

template <typename G>
auto Upwind<G>::interpolate(double m_dot,
                            const mesh::Cell& cell,     // NOLINT
                            const mesh::Cell& neighbor, // NOLINT
                            const mesh::Face& face)     // NOLINT
    -> detail::CoeffsTriplet {
    // in case `cell` is the upstream cell
    const double a_C = std::max(m_dot, 0.0);
    // in case 'neighbor' is the upstream cell
    const double a_N = -std::max(-m_dot, 0.0);

    return {a_C, a_N, 0.0};
}

template <typename G>
auto SecondOrderUpwind<G>::interpolate(double m_dot,
                                       const mesh::Cell& cell,
                                       const mesh::Cell& neighbor,
                                       const mesh::Face& face) -> detail::CoeffsTriplet {
    // in case `cell` is the upstream cell
    const Vector3d face_grad_phi = this->grad_scheme().gradAtFace(face);
    const Vector3d cell_grad_phi = this->grad_scheme().gradAtCell(cell);
    const Vector3d neighbor_grad_phi = this->grad_scheme().gradAtCell(neighbor);

    const Vector3d d_Cf = face.center() - cell.center();
    // auto correction = d_Cf.dot((2 * cell_grad_phi) - face_grad_phi);
    auto correction = cell_grad_phi.dot(d_Cf);

    const double a_C = std::max(m_dot, 0.0);
    const double b1 = -std::max(m_dot, 0.0) * correction;

    // in case 'neighbor' is the upstream cell
    const Vector3d d_Nf = face.center() - neighbor.center();
    // correction = d_Nf.dot((2 * neighbor_grad_phi) - face_grad_phi);
    correction = neighbor_grad_phi.dot(d_Nf);

    const double a_N = -std::max(-m_dot, 0.0);
    const double b2 = std::max(-m_dot, 0.0) * correction;

    return {a_C, a_N, b1 + b2};
}

template <typename G>
auto QUICK<G>::interpolate(double m_dot,
                           const mesh::Cell& cell,
                           const mesh::Cell& neighbor,
                           const mesh::Face& face) -> detail::CoeffsTriplet {
    // in case `cell` is the upstream cell
    const Vector3d face_grad_phi = this->grad_scheme().gradAtFace(face);
    const Vector3d cell_grad_phi = this->grad_scheme().gradAtCell(cell);
    const Vector3d neighbor_grad_phi = this->grad_scheme().gradAtCell(neighbor);

    const Vector3d d_Cf = face.center() - cell.center();
    auto correction = 0.5 * d_Cf.dot(cell_grad_phi + face_grad_phi);

    const double a_C = std::max(m_dot, 0.0);
    const double b1 = -std::max(m_dot, 0.0) * correction;

    // in case 'neighbor' is the upstream cell
    const Vector3d d_Nf = face.center() - neighbor.center();
    correction = 0.5 * d_Nf.dot(neighbor_grad_phi + face_grad_phi);

    const double a_N = -std::max(-m_dot, 0.0);
    const double b2 = std::max(-m_dot, 0.0) * correction;

    return {a_C, a_N, b1 + b2};
}


} // namespace prism::scheme::convection
