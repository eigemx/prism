#pragma once

#include <cassert>

#include "diffusion_boundary.h"
#include "diffusion_coeff.h"
#include "nonortho.h"
#include "prism/boundary.h"
#include "prism/field/scalar.h"
#include "prism/mesh/cell.h"
#include "prism/mesh/face.h"
#include "prism/mesh/utilities.h"
#include "prism/scheme/boundary.h"
#include "prism/scheme/scheme.h"
#include "prism/types.h"


namespace prism::scheme::diffusion {
// Base class for all diffusion schemes with shared methods (apply(), field(), and kappa()).
class IDiffusion : public IFullScheme {
  public:
    IDiffusion(SharedPtr<IDiffusionCoeff> kappa, const SharedPtr<field::Scalar>& phi);
    auto kappa() -> const SharedPtr<IDiffusionCoeff>& { return _kappa; }

  private:
    SharedPtr<IDiffusionCoeff> _kappa;
};

// Base class for all diffusion schemes that do not need non-orthogonal correction.
class INonCorrected : public IDiffusion {
  public:
    using IDiffusion::IDiffusion;
};

// Base class for all diffusion schemes that need non-orthogonal correction.
class ICorrected : public IDiffusion {
  public:
    using IDiffusion::IDiffusion;
};

// Diffusion (laplacian) scheme without non-orthogonal correction.
class NonCorrected
    : public INonCorrected,
      public prism::boundary::BHManagerProvider<boundary::ISchemeBoundaryHandler<NonCorrected>> {
  public:
    template <typename Kappa>
    NonCorrected(const SharedPtr<Kappa>& kappa, const SharedPtr<field::Scalar>& phi);
    auto needsCorrection() const noexcept -> bool override { return false; }

  private:
    void applyInterior(const mesh::Face& face) override;
    void applyBoundary() override;
};

// Diffusion (laplacian) scheme with non-orthogonal correction.
class Corrected
    : public ICorrected,
      public prism::boundary::BHManagerProvider<boundary::ISchemeBoundaryHandler<Corrected>> {
  public:
    template <typename Kappa>
    Corrected(const SharedPtr<Kappa>& kappa,
              const SharedPtr<field::Scalar>& phi,
              SharedPtr<nonortho::INonOrthoCorrector> corrector =
                  std::make_shared<nonortho::OverRelaxedCorrector>());

    auto corrector() const noexcept -> const nonortho::INonOrthoCorrector& { return *_corrector; }
    auto needsCorrection() const noexcept -> bool override { return true; }

    template <typename Corrector>
    auto setCorrector() -> void {
        _corrector = std::make_shared<Corrector>();
    }

  private:
    void applyInterior(const mesh::Face& face) override;
    void applyBoundary() override;

    SharedPtr<nonortho::INonOrthoCorrector> _corrector;
};

//
// diffusion::IDiffusion implementation
//
inline IDiffusion::IDiffusion(SharedPtr<IDiffusionCoeff> kappa, const SharedPtr<field::Scalar>& phi)
    : _kappa(std::move(kappa)), IFullScheme(phi) {}

namespace detail {
template <typename Kappa>
auto wrapKappa(SharedPtr<Kappa> kappa) -> SharedPtr<IDiffusionCoeff> {
    if constexpr (std::is_same_v<Kappa, field::Scalar>) {
        return std::make_shared<ScalarDiffusionCoeff>(kappa);
    } else {
        return std::make_shared<TensorDiffusionCoeff>(kappa);
    }
}
} // namespace detail

///
/// diffusion::NonCorrected implementation
///
template <typename Kappa>
NonCorrected::NonCorrected(const SharedPtr<Kappa>& kappa, const SharedPtr<field::Scalar>& phi)
    : INonCorrected(detail::wrapKappa(kappa), phi) {
    this->boundaryHandlersManager().template addHandler<boundary::Fixed<NonCorrected>>();
    this->boundaryHandlersManager().template addHandler<boundary::Symmetry<NonCorrected>>();
    this->boundaryHandlersManager().template addHandler<boundary::Outlet<NonCorrected>>();
    this->boundaryHandlersManager().template addHandler<boundary::FixedGradient<NonCorrected>>();
    this->boundaryHandlersManager().template addHandler<boundary::ZeroGradient<NonCorrected>>();
    this->boundaryHandlersManager().template addHandler<boundary::NoSlip<NonCorrected>>();
}

inline void NonCorrected::applyBoundary() {
    prism::boundary::detail::applyBoundary("prism::scheme::diffusion::NonCorrected", *this);
}

inline void NonCorrected::applyInterior(const mesh::Face& face) {
    const auto& mesh = this->field()->mesh();
    const mesh::Cell& owner = mesh->cell(face.owner());
    const mesh::Cell& neighbor = mesh->cell(face.neighbor().value());

    // vector joining the centers of the two cells
    auto d_CF = neighbor.center() - owner.center();
    auto d_CF_norm = d_CF.norm();

    const Vector3d& Sf = face.areaVector();

    // The following is based on equation (8.93) from Moukalled et al. (2015).
    // this handles the general case when kappa is a field::Tensor.
    Vector3d Sf_prime = this->kappa()->multiply(Sf, face);

    // geometric diffusion coefficient
    // Taking the norm of Sf_prime discards the sign of kappa, so we use the following instead
    const f64 g_diff = Sf_prime.dot(d_CF) / (d_CF_norm * d_CF_norm + EPSILON);

    const std::size_t owner_id = owner.id();
    const std::size_t neighbor_id = neighbor.id();

    // g_diff * (Φ_C - Φ_N)
    // diagonal coefficients
    // since we are applying the discretized equation per face basis, not per cell basis,
    // we need to insert the coefficients for both owner and neighbor cells when roles are
    // reversed.
    this->insert(owner_id, owner_id, g_diff);
    this->insert(neighbor_id, neighbor_id, g_diff);

    // off-diagonal coefficients
    this->insert(owner_id, neighbor_id, -g_diff);
    this->insert(neighbor_id, owner_id, -g_diff);
}

//
// diffusion::Corrected implementation
//
template <typename Kappa>
Corrected::Corrected(const SharedPtr<Kappa>& kappa,
                     const SharedPtr<field::Scalar>& phi,
                     SharedPtr<nonortho::INonOrthoCorrector> corrector)
    : ICorrected(detail::wrapKappa(kappa), phi), _corrector(std::move(corrector)) {
    this->boundaryHandlersManager().template addHandler<boundary::Fixed<Corrected>>();
    this->boundaryHandlersManager().template addHandler<boundary::Symmetry<Corrected>>();
    this->boundaryHandlersManager().template addHandler<boundary::Outlet<Corrected>>();
    this->boundaryHandlersManager().template addHandler<boundary::FixedGradient<Corrected>>();
    this->boundaryHandlersManager().template addHandler<boundary::ZeroGradient<Corrected>>();
    this->boundaryHandlersManager().template addHandler<boundary::NoSlip<Corrected>>();
}

inline void Corrected::applyBoundary() {
    prism::boundary::detail::applyBoundary("prism::scheme::diffusion::Corrected", *this);
}

inline void Corrected::applyInterior(const mesh::Face& face) {
    const auto& mesh = this->field()->mesh();
    const mesh::Cell& owner = mesh->cell(face.owner());
    const mesh::Cell& neighbor = mesh->cell(face.neighbor().value());
    const std::size_t owner_id = owner.id();
    const std::size_t neighbor_id = neighbor.id();

    // vector joining the centers of the two cells
    const Vector3d d_CF = neighbor.center() - owner.center();
    const f64 d_CF_norm = d_CF.norm();
    const Vector3d e = d_CF / d_CF_norm;

    const auto Sf = mesh::outwardAreaVector(face, owner);
    const Vector3d Sf_prime = this->kappa()->multiply(Sf, face);
    const auto [Ef_prime, Tf_prime] = _corrector->decompose(Sf_prime, e);

    // geometric diffusion coefficient
    const f64 g_diff = Ef_prime.dot(d_CF) / (d_CF_norm * d_CF_norm + EPSILON);

    /// g_diff * (Φ_C - Φ_N)
    // diagonal coefficients
    this->insert(owner_id, owner_id, g_diff);
    this->insert(neighbor_id, neighbor_id, g_diff);

    // off-diagonal coefficients
    this->insert(owner_id, neighbor_id, -g_diff);
    this->insert(neighbor_id, owner_id, -g_diff);

    // update right hand side
    // cross-diffusion term is added to the right hand side of the equation
    // check equation 8.80 - Chapter 8 (Moukalled et al., 2015)
    const Vector3d grad_f = this->field()->gradAtFace(face);
    this->rhs(owner_id) += Tf_prime.dot(grad_f);
    this->rhs(neighbor_id) += -Tf_prime.dot(grad_f);
}


} // namespace prism::scheme::diffusion
