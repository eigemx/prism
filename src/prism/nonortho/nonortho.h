#pragma once

#include "prism/constants.h"
#include "prism/field.h"
#include "prism/gradient/gradient.h"
#include "prism/mesh/pmesh.h"
#include "prism/types.h"

namespace prism::nonortho {

struct NonOrthoTriplet {
    Vector3d Sf {}, Ef {}, Tf {};
};

class AbstractCorrector {
  public:
    AbstractCorrector() = default;
    virtual ~AbstractCorrector() = default;
    AbstractCorrector(const AbstractCorrector&) = default;
    AbstractCorrector(AbstractCorrector&&) noexcept = default;
    auto operator=(const AbstractCorrector&) -> AbstractCorrector& = default;
    auto operator=(AbstractCorrector&&) noexcept -> AbstractCorrector& = default;

    virtual auto interior_triplet(const mesh::Cell& owner,
                                  const mesh::Cell& neighbor,
                                  const mesh::Face& face) -> NonOrthoTriplet = 0;

    virtual auto boundary_triplet(const mesh::Cell& owner, const mesh::Face& face)
        -> NonOrthoTriplet = 0;
};

class NilCorrector : public AbstractCorrector {
  public:
    auto interior_triplet(const mesh::Cell& owner,
                          const mesh::Cell& neighbor,
                          const mesh::Face& face) -> NonOrthoTriplet override;

    auto boundary_triplet(const mesh::Cell& owner, const mesh::Face& face)
        -> NonOrthoTriplet override;
};

template <typename GradScheme = gradient::LeastSquares>
class AbstractCorrectorWithGradient : public AbstractCorrector {
  public:
    AbstractCorrectorWithGradient(ScalarField& field) : _grad_scheme(field) {}
    auto grad_scheme() -> GradScheme& { return _grad_scheme; }

    auto interior_triplet(const mesh::Cell& owner,
                          const mesh::Cell& neighbor,
                          const mesh::Face& face) -> NonOrthoTriplet override;

    auto boundary_triplet(const mesh::Cell& owner, const mesh::Face& face)
        -> NonOrthoTriplet override;


  private:
    virtual auto decompose(const Vector3d& Sf, const Vector3d& e) -> Vector3d = 0;
    GradScheme _grad_scheme;
};

template <typename GradScheme = gradient::LeastSquares>
class OverRelaxedCorrector : public AbstractCorrectorWithGradient<GradScheme> {
  public:
    OverRelaxedCorrector(ScalarField& field) : AbstractCorrectorWithGradient<GradScheme>(field) {}

  private:
    auto decompose(const Vector3d& Sf, const Vector3d& e) -> Vector3d override;
};

template <typename GradScheme = gradient::LeastSquares>
class MinimumCorrector : public AbstractCorrectorWithGradient<GradScheme> {
  public:
    MinimumCorrector(ScalarField& field) : AbstractCorrectorWithGradient<GradScheme>(field) {}

  private:
    auto decompose(const Vector3d& Sf, const Vector3d& e) -> Vector3d override;
};

template <typename GradScheme = gradient::LeastSquares>
class OrthogonalCorrector : public AbstractCorrectorWithGradient<GradScheme> {
  public:
    OrthogonalCorrector(ScalarField& field) : AbstractCorrectorWithGradient<GradScheme>(field) {}

  private:
    auto decompose(const Vector3d& Sf, const Vector3d& e) -> Vector3d override;
};

auto inline NilCorrector::interior_triplet(const mesh::Cell& owner,    // NOLINT
                                           const mesh::Cell& neighbor, // NOLINT
                                           const mesh::Face& face) -> NonOrthoTriplet {
    return boundary_triplet(owner, face);
}

auto inline NilCorrector::boundary_triplet(const mesh::Cell& owner, // NOLINT
                                           const mesh::Face& face)  // NOLINT
    -> NonOrthoTriplet {
    NonOrthoTriplet triplets;

    triplets.Sf = face.area_vector();
    triplets.Ef = triplets.Sf;
    triplets.Tf = Vector3d {0., 0., 0.};

    return triplets;
}

template <typename GradScheme>
auto AbstractCorrectorWithGradient<GradScheme>::interior_triplet(const mesh::Cell& owner,
                                                                 const mesh::Cell& neighbor,
                                                                 const mesh::Face& face)
    -> NonOrthoTriplet {
    const auto& Sf = face.area_vector();

    // vector joining the centers of the two cells
    auto d_CF = neighbor.center() - owner.center();
    auto d_CF_norm = d_CF.norm();

    // unit vector in d_CF direction
    auto e = d_CF / d_CF_norm;

    auto Ef = decompose(Sf, e);

    return {Sf, Ef, Sf - Ef};
}

template <typename GradScheme>
auto AbstractCorrectorWithGradient<GradScheme>::boundary_triplet(const mesh::Cell& owner,
                                                                 const mesh::Face& face)
    -> NonOrthoTriplet {
    const auto& Sf = face.area_vector();

    // vector joining the centers of the cell and the face
    auto d_Cf = face.center() - owner.center();
    auto d_Cf_norm = d_Cf.norm();
    auto e = d_Cf / d_Cf_norm;

    auto Ef = decompose(Sf, e);

    return {Sf, Ef, Sf - Ef};
}

template <typename GradScheme>
auto inline OverRelaxedCorrector<GradScheme>::decompose(const Vector3d& Sf, const Vector3d& e)
    -> Vector3d {
    return ((Sf.dot(Sf) / (e.dot(Sf) + EPSILON))) * e;
}

template <typename GradScheme>
auto inline MinimumCorrector<GradScheme>::decompose(const Vector3d& Sf, const Vector3d& e)
    -> Vector3d {
    return (e.dot(Sf)) * e;
}

template <typename GradScheme>
auto inline OrthogonalCorrector<GradScheme>::decompose(const Vector3d& Sf, const Vector3d& e)
    -> Vector3d {
    return Sf.norm() * e;
}


} // namespace prism::nonortho