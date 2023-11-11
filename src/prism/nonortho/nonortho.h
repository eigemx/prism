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

class ICorrector {
  public:
    ICorrector() = default;
    virtual ~ICorrector() = default;
    ICorrector(const ICorrector&) = default;
    ICorrector(ICorrector&&) noexcept = default;
    auto operator=(const ICorrector&) -> ICorrector& = default;
    auto operator=(ICorrector&&) noexcept -> ICorrector& = default;

    virtual auto interior_triplet(const mesh::Cell& owner,
                                  const mesh::Cell& neighbor,
                                  const mesh::Face& face) -> NonOrthoTriplet = 0;

    virtual auto boundary_triplet(const mesh::Cell& owner, const mesh::Face& face)
        -> NonOrthoTriplet = 0;
};


class NilCorrector : public ICorrector {
  public:
    auto interior_triplet(const mesh::Cell& owner,
                          const mesh::Cell& neighbor,
                          const mesh::Face& face) -> NonOrthoTriplet override;

    auto boundary_triplet(const mesh::Cell& owner, const mesh::Face& face)
        -> NonOrthoTriplet override;
};


template <typename GradScheme = gradient::LeastSquares>
class ICorrectorWithGradScheme : public ICorrector {
  public:
    ICorrectorWithGradScheme(ScalarField& field) : _grad_scheme(field) {}
    auto grad_scheme() -> GradScheme& { return _grad_scheme; }

  private:
    GradScheme _grad_scheme;
};


template <typename GradScheme = gradient::LeastSquares>
class OverRelaxedCorrector : public ICorrectorWithGradScheme<GradScheme> {
  public:
    OverRelaxedCorrector(ScalarField& field) : ICorrectorWithGradScheme<GradScheme>(field) {}
    auto interior_triplet(const mesh::Cell& owner,
                          const mesh::Cell& neighbor,
                          const mesh::Face& face) -> NonOrthoTriplet override;

    auto boundary_triplet(const mesh::Cell& owner, const mesh::Face& face)
        -> NonOrthoTriplet override;
};


auto inline NilCorrector::interior_triplet(const mesh::Cell& owner,    // NOLINT
                                           const mesh::Cell& neighbor, // NOLINT
                                           const mesh::Face& face) -> NonOrthoTriplet {
    NonOrthoTriplet triplets;

    triplets.Sf = face.area_vector();
    triplets.Ef = triplets.Sf;
    triplets.Tf = Vector3d {0., 0., 0.};

    return triplets;
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
auto OverRelaxedCorrector<GradScheme>::interior_triplet(const mesh::Cell& owner,
                                                        const mesh::Cell& neighbor,
                                                        const mesh::Face& face)
    -> NonOrthoTriplet {
    const auto& Sf = face.area_vector();

    // vector joining the centers of the two cells
    auto d_CF = neighbor.center() - owner.center();
    auto d_CF_norm = d_CF.norm();

    // unit vector in d_CF direction
    auto e = d_CF / d_CF_norm;

    // orthogonal-like normal vector E_f using over-relaxed approach
    auto Ef = ((Sf.dot(Sf) / (e.dot(Sf) + EPSILON))) * e;

    return NonOrthoTriplet {Sf, Ef, Sf - Ef};
}


template <typename GradScheme>
auto OverRelaxedCorrector<GradScheme>::boundary_triplet(const mesh::Cell& owner,
                                                        const mesh::Face& face)
    -> NonOrthoTriplet {
    const auto& Sf = face.area_vector();

    // vector joining the centers of the cell and the face
    auto d_Cf = face.center() - owner.center();
    auto d_Cf_norm = d_Cf.norm();
    auto e = d_Cf / d_Cf_norm;

    auto Ef = ((Sf.dot(Sf) / e.dot(Sf))) * e;

    return NonOrthoTriplet {Sf, Ef, Sf - Ef};
}

} // namespace prism::nonortho