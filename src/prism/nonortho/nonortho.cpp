#include "nonortho.h"

#include "prism/constants.h"
#include "prism/types.h"

namespace prism::nonortho {
auto NoneCorrector::interior_triplet(const mesh::Cell& owner,    // NOLINT
                                     const mesh::Cell& neighbor, // NOLINT
                                     const mesh::Face& face) -> NonOrthoTriplet {
    NonOrthoTriplet triplets;

    triplets.Sf = face.area_vector();
    triplets.Ef = triplets.Sf;
    triplets.Tf = Vector3d {0., 0., 0.};

    return triplets;
}

auto NoneCorrector::boundary_triplet(const mesh::Cell& owner, const mesh::Face& face) // NOLINT
    -> NonOrthoTriplet {
    NonOrthoTriplet triplets;

    triplets.Sf = face.area_vector();
    triplets.Ef = triplets.Sf;
    triplets.Tf = Vector3d {0., 0., 0.};

    return triplets;
}

auto OverRelaxedCorrector::interior_triplet(const mesh::Cell& owner,
                                            const mesh::Cell& neighbor,
                                            const mesh::Face& face) -> NonOrthoTriplet {
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

auto OverRelaxedCorrector::boundary_triplet(const mesh::Cell& owner, const mesh::Face& face)
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