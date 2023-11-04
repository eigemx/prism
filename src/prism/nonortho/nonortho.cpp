#include "nonortho.h"

#include "../constants.h"
#include "../types.h"

namespace prism::nonortho {
auto NoneCorrector::triplets_interior(const mesh::Cell& owner,    // NOLINT
                                      const mesh::Cell& neighbor, // NOLINT
                                      const mesh::Face& face) -> NonOrthoTriplets {
    NonOrthoTriplets triplets;

    triplets.Sf = face.area_vector();
    triplets.Ef = triplets.Sf;
    triplets.Tf = Vector3d {0., 0., 0.};

    return triplets;
}

auto NoneCorrector::triplets_boundary(const mesh::Cell& owner, const mesh::Face& face) // NOLINT
    -> NonOrthoTriplets {
    NonOrthoTriplets triplets;

    triplets.Sf = face.area_vector();
    triplets.Ef = triplets.Sf;
    triplets.Tf = Vector3d {0., 0., 0.};

    return triplets;
}

auto OverRelaxed::triplets_interior(const mesh::Cell& owner,
                                    const mesh::Cell& neighbor,
                                    const mesh::Face& face) -> NonOrthoTriplets {
    const auto& Sf = face.area_vector();

    // vector joining the centers of the two cells
    auto d_CF = neighbor.center() - owner.center();
    auto d_CF_norm = d_CF.norm();

    // unit vector in d_CF direction
    auto e = d_CF / d_CF_norm;

    // orthogonal-like normal vector E_f using over-relaxed approach
    auto Ef = ((Sf.dot(Sf) / (e.dot(Sf) + EPSILON))) * e;

    return NonOrthoTriplets {Sf, Ef, Sf - Ef};
}

auto OverRelaxed::triplets_boundary(const mesh::Cell& owner, const mesh::Face& face)
    -> NonOrthoTriplets {
    const auto& Sf = face.area_vector();

    // vector joining the centers of the cell and the face
    auto d_Cf = face.center() - owner.center();
    auto d_Cf_norm = d_Cf.norm();
    auto e = d_Cf / d_Cf_norm;

    auto Ef = ((Sf.dot(Sf) / e.dot(Sf))) * e;

    return NonOrthoTriplets {Sf, Ef, Sf - Ef};
}
} // namespace prism::nonortho