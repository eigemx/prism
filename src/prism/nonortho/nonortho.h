#pragma once

#include "../mesh/pmesh.h"
#include "prism/mesh/face.h"

namespace prism::nonortho {

struct NonOrthoTriplets {
    Vector3d Sf {}, Ef {}, Tf {};
};

class Corrector {
  public:
    Corrector() = default;
    virtual ~Corrector() = default;
    Corrector(const Corrector&) = default;
    Corrector(Corrector&&) = default;
    auto operator=(const Corrector&) -> Corrector& = default;
    auto operator=(Corrector&&) -> Corrector& = default;

    virtual auto triplets_interior(const mesh::Cell& owner,
                                   const mesh::Cell& neighbor,
                                   const mesh::Face& face) -> NonOrthoTriplets = 0;

    virtual auto triplets_boundary(const mesh::Cell& owner, const mesh::Face& face)
        -> NonOrthoTriplets = 0;
};

class NoneCorrector : public Corrector {
  public:
    auto triplets_interior(const mesh::Cell& owner,
                           const mesh::Cell& neighbor,
                           const mesh::Face& face) -> NonOrthoTriplets override;

    auto triplets_boundary(const mesh::Cell& owner, const mesh::Face& face)
        -> NonOrthoTriplets override;
};

class OverRelaxed : public Corrector {
  public:
    auto triplets_interior(const mesh::Cell& owner,
                           const mesh::Cell& neighbor,
                           const mesh::Face& face) -> NonOrthoTriplets override;

    auto triplets_boundary(const mesh::Cell& owner, const mesh::Face& face)
        -> NonOrthoTriplets override;
};

} // namespace prism::nonortho