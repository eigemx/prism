#pragma once

#include "prism/mesh/pmesh.h"

namespace prism::nonortho {

struct NonOrthoTriplet {
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

    virtual auto interior_triplet(const mesh::Cell& owner,
                                  const mesh::Cell& neighbor,
                                  const mesh::Face& face) -> NonOrthoTriplet = 0;

    virtual auto boundary_triplet(const mesh::Cell& owner, const mesh::Face& face)
        -> NonOrthoTriplet = 0;
};

class NilCorrector : public Corrector {
  public:
    auto interior_triplet(const mesh::Cell& owner,
                          const mesh::Cell& neighbor,
                          const mesh::Face& face) -> NonOrthoTriplet override;

    auto boundary_triplet(const mesh::Cell& owner, const mesh::Face& face)
        -> NonOrthoTriplet override;
};

class OverRelaxedCorrector : public Corrector {
  public:
    auto interior_triplet(const mesh::Cell& owner,
                          const mesh::Cell& neighbor,
                          const mesh::Face& face) -> NonOrthoTriplet override;

    auto boundary_triplet(const mesh::Cell& owner, const mesh::Face& face)
        -> NonOrthoTriplet override;
};

} // namespace prism::nonortho