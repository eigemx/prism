#pragma once

#include <stdexcept>

#include "../field.h"
#include "prism/fvscheme.h"
#include "prism/mesh/cell.h"
#include "prism/mesh/face.h"
#include "source.h"

namespace prism::source {

template <SourceSign Sign = SourceSign::Positive>
class Divergence : public FVScheme {
  public:
    Divergence(VectorField U) : _U(U), FVScheme(U.mesh().n_cells()) {}

    void apply() override;
    auto field() -> ScalarField& override;

  private:
    void inline apply_interior(const mesh::Face& face) override {}
    void inline apply_boundary(const mesh::Cell& cell, const mesh::Face& face) override {}
    VectorField _U;
};

template <SourceSign Sign>
void Divergence<Sign>::apply() {
    const auto& mesh = _U.mesh();
    for (const auto& cell : mesh.cells()) {
        double div_sum = 0.0;

        for (auto face_id : cell.faces_ids()) {
            const auto& face = mesh.face(face_id);
        }
    }
}

} // namespace prism::source
