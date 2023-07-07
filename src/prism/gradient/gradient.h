#pragma once

#include <type_traits>

#include "../field.h"
#include "../mesh/pmesh.h"
#include "../types.h"

namespace prism::gradient {

// Base class for gradient schemes for explicity calculating the cell gradient of a scalar field.
// All gradient schemes should inherit from this class and define gradient() function.
class GradientSchemeBase {
  public:
    GradientSchemeBase() = default;
    GradientSchemeBase(const GradientSchemeBase& g) = default;
    GradientSchemeBase(GradientSchemeBase&& g) = default;
    auto operator=(GradientSchemeBase&& g) -> GradientSchemeBase& = default;
    auto operator=(const GradientSchemeBase& g) -> GradientSchemeBase& = default;
    virtual ~GradientSchemeBase() = default;

    virtual auto gradient_at_cell(const mesh::Cell& c) -> Vector3d = 0;
    virtual auto gradient_at_face(const mesh::Face& f) -> Vector3d = 0;
    virtual auto gradient_field() -> VectorField = 0;
};

class GreenGauss : public GradientSchemeBase {
  public:
    GreenGauss(const ScalarField& field);
    auto gradient_at_cell(const mesh::Cell& cell) -> Vector3d override;
    auto gradient_at_face(const mesh::Face& face) -> Vector3d override;
    auto gradient_field() -> VectorField override;

  private:
    const ScalarField& _field;
    MatrixX3d _cell_gradients;
};


template <typename G>
auto create(const ScalarField& field)
    -> std::enable_if_t<std::is_base_of_v<GradientSchemeBase, G>, std::shared_ptr<G>> {
    /**
     * @brief Creates a gradient scheme for a scalar field, and returns a shared pointer to it.
     *
     * @tparam G Gradient scheme type (must inherit from GradientSchemeBase)
     * @param field Scalar field
     * @return std::shared_ptr<G> Gradient scheme
     */
    return std::make_shared<G>(field);
}

} // namespace prism::gradient