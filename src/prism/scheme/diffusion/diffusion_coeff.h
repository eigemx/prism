#pragma once

#include "prism/field/scalar.h"
#include "prism/field/tensor.h"
#include "prism/mesh/cell.h"
#include "prism/mesh/face.h"
#include "prism/types.h"

namespace prism::scheme::diffusion {

/**
 * @brief Interface for diffusion coefficient in the discretized diffusion equation.
 * @details The diffusion coefficient kappa appears in the term ∇·(κ∇φ). This interface
 * abstracts over scalar (isotropic) and tensor (anisotropic) diffusion coefficients.
 * @see ScalarDiffusionCoeff for isotropic diffusion.
 * @see TensorDiffusionCoeff for anisotropic diffusion.
 */
class IDiffusionCoeff {
  public:
    virtual ~IDiffusionCoeff() = default;
    IDiffusionCoeff() = default;
    IDiffusionCoeff(const IDiffusionCoeff&) = default;
    IDiffusionCoeff(IDiffusionCoeff&&) = default;
    auto operator=(const IDiffusionCoeff&) -> IDiffusionCoeff& = default;
    auto operator=(IDiffusionCoeff&&) -> IDiffusionCoeff& = default;

    /**
     * @brief Multiplies a vector by the diffusion coefficient at a face.
     * @param vector The vector to multiply.
     * @param face The face at which to evaluate the coefficient.
     * @return The multiplied vector.
     */
    virtual auto multiply(const Vector3d& vector, const mesh::Face& face) const -> Vector3d = 0;

    /**
     * @brief Multiplies a vector by the diffusion coefficient at a cell.
     * @param vector The vector to multiply.
     * @param cell The cell at which to evaluate the coefficient.
     * @return The multiplied vector.
     */
    virtual auto multiply(const Vector3d& vector, const mesh::Cell& cell) const -> Vector3d = 0;
};

/**
 * @brief Scalar (isotropic) diffusion coefficient.
 * @details Represents an isotropic diffusion coefficient where κ is a scalar field.
 * The multiply operation returns κ_f * vector where κ_f is the scalar value at the face/cell.
 */
class ScalarDiffusionCoeff : public IDiffusionCoeff {
  public:
    /**
     * @brief Constructs a scalar diffusion coefficient.
     * @param scalar The scalar field representing the diffusion coefficient.
     */
    explicit ScalarDiffusionCoeff(SharedPtr<field::Scalar> scalar);

    auto multiply(const Vector3d& vector, const mesh::Face& face) const -> Vector3d override;
    auto multiply(const Vector3d& vector, const mesh::Cell& cell) const -> Vector3d override;

  private:
    SharedPtr<field::Scalar> _scalar;
};

/**
 * @brief Tensor (anisotropic) diffusion coefficient.
 * @details Represents an anisotropic diffusion coefficient where κ is a tensor field.
 * The multiply operation returns κ_f^T * vector where κ_f is the tensor value at the face/cell.
 * The transpose is taken to follow equation (8.93) from Moukalled et al. (2015).
 */
class TensorDiffusionCoeff : public IDiffusionCoeff {
  public:
    /**
     * @brief Constructs a tensor diffusion coefficient.
     * @param tensor The tensor field representing the diffusion coefficient.
     */
    explicit TensorDiffusionCoeff(SharedPtr<field::Tensor> tensor);

    auto multiply(const Vector3d& vector, const mesh::Face& face) const -> Vector3d override;
    auto multiply(const Vector3d& vector, const mesh::Cell& cell) const -> Vector3d override;

  private:
    SharedPtr<field::Tensor> _tensor;
};

} // namespace prism::scheme::diffusion
