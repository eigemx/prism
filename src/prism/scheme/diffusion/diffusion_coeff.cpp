#include "diffusion_coeff.h"

namespace prism::scheme::diffusion {

ScalarDiffusionCoeff::ScalarDiffusionCoeff(SharedPtr<field::Scalar> scalar)
    : _scalar(std::move(scalar)) {}

auto ScalarDiffusionCoeff::multiply(const Vector3d& vector, const mesh::Face& face) const
    -> Vector3d {
    return _scalar->valueAtFace(face) * vector;
}

auto ScalarDiffusionCoeff::multiply(const Vector3d& vector, const mesh::Cell& cell) const
    -> Vector3d {
    return _scalar->valueAtCell(cell) * vector;
}

TensorDiffusionCoeff::TensorDiffusionCoeff(SharedPtr<field::Tensor> tensor)
    : _tensor(std::move(tensor)) {}

auto TensorDiffusionCoeff::multiply(const Vector3d& vector, const mesh::Face& face) const
    -> Vector3d {
    return _tensor->valueAtFace(face).transpose() * vector;
}

auto TensorDiffusionCoeff::multiply(const Vector3d& vector, const mesh::Cell& cell) const
    -> Vector3d {
    // The transpose is to follow equation (8.93) from Moukalled et al. (2015)
    return _tensor->valueAtCell(cell).transpose() * vector;
}

} // namespace prism::scheme::diffusion
