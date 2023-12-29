#pragma once

#include "prism/field.h"
#include "prism/gradient/gradient.h"

namespace prism::ops {

template <typename GradScheme = gradient::LeastSquares>
void rhie_chow_correct(field::Vector& U, const field::Tensor& D, const field::Scalar& P);

void rhie_chow_correct(field::Vector& U,
                       const field::Tensor& D,
                       gradient::AbstractGradient* p_grad_scheme);

template <typename GradScheme>
void inline rhie_chow_correct(field::Vector& U, const field::Tensor& D, const field::Scalar& P) {
    GradScheme p_grad_scheme(P);
    return rhie_chow_correct(U, D, &p_grad_scheme);
}

} // namespace prism::ops