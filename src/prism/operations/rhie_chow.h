#pragma once

#include "prism/field/pressure.h"
#include "prism/field/tensor.h"
#include "prism/field/vector.h"
#include "prism/gradient/gradient.h"

namespace prism::ops {

template <typename GradScheme = gradient::LeastSquares>
void rhie_chow_correct(field::Vector& U, const field::Tensor& D, const field::Pressure& P);


} // namespace prism::ops