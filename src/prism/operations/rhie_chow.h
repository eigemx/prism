#pragma once

#include "prism/field/field.h"
#include "prism/field/pressure.h"
#include "prism/gradient/gradient.h"

namespace prism::ops {

template <typename GradScheme = gradient::LeastSquares>
void rhie_chow_correct(field::Vector& U, const field::Tensor& D, const field::Pressure& P);


} // namespace prism::ops