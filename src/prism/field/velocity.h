#pragma once

#include "prism/field/field.h"

namespace prism::field {

class VelocityComponent : public Scalar {};
class Velocity : public detail::Vector<VelocityComponent> {};

} // namespace prism::field