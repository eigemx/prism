#pragma once

#include "prism/field/field.h"

namespace prism {
void export_field_vtu(const field::Scalar& field, const std::string& file_name);
}