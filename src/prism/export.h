#pragma once

#include "field.h"

namespace prism {
void export_field_vtu(const ScalarField& field, const std::string& file_name);
}