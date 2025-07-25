#pragma once

#include "prism/field/ifield.h"

namespace prism {
void exportToVTU(const field::IScalar& field, const std::string& file_name);
}
