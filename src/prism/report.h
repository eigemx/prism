#pragma once

#include "prism/types.h"

namespace prism::report {

struct Entry {
    std::string field_name;
    size_t n_iterations;
    f64 initial_residual;
    f64 final_residual;
    bool converged;
};

} // namespace prism::report
