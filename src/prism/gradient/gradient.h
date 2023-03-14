#pragma once

#include <string>

#include "../mesh/pmesh.h"
#include "../types.h"

namespace prism::gradient {
class GradientSchemeBase {
  public:
    virtual auto gradient(const mesh::Cell& c,
                          const VectorXd& scalar_field,
                          const mesh::PMesh& mesh,
                          const std::string& scalar_name) -> Vector3d = 0;
};

class GreenGauss : public GradientSchemeBase {
  public:
    auto gradient(const mesh::Cell& c,
                  const VectorXd& scalar_field,
                  const mesh::PMesh& mesh,
                  const std::string& scalar_name) -> Vector3d override;
};
} // namespace prism::gradient