#pragma once

#include <memory>
#include <string>

#include "diffusion/diffusion.h"
#include "mesh/pmesh.h"
#include "types.h"

namespace prism {

class LinearSystem {
  public:
    virtual auto coeff_matrix() const -> const SparseMatrix& { return _coeff_matrix; }
    virtual auto coeff_matrix() -> SparseMatrix& { return _coeff_matrix; }

    virtual auto lhs_vector() const -> const VectorXd& { return _b; }
    virtual auto lhs_vector() -> VectorXd& { return _b; }

    virtual void update_coeffs() = 0;

  private:
    SparseMatrix _coeff_matrix;
    VectorXd _b;
};

class SchemeCollector {
  public:
    template <typename Scheme>
    void add_scheme(Scheme& scheme) {
        _schemes.emplace_back(std::make_shared<Scheme>(scheme));
    }
    auto schemes() const -> const std::vector<std::shared_ptr<FVScheme>>& { return _schemes; }

  private:
    std::vector<std::shared_ptr<FVScheme>> _schemes {};
};

class SteadyConservedScalar : public LinearSystem, public SchemeCollector {
  public:
    SteadyConservedScalar(std::string scalar_name, mesh::PMesh& mesh);
    void update_coeffs() override;

  private:
    mesh::PMesh& _mesh;
    std::string _scalar_name;
};
} // namespace prism