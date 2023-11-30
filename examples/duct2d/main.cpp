#include <prism/prism.h>

#include <iostream>

#include "fmt/core.h"
#include "prism/constants.h"
#include "prism/equation.h"
#include "prism/field.h"
#include "prism/mesh/unv.h"
#include "prism/mesh/utilities.h"
#include "prism/schemes/convection.h"
#include "prism/schemes/diffusion.h"
#include "prism/schemes/source.h"
#include "prism/types.h"

auto main(int argc, char* argv[]) -> int {
    using namespace prism;
    if (argc < 2) {
        fmt::println("usage: {} [mesh-file]", argv[1]); // NOLINT
        return 1;
    }

    auto mesh = mesh::UnvToPMeshConverter(argv[1]).to_pmesh(); // NOLINT

    auto U = VectorField("velocity", mesh, 0.0);
    auto P = ScalarField("pressure", mesh, 0.0);
    auto rho = ScalarField("density", mesh, 1.18);

    auto uEqn =
        TransportEquation(convection::Upwind(rho, U, U.x()),         // ∇.(ρUu)
                          diffusion::AbstractDiffusion(1e-6, U.x()), // - ∇.(μ∇u)
                          source::Gradient<source::SourceSign::Negative>(P, Coord::X), // ∂p/∂x
                          source::Laplacian(1e-6, U.x()) // - ∇.(μ∇u^T)
        );

    auto vEqn = TransportEquation(convection::Upwind(rho, U, U.y()),
                                  diffusion::AbstractDiffusion(1e-6, U.y()),
                                  source::Gradient<source::SourceSign::Negative>(P, Coord::Y),
                                  source::Laplacian(1e-6, U.y()));


    uEqn.update_coeffs();
    vEqn.update_coeffs();

    // calculate coefficients for the pressure equation
    const auto vol_vec = mesh::cells_volume_vec(mesh);
    const auto& uEqn_diag = uEqn.matrix().diagonal();
    const auto& vEqn_diag = vEqn.matrix().diagonal();

    auto D_data = std::vector<prism::Matrix3d>();
    D_data.reserve(mesh.n_cells());

    auto Du = vol_vec.array() / (uEqn_diag.array() + prism::EPSILON);
    auto Dv = vol_vec.array() / (vEqn_diag.array() + prism::EPSILON);

    for (std::size_t i = 0; i < mesh.n_cells(); ++i) {
        Matrix3d Di;
        Di << Du[i], 0, 0, 0, Dv[i], 0, 0, 0, 0;
        D_data.emplace_back(std::move(Di));
    }

    auto D = prism::TensorField("D", mesh, D_data);
}