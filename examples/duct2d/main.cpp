#include <prism/prism.h>

#include "fmt/core.h"
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
        TransportEquation(convection::Upwind(rho, U, U.x()), // ∇.(ρUu)
                          diffusion::Diffusion(1e-6, U.x()), // - ∇.(μ∇u)
                          source::Gradient<source::SourceSign::Negative>(P, Coord::X), // ∂p/∂x
                          source::Laplacian(1e-6, U.x()) // - ∇.(μ∇u^T)
        );

    auto vEqn = TransportEquation(convection::Upwind(rho, U, U.y()),
                                  diffusion::Diffusion(1e-6, U.y()),
                                  source::Gradient<source::SourceSign::Negative>(P, Coord::Y),
                                  source::Laplacian(1e-6, U.y()));


    uEqn.update_coeffs();
    vEqn.update_coeffs();

    // calculate coefficients for the pressure equation
    const auto vol_vec = mesh::cells_volume_vec(mesh);
    const auto& uEqn_diag = uEqn.matrix().diagonal();
    const auto& vEqn_diag = vEqn.matrix().diagonal();

    auto g = ConstScalarField("gravity", mesh, 9.81);

    fmt::println("Value at cell[10] = {}", g.value_at_cell(10));
    auto iface = *(mesh.interior_faces().begin());
    fmt::println("Value at face[{}] = {}", iface.id(), g.value_at_face(iface));
}