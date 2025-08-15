#include "SIMPLE.h"

#include <fmt/format.h>

#include <stdexcept>

#include "prism/equation/transport.h"
#include "prism/log.h"
#include "prism/numerics/solver.h"
#include "prism/operations/rhie_chow.h"
#include "prism/scheme/source/divergence.h"
#include "prism/types.h"


namespace prism::algo {
void solveMomentumImplicitly(SIMPLEParameters params,
                             std::span<eqn::Momentum*> momentum_predictors) {
    // solve momentum equations implicitly
    auto momentum_solver = solver::BiCGSTAB<field::VelocityComponent>();
    log::debug("prism::algo::solveMomentumImplicitly(): solving momentum equations");
    for (auto* eqn : momentum_predictors) {
        eqn->updateCoeffs();
        eqn->relax();
        momentum_solver.solve(*eqn, params.momentum_max_iter, params.momentum_residual);
        eqn->zeroOutCoeffs();
    }
}

void constrainPPrime(field::Pressure& pprime) {
    // we need to reset the _pprime field to zero at the boundaries where a Dirichlet condition is
    // applied.
    VectorXd face_values = VectorXd::Zero(pprime.mesh()->faceCount());

    for (const auto& patch : pprime.mesh()->boundaryPatches()) {
        if (patch.isEmpty()) {
            continue; // skip empty patches
        }

        const auto& bc = patch.getBoundaryCondition("P");
        const auto& handler = pprime.boundaryHandlersManager().getHandler(bc.kindString());

        if (handler->isDirichlet() || handler->name() == "symmetry") {
            for (const auto& face_id : patch.facesIds()) {
                face_values[face_id] = 0.0;
            }
            continue;
        }

        // keep the values at the faces for other boundary conditions
        for (const auto& face_id : patch.facesIds()) {
            const auto& face = pprime.mesh()->face(face_id);
            face_values[face_id] = pprime.valueAtFace(face);
        }
    }
    pprime.setFaceValues(face_values);
}

auto pressureEquationCoeffsTensor(std::span<eqn::Momentum*> momentum_predictors,
                                  const field::Pressure& p) -> field::Tensor {
    for (auto* eqn : momentum_predictors) {
        eqn->updateCoeffs();
        eqn->relax();
    }

    // calculate coefficients for the pressure equation
    const auto& mesh = p.mesh();
    const VectorXd& vol_vec = mesh->cellsVolumeVector();
    const VectorXd& uEqn_diag = momentum_predictors[0]->matrix().diagonal();
    const VectorXd& vEqn_diag = momentum_predictors[1]->matrix().diagonal();
    auto D_data = std::vector<Matrix3d>(mesh->cellCount(), Matrix3d::Zero());

    VectorXd Du = vol_vec.array() / (uEqn_diag.array() + EPSILON);
    VectorXd Dv = vol_vec.array() / (vEqn_diag.array() + EPSILON);
    VectorXd Dw = VectorXd::Ones(mesh->cellCount());

    if (momentum_predictors.size() == 3) {
        // 3D case
        const VectorXd& wEqn_diag = momentum_predictors[2]->matrix().diagonal();
        Dw = vol_vec.array() / (wEqn_diag.array() + EPSILON);
    }

    for (const auto& cell : mesh->cells()) {
        auto i = cell.id();
        // clang-format off
            D_data[i]  << Du[i], 0,     0,
                          0,     Dv[i], 0,
                          0,     0,     Dw[i];
        // clang-format on
    }

    for (auto* eqn : momentum_predictors) {
        eqn->zeroOutCoeffs();
    }

    return {"D", mesh, std::move(D_data)};
}

auto solvePressureEquation(SIMPLEParameters params,
                           std::span<eqn::Momentum*> momentum_predictors,
                           field::Velocity& U,
                           field::Velocity& mdot,
                           const field::Pressure& p)
    -> std::pair<field::Pressure, field::Tensor> {
    auto D = pressureEquationCoeffsTensor(momentum_predictors, p);

    // Rhie-Chow interpolation for velocity face values
    log::debug(
        "prism::algo::solvePressureEquation(): applying Rhie-Chow correction to faces velocity");
    // first, we update convective flux with latest values of velocity field before applying the
    // Rhie-Chow correction.
    mdot.updateFaces([&](const mesh::Face& face) { return U.valueAtFace(face); });
    mdot.updateInteriorFaces(
        [&](const mesh::Face& face) { return ops::rhieChowCorrectFace(face, U, D, p); });

    // pressure correction field created with same name as pressure field to get same boundary
    // conditions without having to define _pprime in fields.json file.
    auto pprime = field::Pressure(p.name(), p.mesh(), 0.0);

    // The corrector should reset to zero the correction field at every iteration and should
    // also apply a zero value at all boundaries for which a Dirichlet (fixed) boundary
    // condition is used for the pressure.
    constrainPPrime(pprime);

    /// TODO: based on number of non-orhogonal corrections in _params, we should check if we need
    /// diffusion::Corrected or diffusion::NonCorrected
    using laplacian_p =
        scheme::diffusion::Corrected<field::Tensor,
                                     scheme::diffusion::nonortho::OverRelaxedCorrector,
                                     field::Pressure>;
    using div_U = scheme::source::Divergence<Sign::Negative, field::Velocity>;

    auto pEqn = eqn::Transport<field::Pressure>(laplacian_p(D, pprime), // - ∇.(D ∇P')
                                                div_U(mdot)             // == - (∇.U)
    );
    log::debug("prism::algo::solvePressureEquation(): solving pressure correction equation");
    auto p_solver = solver::BiCGSTAB<field::Pressure>();
    p_solver.solve(pEqn, params.pressure_max_iter, params.pressure_residual);

    // non-orthogonal correction
    for (auto i = 0; i < params.non_ortho_correctors; ++i) {
        p_solver.solve(pEqn, params.pressure_max_iter, params.pressure_residual);
    }

    return {pprime, D};
}

IncompressibleSIMPLE::IncompressibleSIMPLE(SIMPLEParameters parameters) : _params(parameters) {}

void IncompressibleSIMPLE::step(std::span<eqn::Momentum*> momentum_predictors,
                                field::Velocity& U,
                                field::Velocity& mdot,
                                field::Pressure& p) {
    if (momentum_predictors.size() != 2 && momentum_predictors.size() != 3) {
        throw std::runtime_error(
            fmt::format("prism::algo::IncompressibleSIMPLE::step() expects 2 or 3 momentum "
                        "predictors, not {}",
                        momentum_predictors.size()));
    }

    // solve momentum equations
    auto momentum_solver = solver::BiCGSTAB<field::VelocityComponent>();
    log::debug("prism::algo::IncompressibleSIMPLE::step(): solving momentum equations");
    for (auto* eqn : momentum_predictors) {
        momentum_solver.solve(*eqn, _params.momentum_max_iter, _params.momentum_residual);
    }
    auto [pprime, D] = solvePressureEquation(_params, momentum_predictors, U, mdot, p);

    correctFields(U, mdot, p, D, pprime, _params.pressure_urf);
}

void correctFields(field::Velocity& U,
                   field::Velocity& mdot,
                   field::Pressure& p,
                   const field::Tensor& D,
                   const field::Pressure& pprime,
                   double pressure_urf) {
    // update velocity field
    U.updateCells([&](const mesh::Cell& cell) {
        return U.valueAtCell(cell) - (D.valueAtCell(cell) * pprime.gradAtCell(cell));
    });

    // update mass flow rate at faces
    mdot.updateInteriorFaces([&](const mesh::Face& face) {
        return mdot.valueAtFace(face) - (D.valueAtFace(face) * pprime.gradAtFace(face));
    });

    // update pressure
    p.update(p.values() + (pressure_urf * pprime.values()));
}
} // namespace prism::algo
