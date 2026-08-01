#include "SIMPLE.h"

#include <fmt/format.h>

#include <cstddef>
#include <stdexcept>

#include "prism/equation/transport.h"
#include "prism/log.h"
#include "prism/numerics/bicgstab.h"
#include "prism/operations/rhie_chow.h"
#include "prism/scheme/diffusion/diffusion.h"
#include "prism/scheme/source/divergence.h"
#include "prism/types.h"

namespace prism::algo {
using field::Pressure;

auto solveImplicitMomentum(SIMPLEControls controls, std::span<eqn::Momentum*> momentum_predictors)
    -> std::vector<report::Entry> {
    auto reports = std::vector<report::Entry>();
    log::debug("prism::algo::solveMomentumImplicitly(): solving momentum equations");

    auto momentum_solver = solver::BiCGSTAB();
    for (auto* eqn : momentum_predictors) {
        eqn->setUnderRelaxFactor(controls.momentum_urf);

        // Keep the relaxed matrix after solving so pressureEquationCoeffsTensor() can reuse its
        // diagonal to build D without re-assembling the equation. So, we run solve() with
        // keep_matrix = true, and we don't need to provide a reference cell or value for momentum
        // equations.
        auto result = momentum_solver.solve(
            *eqn, controls.momentum_max_iter, controls.momentum_residual, NullOption, 0.0, false);
        reports.push_back({.field_name = eqn->field()->name(),
                           .n_iterations = result.iteration(),
                           .initial_residual = result.initialResidual(),
                           .final_residual = result.finalResidual(),
                           .converged = result.hasConverged()});
    }
    return reports;
}

void constrainPPrime(SharedPtr<field::Pressure>& pprime) {
    // we need to reset the pprime field to zero at the boundaries where a Dirichlet or a symmetry
    // condition is applied.
    VectorXd face_values = VectorXd::Zero(pprime->mesh()->faceCount());

    for (const auto& patch : pprime->mesh()->boundaryPatches()) {
        if (patch.isEmpty()) {
            continue; // skip empty patches
        }

        const auto& bc = patch.getBoundaryCondition("P");
        const auto& handler = pprime->boundaryHandlersManager().getHandler(bc.kindString());

        if (handler == nullptr) {
            throw error::NonImplementedBoundaryCondition(
                "prism::algo::constrainPPrime()", patch.name(), bc.kindString());
        }

        if (handler->isDirichlet() || handler->name() == "symmetry") {
            for (const auto& face_id : patch.facesIds()) {
                face_values[face_id] = 0.0;
            }
            continue;
        }

        // keep the values at the faces for other boundary conditions
        for (const auto& face_id : patch.facesIds()) {
            const auto& face = pprime->mesh()->face(face_id);
            face_values[face_id] = pprime->valueAtFace(face);
        }
    }
    pprime->setFaceValues(face_values);
}

auto pressureEquationCoeffsTensor(std::span<eqn::Momentum*> momentum_predictors,
                                  const SharedPtr<field::Pressure>& p) -> SharedPtr<field::Tensor> {
    // The momentum equations were already assembled and relaxed by the momentum phase
    // (solveImplicitMomentum / solveExplicitMomentum, which keep the matrix for reuse here);
    // we only need their (relaxed) diagonal to build D = V / aP.
    const auto& mesh = p->mesh();
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
        D_data[i].diagonal() << Du[i], Dv[i], Dw[i];
    }

    for (auto* eqn : momentum_predictors) {
        eqn->zeroOutCoeffs();
    }

    return makeShared<field::Tensor>("D", mesh, std::move(D_data));
}

auto needsPressureReference(const SharedPtr<field::Pressure>& pprime) -> bool {
    for (const auto& patch : pprime->mesh()->boundaryPatches()) {
        if (patch.isEmpty()) {
            continue;
        }
        const auto& bc = patch.getBoundaryCondition(pprime->name());
        const auto& handler = pprime->boundaryHandlersManager().getHandler(bc.kindString());
        if (handler != nullptr && handler->isDirichlet()) {
            return false;
        }
    }
    return true;
}

auto solvePressureEquation(SIMPLEControls controls,
                           std::span<eqn::Momentum*> momentum_predictors,
                           SharedPtr<field::Velocity>& U,
                           SharedPtr<field::Velocity>& mdot,
                           SharedPtr<field::Pressure>& p) -> PressureEquationResult {
    using namespace scheme::diffusion;
    SharedPtr<field::Tensor> D = pressureEquationCoeffsTensor(momentum_predictors, p);

    // Rhie-Chow interpolation for velocity face values
    log::debug(
        "prism::algo::solvePressureEquation(): applying Rhie-Chow correction to faces velocity");

    // first, we update convective flux with latest values of velocity field before applying the
    // Rhie-Chow correction.
    mdot->mapFaceValues([&](const mesh::Face& face) { return U->valueAtFace(face); });
    mdot->mapInteriorFaceValues(
        [&](const mesh::Face& face) { return ops::rhieChowCorrectFace(face, *U, *D, *p); });

    // pressure correction field created with same name as pressure field to get same boundary
    // conditions without having to define _pprime in fields.json file.
    auto pprime = makeShared<field::Pressure>(Pressure(p->name(), p->mesh(), 0.0));

    // The corrector should reset to zero the correction field at every iteration and should
    // also apply a zero value at all boundaries for which a Dirichlet (fixed) boundary
    // condition is used for the pressure.
    constrainPPrime(pprime);

    // If no reference cell was explicitly set and all pressure BCs are Neumann (closed
    // domain), auto-select cell 0 to break the matrix singularity (OpenFOAM compatibility).
    if (!controls.p_ref_cell.has_value() && needsPressureReference(pprime)) {
        controls.p_ref_cell = 0;
        log::debug(
            "prism::algo::solvePressureEquation(): all-Neumann pressure BCs "
            "detected, auto-selecting p_ref_cell = 0");
    }

    /// TODO: based on number of non-orthogonal corrections in _controls, we should check if we need
    /// diffusion::Corrected or diffusion::NonCorrected
    using laplacian_p = Corrected;
    using div_U = scheme::source::Divergence<Sign::Negative>;

    auto pEqn = eqn::Transport(laplacian_p(D, pprime), // - ∇.(D ∇P')
                               div_U(mdot)             // == - (∇.U)
    );

    auto reports = std::vector<report::Entry>();
    log::debug("prism::algo::solvePressureEquation(): solving pressure equation");
    auto p_solver = solver::BiCGSTAB();

    for (std::size_t i = 0; i <= controls.non_ortho_correctors; ++i) {
        auto result = p_solver.solve(pEqn,
                                     controls.pressure_max_iter,
                                     controls.pressure_residual,
                                     controls.p_ref_cell,
                                     controls.p_ref_value);
        reports.push_back({fmt::format("{}_corr_{}", pprime->name(), i),
                           result.iteration(),
                           result.initialResidual(),
                           result.finalResidual(),
                           result.hasConverged()});
    }

    return {.pprime = pprime, .D = D, .reports = std::move(reports)};
}

IncompressibleSIMPLE::IncompressibleSIMPLE(SIMPLEControls controls) : _controls(controls) {}

auto IncompressibleSIMPLE::step(std::span<eqn::Momentum*> momentum_predictors,
                                SharedPtr<field::Velocity>& U,
                                SharedPtr<field::Velocity>& mdot,
                                SharedPtr<field::Pressure>& p) -> std::vector<report::Entry> {
    if (momentum_predictors.size() != 2 && momentum_predictors.size() != 3) {
        throw std::runtime_error(
            fmt::format("prism::algo::IncompressibleSIMPLE::step() expects 2 or 3 momentum "
                        "predictors, not {}",
                        momentum_predictors.size()));
    }

    auto reports = solveImplicitMomentum(_controls, momentum_predictors);
    auto result = solvePressureEquation(_controls, momentum_predictors, U, mdot, p);
    reports.insert(reports.end(), result.reports.begin(), result.reports.end());
    correctFields(U, mdot, p, result.D, result.pprime, _controls.pressure_urf);
    return reports;
}

void correctFields(SharedPtr<field::Velocity>& U,
                   SharedPtr<field::Velocity>& mdot,
                   SharedPtr<field::Pressure>& p,
                   const SharedPtr<field::Tensor>& D,
                   SharedPtr<field::Pressure>& pprime,
                   double pressure_urf) {
    // update velocity field
    U->mapCellValues([&](const mesh::Cell& cell) -> Vector3d {
        // Eigen's lazy expression templates hold references to temporaries.
        // Without `-> Vector3d`, the lambda would return an expression
        // template with dangling references after temporaries are destroyed,
        // causing stack-use-after-return errors under AddressSanitizer.
        return U->valueAtCell(cell) - (D->valueAtCellRef(cell) * pprime->gradAtCell(cell));
    });

    // update mass flow rate at interior faces
    mdot->mapInteriorFaceValues([&](const mesh::Face& face) -> Vector3d {
        // Same as above, we need to return a Vector3d to avoid dangling references.
        return mdot->valueAtFace(face) - (D->valueAtFace(face) * pprime->gradAtFace(face));
    });

    // update pressure
    p->update(p->values() + (pressure_urf * pprime->values()));
}
} // namespace prism::algo
