#pragma once

#include <cstddef>

#include "prism/field/ifield.h"
#include "prism/field/pressure.h"
#include "prism/field/tensor.h"
#include "prism/mesh/face.h"
#include "prism/types.h"

namespace prism::ops {

template <field::IVectorBased Vector>
auto rhieChowCorrectFace(const mesh::Face& face,
                         Vector& U,
                         const field::Tensor& D,
                         const field::Pressure& P) -> Vector3d;

template <field::IVectorBased Vector>
auto rhieChowCorrectBoundaryFace(const mesh::Face& face,
                                 Vector& U,
                                 const field::Tensor& D,
                                 const field::Pressure& P) -> Vector3d;

namespace detail {
auto pressureGradCalculated(const mesh::Face& face,
                            const field::Pressure& P,
                            const Vector3d& gradp_avg) -> Vector3d;
}

template <field::IVectorBased Vector>
auto rhieChowCorrectBoundaryFace(const mesh::Face& face,
                                 Vector& U,
                                 const field::Tensor& D,
                                 const field::Pressure& P) -> Vector3d {
    /// TODO: correcting boundary faces returns wrong solution near boundaries.
    const std::size_t face_id = face.id();
    const std::size_t owner_id = face.owner();
    const Vector3d& Uc = U.valueAtCell(owner_id);
    const Matrix3d& Dc = D.valueAtCell(owner_id);

    // Equation 15.60
    Vector3d Uf_corrected =
        Uc - (Dc * (P.gradAtFace(face) - P.gradAtCell(U.mesh()->cell(owner_id))));
    return Uf_corrected;
}

template <field::IVectorBased Vector>
auto rhieChowCorrectFace(const mesh::Face& face,
                         Vector& U,
                         const field::Tensor& D,
                         const field::Pressure& P) -> Vector3d {
    if (face.isBoundary()) {
        return rhieChowCorrectBoundaryFace(face, U, D, P);
    }
    const Vector3d& Uf = U.valueAtFace(face);
    const Matrix3d& Df = D.valueAtFace(face);
    const Vector3d gradp_avg = P.gradAtFace(face);
    const Vector3d gradp_calculated = detail::pressureGradCalculated(face, P, gradp_avg);

    // Equation 15.60
    Vector3d Uf_corrected = Uf - (Df * (gradp_calculated - gradp_avg));
    return Uf_corrected;
}
} // namespace prism::ops
