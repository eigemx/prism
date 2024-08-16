#pragma once

#include "prism/field/pressure.h"
#include "prism/field/tensor.h"
#include "prism/gradient/gradient.h"
#include "prism/mesh/face.h"
#include "prism/mesh/pmesh.h"
#include "prism/types.h"

namespace prism::ops {

template <typename Vector, typename GradScheme = gradient::LeastSquares<field::Pressure>>
void correctRhieChow(Vector& U, const field::Tensor& D, const field::Pressure& P);

namespace detail {
auto correctGrad(const mesh::PMesh& mesh,
                 const mesh::Face& face,
                 const field::Pressure& P,
                 const Vector3d& grad_p_f) -> Vector3d;
}

template <typename Vector, typename GradScheme>
void correctRhieChow(Vector& U, const field::Tensor& D, const field::Pressure& P) {
    const auto& mesh = U.mesh();
    GradScheme p_grad_scheme(P);

    VectorXd u_face_data;
    VectorXd v_face_data;
    VectorXd w_face_data;
    u_face_data.resize(mesh.nFaces());
    v_face_data.resize(mesh.nFaces());
    w_face_data.resize(mesh.nFaces());

    for (const auto& face : mesh.interiorFaces()) {
        const std::size_t face_id = face.id();
        const Vector3d& Uf = U.valueAtFace(face);
        const Matrix3d& Df = D.valueAtFace(face);
        const Vector3d& grad_p_f = p_grad_scheme.gradAtFace(face);

        Vector3d grad_p_f_corr = grad_p_f + detail::correctGrad(mesh, face, P, grad_p_f);

        // Equation 15.60
        Vector3d Uf_corrected = Uf - (Df * (grad_p_f_corr - grad_p_f));
        u_face_data[face_id] = Uf_corrected.x();
        v_face_data[face_id] = Uf_corrected.y();
        w_face_data[face_id] = Uf_corrected.z();
    }

    for (const auto& face : mesh.boundaryFaces()) {
        const std::size_t face_id = face.id();
        const auto& owner = mesh.cell(face.owner());

        // Equation (15.110)
        const Vector3d Ub = U.valueAtCell(owner);
        const Matrix3d& Df = D.valueAtCell(owner);
        const Vector3d grad_p_f = p_grad_scheme.gradAtFace(face);
        const Vector3d grad_p_C = p_grad_scheme.gradAtCell(owner);

        // Equation (15.111)
        Vector3d Ub_corrected = Ub - (Df * (grad_p_f - grad_p_C));
        u_face_data[face_id] = Ub_corrected.x();
        v_face_data[face_id] = Ub_corrected.y();
        w_face_data[face_id] = Ub_corrected.z();
    }

    U.x().setFaceValues(u_face_data);
    U.y().setFaceValues(v_face_data);
    U.z().setFaceValues(w_face_data);
}

} // namespace prism::ops