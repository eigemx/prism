#pragma once

#include "prism/field/pressure.h"
#include "prism/field/tensor.h"
#include "prism/mesh/face.h"
#include "prism/mesh/pmesh.h"
#include "prism/types.h"

namespace prism::ops {

template <typename Vector>
void correctRhieChow(Vector& U, const field::Tensor& D, const field::Pressure& P);

namespace detail {
auto pressureGradCalculated(const mesh::PMesh& mesh,
                            const mesh::Face& face,
                            const field::Pressure& P,
                            const Vector3d& grad_p_f) -> Vector3d;
}

template <typename Vector>
void correctRhieChow(Vector& U, const field::Tensor& D, const field::Pressure& P) {
    const auto& mesh = U.mesh();
    VectorXd u_face_data;
    VectorXd v_face_data;
    VectorXd w_face_data;
    u_face_data.resize(mesh.faceCount());
    v_face_data.resize(mesh.faceCount());
    w_face_data.resize(mesh.faceCount());


    for (const auto& face : mesh.interiorFaces()) {
        const std::size_t face_id = face.id();
        const Vector3d& Uf = U.valueAtFace(face);
        const Matrix3d& Df = D.valueAtFace(face);
        const Vector3d& grad_p_f_bar = P.gradAtFace(face);
        const Vector3d grad_p_f = detail::pressureGradCalculated(mesh, face, P, grad_p_f_bar);

        // Equation 15.60
        Vector3d Uf_corrected = Uf - (Df * ((grad_p_f - grad_p_f_bar)));
        u_face_data[face_id] = Uf_corrected.x();
        v_face_data[face_id] = Uf_corrected.y();
        w_face_data[face_id] = Uf_corrected.z();
    }

    for (const auto& face : mesh.boundaryFaces()) {
        const std::size_t face_id = face.id();
        const auto& owner = mesh.cell(face.owner());

        // Equation (15.110)
        //const Vector3d Uc = U.valueAtCell(owner);
        //const Matrix3d& Df = D.valueAtCell(owner);
        //const Vector3d grad_p_f = P.gradAtFace(face);
        //const Vector3d grad_p_C = P.gradAtCell(owner);

        // Equation (15.111)
        //Vector3d Ub_corrected = Uc - (Df * (grad_p_f - grad_p_C));
        // u_face_data[face_id] = Ub_corrected.x();
        // v_face_data[face_id] = Ub_corrected.y();
        // w_face_data[face_id] = Ub_corrected.z();
        u_face_data[face_id] = U.x().valueAtFace(face);
        v_face_data[face_id] = U.y().valueAtFace(face);
        w_face_data[face_id] = U.z().valueAtFace(face);
    }

    U.x().setFaceValues(u_face_data);
    U.y().setFaceValues(v_face_data);
    U.z().setFaceValues(w_face_data);
}
} // namespace prism::ops
