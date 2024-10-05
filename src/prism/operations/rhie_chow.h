#pragma once

#include "prism/field/ifield.h"
#include "prism/field/pressure.h"
#include "prism/field/tensor.h"
#include "prism/mesh/face.h"
#include "prism/mesh/pmesh.h"
#include "prism/types.h"

namespace prism::ops {

template <field::IVectorBased Vector>
auto rhieChowCorrect(Vector& U, const field::Tensor& D, const field::Pressure& P) -> Vector;

namespace detail {
auto pressureGradCalculated(const mesh::PMesh& mesh,
                            const mesh::Face& face,
                            const field::Pressure& P,
                            const Vector3d& gradp_avg) -> Vector3d;
}

template <field::IVectorBased Vector>
auto rhieChowCorrect(Vector& U, const field::Tensor& D, const field::Pressure& P) -> Vector {
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
        const Vector3d gradp_avg = P.gradAtFace(face);
        const Vector3d gradp = detail::pressureGradCalculated(mesh, face, P, gradp_avg);

        // Equation 15.60
       Vector3d Uf_corrected = Uf - (Df * (gradp - gradp_avg));
        u_face_data[face_id] = Uf_corrected.x();
        v_face_data[face_id] = Uf_corrected.y();
        w_face_data[face_id] = Uf_corrected.z();
    }

    for (const auto& face : mesh.boundaryFaces()) {
        const std::size_t face_id = face.id();
        u_face_data[face_id] = U.x().valueAtFace(face);
        v_face_data[face_id] = U.y().valueAtFace(face);
        w_face_data[face_id] = U.z().valueAtFace(face);
    }

    using Component = typename Vector::ComponentType;
    auto x_rh = Component(U.x().name(), mesh, U.x().values(), Coord::X);
    auto y_rh = Component(U.y().name(), mesh, U.y().values(), Coord::Y);
    auto z_rh = Component(U.z().name(), mesh, U.z().values(), Coord::Z);
    x_rh.setFaceValues(u_face_data);
    y_rh.setFaceValues(v_face_data);
    z_rh.setFaceValues(w_face_data);

    auto components = std::array {x_rh, y_rh, z_rh};
    return Vector {U.name(), mesh, components};
}
} // namespace prism::ops
