#pragma once

#include "prism/field/ifield.h"
#include "prism/field/pressure.h"
#include "prism/field/tensor.h"
#include "prism/mesh/face.h"
#include "prism/types.h"

namespace prism::ops {

template <field::IVectorBased Vector>
auto rhieChowCorrect(Vector& U, const field::Tensor& D, const field::Pressure& P) -> Vector;

template <field::IVectorBased Vector>
auto rhieChowCorrectFace(const mesh::Face& face,
                         Vector& U,
                         const field::Tensor& D,
                         const field::Pressure& P) -> Vector3d;

namespace detail {
auto pressureGradCalculated(const mesh::Face& face,
                            const field::Pressure& P,
                            const Vector3d& gradp_avg) -> Vector3d;
}

template <field::IVectorBased Vector>
auto rhieChowCorrectFace(const mesh::Face& face,
                         Vector& U,
                         const field::Tensor& D,
                         const field::Pressure& P) -> Vector3d {
    const std::size_t face_id = face.id();
    const Vector3d& Uf = U.valueAtFace(face);
    const Matrix3d& Df = D.valueAtFace(face);
    const Vector3d gradp_avg = P.gradAtFace(face);
    const Vector3d gradp_calculated = detail::pressureGradCalculated(face, P, gradp_avg);

    // Equation 15.60
    Vector3d Uf_corrected = Uf - (Df * (gradp_calculated - gradp_avg));
    return Uf_corrected;
}

/// TODO: delete this
template <field::IVectorBased Vector>
auto rhieChowCorrect(Vector& U, const field::Tensor& D, const field::Pressure& P) -> Vector {
    const auto& mesh = U.mesh();
    VectorXd u_face_data;
    VectorXd v_face_data;
    VectorXd w_face_data;
    u_face_data.resize(mesh->faceCount());
    v_face_data.resize(mesh->faceCount());
    w_face_data.resize(mesh->faceCount());

    for (const auto& face : mesh->interiorFaces()) {
        const std::size_t face_id = face.id();
        const Vector3d& Uf = U.valueAtFace(face);
        const Matrix3d& Df = D.valueAtFace(face);
        const Vector3d gradp_avg = P.gradAtFace(face);
        const Vector3d gradp_calculated = detail::pressureGradCalculated(face, P, gradp_avg);

        // Equation 15.60
        Vector3d Uf_corrected = Uf - (Df * (gradp_calculated - gradp_avg));
        u_face_data[face_id] = Uf_corrected.x();
        v_face_data[face_id] = Uf_corrected.y();
        w_face_data[face_id] = Uf_corrected.z();
    }


    for (const auto& face : mesh->boundaryFaces()) {
        const std::size_t face_id = face.id();
        const auto& patch = mesh->faceBoundaryPatch(face);

        if (patch.isEmpty()) {
            /// TODO: this is a hack to avoid calling valueAtFace() on empty faces. We need to fix
            // this.
            u_face_data[face_id] = 0.0;
            v_face_data[face_id] = 0.0;
            w_face_data[face_id] = 0.0;
            continue;
        }
        u_face_data[face_id] = U.x().valueAtFace(face);
        v_face_data[face_id] = U.y().valueAtFace(face);
        w_face_data[face_id] = U.z().valueAtFace(face);
    }

    using Component = typename Vector::ComponentType;
    auto x_rh = Component(U.x().name(), mesh, U.x().values(), u_face_data, Coord::X);
    auto y_rh = Component(U.y().name(), mesh, U.y().values(), v_face_data, Coord::Y);
    auto z_rh = Component(U.z().name(), mesh, U.z().values(), w_face_data, Coord::Z);
    auto components = std::array {x_rh, y_rh, z_rh};
    return Vector {U.name(), mesh, components};
}
} // namespace prism::ops
