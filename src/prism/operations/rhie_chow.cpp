#include "rhie_chow.h"

#include <cstddef>

#include "prism/field.h"
#include "prism/mesh/utilities.h"

namespace prism::ops {

void rhie_chow_correct(field::Vector& U,
                       const field::Tensor& D,
                       gradient::AbstractGradient* p_grad_scheme) {
    const auto& mesh = U.mesh();

    VectorXd u_face_data;
    VectorXd v_face_data;
    VectorXd w_face_data;
    u_face_data.resize(mesh.n_faces());
    v_face_data.resize(mesh.n_faces());
    w_face_data.resize(mesh.n_faces());

    for (const auto& face : mesh.interior_faces()) {
        const std::size_t face_id = face.id();
        const Vector3d Uf = U.value_at_face(face);
        const Matrix3d Df = D.value_at_face(face);
        const Vector3d grad_p_f = p_grad_scheme->gradient_at_face(face);

        // TODO: owner and neighbor cells finders are all over the codebase
        // this should be a utility function
        const auto& owner = mesh.cell(face.owner());
        const auto& neighbor = mesh.cell(face.neighbor().value());
        const double gc = mesh::geo_weight(owner, neighbor, face);
        Vector3d avg_grad = gc * p_grad_scheme->gradient_at_cell(owner);
        avg_grad += (1 - gc) * p_grad_scheme->gradient_at_cell(neighbor);

        // Equation 15.60
        Vector3d Uf_corrected = Uf - (Df * (grad_p_f - avg_grad));
        u_face_data[face_id] = Uf_corrected.x();
        v_face_data[face_id] = Uf_corrected.y();
        w_face_data[face_id] = Uf_corrected.z();
    }

    for (const auto& face : mesh.boundary_faces()) {
        const std::size_t face_id = face.id();
        const auto& owner = mesh.cell(face.owner());

        // Equation (15.110)
        const Vector3d Ub = U.value_at_cell(owner);
        const Matrix3d& Df = D.value_at_cell(owner);
        const Vector3d grad_p_f = p_grad_scheme->gradient_at_face(face);
        const Vector3d grad_p_C = p_grad_scheme->gradient_at_cell(owner);

        // Equation (15.111)
        Vector3d Ub_corrected = Ub - (Df * (grad_p_f - grad_p_C));
        u_face_data[face_id] = Ub_corrected.x();
        v_face_data[face_id] = Ub_corrected.y();
        w_face_data[face_id] = Ub_corrected.z();
    }

    U.x().set_face_values(u_face_data);
    U.y().set_face_values(v_face_data);
    U.z().set_face_values(w_face_data);
}
} // namespace prism::ops