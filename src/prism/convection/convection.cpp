#include "convection.h"


namespace prism::convection {
auto CentralDifference::interpolate(
    double m_dot,
    const mesh::Cell& cell,
    const mesh::Cell& neighbor,
    const mesh::Face& face,
    const std::shared_ptr<gradient::GradientSchemeBase>& grad_scheme) -> InterpolationCoeffs {
    // in case `cell` is the upstream cell
    auto face_grad_phi = grad_scheme->gradient_at_face(face);
    auto d_Cf = face.center() - cell.center();
    auto a_C = std::max(m_dot, 0.0);
    auto b = -std::max(m_dot, 0.0) * d_Cf.dot(face_grad_phi);

    // in case 'neighbor' is the upstream cell
    auto d_Nf = face.center() - neighbor.center();
    auto a_D = std::max(-m_dot, 0.0);
    b -= std::max(-m_dot, 0.0) * d_Nf.dot(face_grad_phi);

    return {a_C, a_D, b};
}

auto Upwind::interpolate(double m_dot,
                         const mesh::Cell& cell,
                         const mesh::Cell& neighbor,
                         const mesh::Face& face,
                         const std::shared_ptr<gradient::GradientSchemeBase>& grad_scheme)
    -> InterpolationCoeffs {
    // in case `cell` is the upstream cell
    auto a_C = std::max(m_dot, 0.0);
    // in case 'neighbor' is the upstream cell
    auto a_D = -std::max(-m_dot, 0.0);

    auto b = 0.0;

    return {a_C, a_D, b};
}

auto SecondOrderUpwind::interpolate(
    double m_dot,
    const mesh::Cell& cell,
    const mesh::Cell& neighbor,
    const mesh::Face& face,
    const std::shared_ptr<gradient::GradientSchemeBase>& grad_scheme) -> InterpolationCoeffs {
    // in case `cell` is the upstream cell
    const auto& face_grad_phi = grad_scheme->gradient_at_face(face);
    const auto& cell_grad_phi = grad_scheme->gradient_at_cell(cell);
    const auto& neighbor_grad_phi = grad_scheme->gradient_at_cell(neighbor);

    auto d_Cf = face.center() - cell.center();
    auto correction = d_Cf.dot(2 * cell_grad_phi - face_grad_phi);
    auto a_C = std::max(m_dot, 0.0);
    auto b = -std::max(m_dot, 0.0) * correction;

    // in case 'neighbor' is the upstream cell
    auto d_Nf = face.center() - neighbor.center();
    correction = d_Nf.dot(2 * neighbor_grad_phi - face_grad_phi);
    auto a_D = std::max(-m_dot, 0.0);
    b -= std::max(-m_dot, 0.0) * correction;

    return {a_C, a_D, b};
}
} // namespace prism::convection
