#include "gradient.h"

#include <map>

#include "green_gauss.h"
#include "least_squares.h"

namespace prism::gradient {
namespace {
class GradSchemeRegistry {
  public:
    static auto instance() -> GradSchemeRegistry&;

    void registerScheme(std::initializer_list<std::string> names, GradSchemeCreator creator);

    template <typename Scheme>
    void registerScheme(std::initializer_list<std::string> names) {
        registerScheme(names, [](const SharedPtr<mesh::PMesh>& mesh) -> SharedPtr<IGradient> {
            return makeShared<Scheme>(mesh);
        });
    }

    auto createScheme(const std::string& name, const SharedPtr<mesh::PMesh>& mesh) const
        -> SharedPtr<IGradient>;

  private:
    GradSchemeRegistry();

    std::map<std::string, GradSchemeCreator> _creators;
};

auto GradSchemeRegistry::instance() -> GradSchemeRegistry& {
    static GradSchemeRegistry registry;
    return registry;
}

GradSchemeRegistry::GradSchemeRegistry() {
    // add default gradient schemes implemented by prism
    registerScheme<LeastSquares>({"least-squares", "leastSquares"});
    registerScheme<GreenGauss>({"green-gauss", "greenGauss"});
    registerScheme<SkewCorrectedGreenGauss>(
        {"green-gauss-skew", "greenGaussSkew", "skewCorrectedGreenGauss"});
}

void GradSchemeRegistry::registerScheme(std::initializer_list<std::string> names,
                                        GradSchemeCreator creator) {
    for (const auto& name : names) {
        _creators.insert_or_assign(name, creator);
    }
}

auto GradSchemeRegistry::createScheme(const std::string& name,
                                      const SharedPtr<mesh::PMesh>& mesh) const
    -> SharedPtr<IGradient> {
    const auto it = _creators.find(name);
    if (it == _creators.end()) {
        return nullptr;
    }
    return it->second(mesh);
}
} // namespace

auto createGradScheme(const std::string& name, const SharedPtr<mesh::PMesh>& mesh)
    -> SharedPtr<IGradient> {
    return GradSchemeRegistry::instance().createScheme(name, mesh);
}

void registerGradScheme(std::initializer_list<std::string> names, GradSchemeCreator creator) {
    GradSchemeRegistry::instance().registerScheme(names, std::move(creator));
}

IGradient::IGradient() {
    log::debug("prism::gradient::IGradient() adding default boundary handlers for IGradient");
    this->boundaryHandlersManager().template addHandler<boundary::Fixed>();
    this->boundaryHandlersManager().template addHandler<boundary::Outlet>();
    this->boundaryHandlersManager().template addHandler<boundary::Symmetry>();
    this->boundaryHandlersManager().template addHandler<boundary::ZeroGradient>();
    this->boundaryHandlersManager().template addHandler<boundary::NoSlip>();
}

auto IGradient::gradAtFace(const mesh::Face& face, const field::IScalar& field) -> Vector3d {
    // interpolate gradient at surrounding cells to the face center
    if (face.isInterior()) {
        // interior face
        const auto& mesh = field.mesh();
        const auto& owner_cell = mesh->cell(face.owner());
        auto owner_grad = gradAtCell(owner_cell, field);

        const auto& neighbor_cell = mesh->cell(face.neighbor().value());
        auto neighbor_grad = gradAtCell(neighbor_cell, field);

        auto gc = mesh->faceGeometricWeight(face.id());

        // Equation 9.33 without the correction part, a simple linear interpolation.
        return (gc * owner_grad) + ((1. - gc) * neighbor_grad);
    }

    // boundary face
    return gradAtBoundaryFace(face, field);
}

auto IGradient::gradAtBoundaryFace(const mesh::Face& face, const field::IScalar& field)
    -> Vector3d {
    const auto& boundary_patch = field.mesh()->boundaryPatch(face);
    const auto& boundary_condition = boundary_patch.getBoundaryCondition(field.name());

    auto handler = this->boundaryHandlersManager().getHandler(boundary_condition.kindString());

    if (handler == nullptr) {
        throw prism::error::NonImplementedBoundaryCondition(
            "prism::gradient::IGradient::gradAtBoundaryFace()",
            boundary_patch.name(),
            boundary_condition.kindString());
    }

    return handler->get(field, face);
}
} // namespace prism::gradient
