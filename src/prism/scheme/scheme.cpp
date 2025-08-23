#include "scheme.h"

#include "prism/field/scalar.h"

namespace prism::scheme {

IPartialScheme::IPartialScheme(std::size_t n_cells) : RHSProvider(n_cells) {}

IFullScheme::IFullScheme(const SharedPtr<field::Scalar>& field)
    : LinearSystem(field->mesh()->cellCount()), _field(field) {}

void IFullScheme::apply() {
    applyBoundary();

    const auto& interior_faces = _field->mesh()->interiorFaces();
    std::for_each(interior_faces.begin(), interior_faces.end(), [this](const mesh::Face& face) {
        applyInterior(face);
    });

    // we've inserted all the triplets, now we can collect them into the matrix
    this->collect();
}

auto IFullScheme::field() -> SharedPtr<field::Scalar>& {
    return _field;
}


auto IFullScheme::field() const -> const SharedPtr<field::Scalar>& {
    return _field;
}

} // namespace prism::scheme
