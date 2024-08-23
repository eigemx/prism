#include "nonortho.h"

namespace prism::nonortho {

auto OverRelaxedCorrector::decompose(const Vector3d& Sf, const Vector3d& e) const
    -> NonOrthoPair {
    Vector3d Ef = ((Sf.dot(Sf) / (e.dot(Sf) + EPSILON))) * e;
    return {Ef, Sf - Ef};
}

auto MinimumCorrector::decompose(const Vector3d& Sf, const Vector3d& e) const -> NonOrthoPair {
    Vector3d Ef = (e.dot(Sf)) * e;
    return {Ef, Sf - Ef};
}

auto OrthogonalCorrector::decompose(const Vector3d& Sf, const Vector3d& e) const -> NonOrthoPair {
    Vector3d Ef = Sf.norm() * e;
    return {Ef, Sf - Ef};
}

} // namespace prism::nonortho