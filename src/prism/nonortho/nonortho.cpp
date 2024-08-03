#include "nonortho.h"

namespace prism::nonortho {

auto inline OverRelaxedCorrector::decompose(const Vector3d& Sf, const Vector3d& e) const
    -> NonOrthoPair {
    Vector3d Ef = ((Sf.dot(Sf) / (e.dot(Sf) + EPSILON))) * e;
    return {Ef, Sf - Ef};
}

auto inline MinimumCorrector::decompose(const Vector3d& Sf, const Vector3d& e) const
    -> NonOrthoPair {
    Vector3d Ef = (e.dot(Sf)) * e;
    return {Ef, Sf - Ef};
}

auto inline OrthogonalCorrector::decompose(const Vector3d& Sf, const Vector3d& e) const
    -> NonOrthoPair {
    Vector3d Ef = Sf.norm() * e;
    return {Ef, Sf - Ef};
}

} // namespace prism::nonortho