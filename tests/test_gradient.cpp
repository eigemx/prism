#include <prism/prism.h>

#include <catch2/catch_test_macros.hpp>
#include <cmath>

#include "test_utils.h"

using namespace prism;
using namespace prism::test;

namespace {

class TestGradient : public gradient::GreenGauss {
  public:
    using GreenGauss::GreenGauss;
};

/** @brief RMS error of the cell gradients against a constant exact gradient. */
auto gradientRmsError(const field::Scalar& phi, const mesh::PMesh& mesh, const Vector3d& exact)
    -> f64 {
    f64 sum = 0.0;
    size_t count = 0;

    for (const auto& cell : mesh.cells()) {
        sum += (phi.gradAtCell(cell) - exact).squaredNorm();
        ++count;
    }

    REQUIRE(count > 0);
    return std::sqrt(sum / static_cast<f64>(count));
}

} // namespace

TEST_CASE("custom gradient schemes can be registered", "[gradient]") {
    auto mesh = loadMesh("tests/cases/cavity/mesh.unv");
    gradient::registerGradScheme<TestGradient>({"test-gradient", "testGradient"});

    auto scheme = gradient::createGradScheme("testGradient", mesh);

    REQUIRE(scheme != nullptr);
    REQUIRE(dynamic_cast<TestGradient*>(scheme.get()) != nullptr);
    REQUIRE(gradient::createGradScheme("does-not-exist", mesh) == nullptr);
}

TEST_CASE("skew-corrected Green-Gauss is more accurate than plain on a skewed mesh", "[gradient]") {
    auto mesh = loadMesh("tests/cases/cavity/mesh.unv");
    const Vector3d exact_grad {1.3, -0.7, 0.0};

    auto phi = makeScalarField(
        mesh, "P", [&exact_grad](const auto& cell) { return exact_grad.dot(cell.center()); });

    phi.setGradScheme(makeShared<gradient::GreenGauss>(mesh));
    const f64 plain_error = gradientRmsError(phi, *mesh, exact_grad);

    phi.setGradScheme(makeShared<gradient::SkewCorrectedGreenGauss>(mesh));
    const f64 corrected_error = gradientRmsError(phi, *mesh, exact_grad);

    REQUIRE(plain_error > 0.0);
    REQUIRE(corrected_error < plain_error);
}
