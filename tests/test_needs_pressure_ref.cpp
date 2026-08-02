#include <prism/algorithm/SIMPLE.h>
#include <prism/prism.h>

#include <catch2/catch_test_macros.hpp>

#include "test_utils.h"

using namespace prism;

TEST_CASE("needsPressureReference cavity_hex (all zero-gradient) returns true", "[algo]") {
    auto mesh = test::loadMesh("tests/cases/cavity_hex/mesh.unv");

    auto p = makeShared<field::Pressure>("P", mesh, 0.0);
    REQUIRE(algo::needsPressureReference(p));
}

TEST_CASE("needsPressureReference coarsePipeHex (fixed outlet) returns false", "[algo]") {
    auto mesh = test::loadMesh("tests/cases/coarsePipeHexSymmetrical/mesh.unv");

    auto p = makeShared<field::Pressure>("P", mesh, 0.0);
    REQUIRE_FALSE(algo::needsPressureReference(p));
}

TEST_CASE("needsPressureReference ductSIMPLE (fixed outlet) returns false", "[algo]") {
    auto mesh = test::loadMesh("tests/cases/ductSIMPLE/mesh.unv");

    auto p = makeShared<field::Pressure>("P", mesh, 0.0);
    REQUIRE_FALSE(algo::needsPressureReference(p));
}
