#include <prism/prism.h>

#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <fstream>
#include <string>
#include <string_view>

#include "test_utils.h"

using namespace prism;
using namespace prism::test;

namespace {

auto readFile(const std::filesystem::path& path) -> std::string {
    std::ifstream file(path);
    REQUIRE(file.good());
    return std::string(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>());
}

auto countOccurrences(const std::string& text, std::string_view needle) -> size_t {
    size_t count = 0;
    size_t pos = 0;
    while ((pos = text.find(needle, pos)) != std::string::npos) {
        ++count;
        pos += needle.size();
    }
    return count;
}

auto tempPath(std::string_view name) -> std::filesystem::path {
    return std::filesystem::temp_directory_path() / name;
}

} // namespace

TEST_CASE("toVtu writes CellData and can skip PointData", "[export]") {
    auto mesh = loadMesh("tests/cases/cavity_hex/mesh.unv");
    auto p = makeShared<field::Scalar>("P", mesh, 1.0);

    const auto cell_file = tempPath("prism_export_cell_data.vtu");
    std::filesystem::remove(cell_file);
    output::toVtu(
        *p, cell_file.string(), {.mode = output::WriteMode::Ascii, .write_point_data = false});

    REQUIRE(std::filesystem::exists(cell_file));
    const auto cell_text = readFile(cell_file);
    CHECK(cell_text.find("<CellData>") != std::string::npos);
    CHECK(countOccurrences(cell_text, "Name=\"P\"") == 1);

    const auto point_file = tempPath("prism_export_point_data.vtu");
    std::filesystem::remove(point_file);
    output::toVtu(
        *p, point_file.string(), {.mode = output::WriteMode::Ascii, .write_point_data = true});

    REQUIRE(std::filesystem::exists(point_file));
    const auto point_text = readFile(point_file);
    CHECK(point_text.find("<PointData>") != std::string::npos);
    CHECK(countOccurrences(point_text, "Name=\"P\"") == 2);

    std::filesystem::remove(cell_file);
    std::filesystem::remove(point_file);
}

TEST_CASE("VtuTimeWriter writes the initial state and follows its interval", "[export]") {
    auto mesh = loadMesh("tests/cases/cavity_hex/mesh.unv");
    auto p = makeShared<field::Scalar>("P", mesh, 1.0);

    const auto prefix = tempPath("prism_time_writer");
    const auto step0 = std::filesystem::path(prefix.string() + "_0000.vtu");
    const auto step1 = std::filesystem::path(prefix.string() + "_0001.vtu");
    std::filesystem::remove(step0);
    std::filesystem::remove(step1);

    output::VtuTimeWriter writer({.prefix = prefix.string(),
                                  .interval = 1,
                                  .write_pvd = false,
                                  .mode = output::WriteMode::Ascii,
                                  .write_point_data = false,
                                  .write_initial = true});
    writer.add(p);
    writer.write(); // initial state at step 0
    writer.writeAndAdvance(0.1);
    writer.writeAndAdvance(0.1);

    CHECK(std::filesystem::exists(step0));
    CHECK(std::filesystem::exists(step1));
    CHECK(writer.currentStep() == 2);

    std::filesystem::remove(step0);
    std::filesystem::remove(step1);
}

TEST_CASE("VtuTimeWriter skips the initial state by default", "[export]") {
    auto mesh = loadMesh("tests/cases/cavity_hex/mesh.unv");
    auto p = makeShared<field::Scalar>("P", mesh, 1.0);

    const auto prefix = tempPath("prism_time_writer_default");
    const auto step0 = std::filesystem::path(prefix.string() + "_0000.vtu");
    const auto step1 = std::filesystem::path(prefix.string() + "_0001.vtu");
    std::filesystem::remove(step0);
    std::filesystem::remove(step1);

    output::VtuTimeWriter writer({.prefix = prefix.string(),
                                  .interval = 1,
                                  .write_pvd = false,
                                  .mode = output::WriteMode::Ascii,
                                  .write_point_data = false});
    writer.add(p);
    writer.write(); // no-op: write_initial is false by default
    writer.writeAndAdvance(0.1);

    CHECK_FALSE(std::filesystem::exists(step0));
    CHECK(std::filesystem::exists(step1));
    CHECK(writer.currentStep() == 1);

    std::filesystem::remove(step1);
}
