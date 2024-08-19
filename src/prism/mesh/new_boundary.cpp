#include <fmt/format.h>
#include <spdlog/spdlog.h>
#include <toml++/toml.h>

#include <filesystem>
#include <optional>
#include <set>
#include <stdexcept>
#include <string_view>

#include "boundary.h"


namespace prism::mesh {
auto readTOML(const std::filesystem::path& path) -> toml::table {
    auto fstream {std::ifstream {path}};

    if (!fstream) {
        throw std::runtime_error(fmt::format(
            "prism::mesh::readBoundaryFile(): Failed to open boundary conditions file `{}`",
            path.string()));
    }

    toml::table doc;

    // parse boundary TOML file
    try {
        doc = toml::parse(fstream);
    } catch (const toml::parse_error& e) {
        // file cannot be parsed
        throw std::runtime_error(
            fmt::format("prism::mesh::readBoundaryFile(): Failed to parse boundary condition "
                        "file: `{}`, complete error: `{}`",
                        path.string(),
                        e.what()));
    }
    return doc;
}

auto readFieldsTable(const toml::table& doc) -> std::set<std::string> {
    // read fields table
    auto fields_table {doc["fields"]["names"]};
    if (!fields_table) {
        throw std::runtime_error(
            fmt::format("prism::mesh::readBoundaryFile(): Couldn't find definition for "
                        "fields in boundary conditions file"));
    }

    if (!fields_table.is_array()) {
        throw std::runtime_error(
            fmt::format("prism::mesh::readBoundaryFile(): `fields.names` must be an array"));
    }

    const auto* names_array = fields_table.as_array();
    std::set<std::string> names;

    for (const auto& name : *names_array) {
        if (!name.is_string()) {
            throw std::runtime_error(
                fmt::format("prism::mesh::readBoundaryFile(): `fields.names` array must contain "
                            "only strings"));
        }
        names.insert(name.value<std::string>().value());
    }
    return names;
}

auto readBoundaryFile(const std::filesystem::path& path) -> std::vector<BoundaryPatch> {
    // read fie
    std::vector<BoundaryPatch> boundary_patches;
    auto file = readTOML(path);
    auto fields_names = readFieldsTable(file);

    // print names
    for (const auto& name : fields_names) {
        spdlog::debug("Field name: {}", name);
    }

    return {};
}
} // namespace prism::mesh