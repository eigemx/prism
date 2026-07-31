#include "export.h"

#include <fmt/format.h>

#include <array>
#include <fstream>
#include <stdexcept>
#include <vtu11/vtu11.hpp>

#include "prism/constants.h"
#include "prism/exceptions.h"
#include "prism/log.h"
#include "prism/types.h"

namespace prism::output {

namespace {

/** @brief Caches vertex-to-cell adjacency and cell centres for a mesh.
 *
 *  Populated once per mesh by ensureCache(). Used by cellToPoint() to
 *  accelerate inverse distance weighting (IDW) interpolation.
 *
 *  The cache is keyed on the mesh address, so it assumes the mesh is not
 *  mutated in place (PMesh::cells() has a non-const overload and
 *  CuthillMckee::reorder rewrites cells) and is not freed then reallocated at
 *  the same address. It is also not thread-safe. */
struct MeshCache {
    std::vector<std::vector<size_t>> vert_to_cells;
    std::vector<Vector3d> cell_centres;
    const mesh::PMesh* mesh_key = nullptr;
};

/** @brief Returns the static MeshCache singleton. */
auto getMeshCache() -> MeshCache& {
    static MeshCache cache;
    return cache;
}

/** @brief Ensures the mesh cache is populated for the given mesh.
 *
 *  Rebuilds the cache only when the mesh pointer differs from the cached
 *  one. Populates vertex-to-cell adjacency and cell centre data. */
void ensureCache(const mesh::PMesh& pmesh) {
    auto& cache = getMeshCache();
    const auto* key = &pmesh;

    if (cache.mesh_key == key) {
        return;
    }

    const auto& cells = pmesh.cells();
    const auto n_cells = cells.size();

    cache.vert_to_cells.assign(pmesh.vertices().size(), {});
    for (size_t ci = 0; ci < n_cells; ++ci) {
        for (auto vi : cells[ci].verticesIds()) {
            cache.vert_to_cells[vi].push_back(ci);
        }
    }

    cache.cell_centres.resize(n_cells);
    for (size_t ci = 0; ci < n_cells; ++ci) {
        cache.cell_centres[ci] = cells[ci].center();
    }

    cache.mesh_key = key;
}

/** @brief Interpolates cell-centred data to vertices using inverse distance weighting.
 *  @param cell_data Flat cell-centred data (n_cells * n_components).
 *  @param pmesh The mesh.
 *  @param n_components Number of components per value (1 for scalar, 3 for vector, etc.).
 *  @return Vertex data (n_vertices * n_components). */
auto cellToPoint(const std::vector<f64>& cell_data, const mesh::PMesh& pmesh, size_t n_components)
    -> std::vector<f64> {
    ensureCache(pmesh);
    const auto& cache = getMeshCache();
    const auto& verts = pmesh.vertices();
    const auto n_verts = verts.size();
    auto result = std::vector<f64>(n_verts * n_components, 0.0);

    for (size_t vertex_id = 0; vertex_id < n_verts; ++vertex_id) {
        const auto& vertex_adj_cells = cache.vert_to_cells[vertex_id];
        if (vertex_adj_cells.empty()) {
            throw error::InvalidMesh(fmt::format(
                "Vertex {} has no adjacent cells; point data cannot be computed.", vertex_id));
        }

        auto weights = std::vector<f64>(vertex_adj_cells.size());
        f64 sum_weights = 0.0;
        for (size_t i = 0; i < vertex_adj_cells.size(); ++i) {
            f64 distance_norm = (verts[vertex_id] - cache.cell_centres[vertex_adj_cells[i]]).norm();
            weights[i] = 1.0 / std::max(distance_norm, EPSILON);
            sum_weights += weights[i];
        }
        for (size_t i = 0; i < vertex_adj_cells.size(); ++i) {
            f64 weight_i = weights[i] / sum_weights;
            auto cell_id = vertex_adj_cells[i];
            for (size_t component_i = 0; component_i < n_components; ++component_i) {
                auto vertex_pos = (vertex_id * n_components) + component_i;
                auto cell_pos = (cell_id * n_components) + component_i;
                result[vertex_pos] += weight_i * cell_data[cell_pos];
            }
        }
    }
    return result;
}

/** @brief Converts a WriteMode enum to the VTU11 format string.
 *  @param mode The write mode.
 *  @return "Ascii" or "RawBinaryCompressed". */
auto toVtuMode(WriteMode mode) -> std::string {
    switch (mode) {
        case WriteMode::Ascii: return "Ascii";
        case WriteMode::RawBinaryCompressed: return "RawBinaryCompressed";
        default:
            throw std::logic_error(
                fmt::format("prism::output::toVtuMode(mode): unsupported WriteMode value {}",
                            static_cast<int>(mode)));
    }
}

/** @brief Flattens a scalar field's cell-centred values into a vector.
 *  @param field The scalar field.
 *  @param pmesh The mesh.
 *  @return Vector of cell-centred values (size = n_cells). */
auto flattenScalar(const field::IScalar& field, const mesh::PMesh& pmesh) -> std::vector<f64> {
    auto data = std::vector<f64>();
    const auto n_cells = pmesh.cells().size();
    data.reserve(n_cells);
    for (std::size_t i = 0; i < n_cells; ++i) {
        data.push_back(field.valueAtCell(i));
    }
    return data;
}

/** @brief Flattens a vector field's cell-centred values into an interleaved vector.
 *  @param field The vector field.
 *  @param pmesh The mesh.
 *  @return Interleaved vector data (x, y, z repeating; size = n_cells * 3). */
auto flattenVector(const field::IVector& field, const mesh::PMesh& pmesh) -> std::vector<f64> {
    auto data = std::vector<f64>();
    const auto n_cells = pmesh.cells().size();
    data.reserve(n_cells * 3);
    for (std::size_t i = 0; i < n_cells; ++i) {
        auto v = field.valueAtCell(i);
        data.push_back(v.x());
        data.push_back(v.y());
        data.push_back(v.z());
    }
    return data;
}

/** @brief Flattens a tensor field's cell-centred values into a row-major vector.
 *  @param field The tensor field.
 *  @param pmesh The mesh.
 *  @return Row-major tensor data (size = n_cells * 9). */
auto flattenTensor(const field::ITensor& field, const mesh::PMesh& pmesh) -> std::vector<f64> {
    const auto n_cells = pmesh.cells().size();
    auto data = std::vector<f64>();
    data.reserve(n_cells * 9);
    for (std::size_t i = 0; i < n_cells; ++i) {
        Tensor3d tensor = field.valueAtCell(i);
        data.push_back(tensor(0, 0));
        data.push_back(tensor(0, 1));
        data.push_back(tensor(0, 2));
        data.push_back(tensor(1, 0));
        data.push_back(tensor(1, 1));
        data.push_back(tensor(1, 2));
        data.push_back(tensor(2, 0));
        data.push_back(tensor(2, 1));
        data.push_back(tensor(2, 2));
    }
    return data;
}

/** @brief Low-level VTU writer. Builds the VTU mesh from a PMesh and delegates to vtu11::writeVtu.
 *  @param pmesh The mesh (vertices + cells).
 *  @param info Per-dataset metadata (name, type, components).
 *  @param data Per-dataset data (CellData then PointData, in same order as info).
 *  @param file_name Output .vtu file path.
 *  @param mode Write mode (Ascii or RawBinaryCompressed). */
void writeVtu(const mesh::PMesh& pmesh,
              const std::vector<vtu11::DataSetInfo>& info,
              const std::vector<vtu11::DataSetData>& data,
              const std::string& file_name,
              WriteMode mode) {
    auto points = std::vector<f64>();
    auto connectivity = std::vector<vtu11::VtkIndexType>();
    auto offsets = std::vector<vtu11::VtkIndexType>();
    auto types = std::vector<vtu11::VtkCellType>();

    points.reserve(pmesh.vertices().size() * 3);
    for (const auto& vertex : pmesh.vertices()) {
        points.push_back(vertex.x());
        points.push_back(vertex.y());
        points.push_back(vertex.z());
    }

    offsets.reserve(pmesh.cells().size());
    types.reserve(pmesh.cells().size());
    for (const auto& cell : pmesh.cells()) {
        connectivity.insert(
            connectivity.end(), cell.verticesIds().begin(), cell.verticesIds().end());
        offsets.push_back(static_cast<vtu11::VtkIndexType>(connectivity.size()));

        switch (cell.verticesIds().size()) {
            /// TODO: This supports only hexahedra, tetrahedra, and prisms. General polygons are not
            /// supported by vtu11.
            case 8: types.push_back(12); break; // hexahedron
            case 4: types.push_back(10); break; // tetrahedron
            case 6: types.push_back(13); break; // prism

            default: throw std::runtime_error("Unsupported cell type for vtu export.");
        }
    }

    auto mesh = vtu11::Vtu11UnstructuredMesh {
        .points_ = points, .connectivity_ = connectivity, .offsets_ = offsets, .types_ = types};
    vtu11::writeVtu(file_name, mesh, info, data, toVtuMode(mode));
}

/** @brief Flattens a field, interpolates it to points, and appends CellData and PointData.
 *  @param f Field pointer.
 *  @param n_components Number of components per value (1 scalar, 3 vector, 9 tensor).
 *  @param flatten Flatten function taking the field and its mesh.
 *  @param expected_mesh Mesh all exported fields must share.
 *  @param info Appended with CellData and PointData metadata.
 *  @param data Appended with the flattened cell and interpolated point values. */
template <typename FieldPtr, typename FlattenFn>
void appendFieldData(FieldPtr f,
                     size_t n_components,
                     FlattenFn flatten,
                     const mesh::PMesh* expected_mesh,
                     std::vector<vtu11::DataSetInfo>& info,
                     std::vector<vtu11::DataSetData>& data) {
    if (f->mesh().get() != expected_mesh) {
        throw error::InvalidMesh(fmt::format(
            "Field '{}' does not share the mesh of the other exported fields.", f->name()));
    }
    auto cell_data = flatten(*f, *expected_mesh);
    auto point_data = cellToPoint(cell_data, *expected_mesh, n_components);

    info.emplace_back(f->name(), vtu11::DataSetType::CellData, n_components);
    data.push_back(std::move(cell_data));
    info.emplace_back(f->name(), vtu11::DataSetType::PointData, n_components);
    data.push_back(std::move(point_data));
}

} // namespace

void toVtu(const field::IScalar& field, const std::string& file_name, WriteMode mode) {
    std::array<const field::IScalar*, 1> scalars {&field};
    toVtu(scalars, {}, {}, file_name, mode);
}

void toVtu(const field::IVector& field, const std::string& file_name, WriteMode mode) {
    std::array<const field::IVector*, 1> vectors {&field};
    toVtu({}, vectors, {}, file_name, mode);
}

void toVtu(const field::ITensor& field, const std::string& file_name, WriteMode mode) {
    std::array<const field::ITensor*, 1> tensors {&field};
    toVtu({}, {}, tensors, file_name, mode);
}

void toVtu(std::span<const field::IScalar*> scalars,
           std::span<const field::IVector*> vectors,
           std::span<const field::ITensor*> tensors,
           const std::string& file_name,
           WriteMode mode) {
    const mesh::PMesh* pmesh_ptr = nullptr;
    if (!scalars.empty()) {
        pmesh_ptr = scalars.front()->mesh().get();
    } else if (!vectors.empty()) {
        pmesh_ptr = vectors.front()->mesh().get();
    } else if (!tensors.empty()) {
        pmesh_ptr = tensors.front()->mesh().get();
    } else {
        log::debug("prism::output::toVtu(): no fields provided; skipping export.");
        return;
    }

    const auto& pmesh = *pmesh_ptr;

    const size_t n_scalars = scalars.size();
    const size_t n_vectors = vectors.size();
    const size_t n_tensors = tensors.size();
    auto info = std::vector<vtu11::DataSetInfo>();
    auto data = std::vector<vtu11::DataSetData>();
    info.reserve((n_scalars * 2) + (n_vectors * 2) + (n_tensors * 2));
    data.reserve((n_scalars * 2) + (n_vectors * 2) + (n_tensors * 2));

    for (const auto* f : scalars) {
        appendFieldData(f, 1, flattenScalar, pmesh_ptr, info, data);
    }
    for (const auto* f : vectors) {
        appendFieldData(f, 3, flattenVector, pmesh_ptr, info, data);
    }
    for (const auto* f : tensors) {
        appendFieldData(f, 9, flattenTensor, pmesh_ptr, info, data);
    }

    writeVtu(pmesh, info, data, file_name, mode);
}

void FieldRegistry::add(SharedPtr<field::IScalar> f) {
    _scalars.push_back(std::move(f));
}

void FieldRegistry::add(SharedPtr<field::IVector> f) {
    _vectors.push_back(std::move(f));
}

void FieldRegistry::add(SharedPtr<field::ITensor> f) {
    _tensors.push_back(std::move(f));
}

void FieldRegistry::exportAll(const std::string& file_name, WriteMode mode) const {
    auto scalar_ptrs = std::vector<const field::IScalar*>();
    scalar_ptrs.reserve(_scalars.size());
    for (const auto& f : _scalars) {
        scalar_ptrs.push_back(f.get());
    }

    auto vector_ptrs = std::vector<const field::IVector*>();
    vector_ptrs.reserve(_vectors.size());
    for (const auto& f : _vectors) {
        vector_ptrs.push_back(f.get());
    }

    auto tensor_ptrs = std::vector<const field::ITensor*>();
    tensor_ptrs.reserve(_tensors.size());
    for (const auto& f : _tensors) {
        tensor_ptrs.push_back(f.get());
    }

    toVtu(scalar_ptrs, vector_ptrs, tensor_ptrs, file_name, mode);
}

TimeWriter::TimeWriter(String prefix, size_t write_interval, WriteMode mode, bool write_pvd)
    : _prefix(std::move(prefix)), _interval(write_interval), _mode(mode), _write_pvd(write_pvd) {
    if (write_interval == 0) {
        throw std::invalid_argument(fmt::format(
            "prism::output::TimeWriter(prefix, write_interval, ...): write_interval must be "
            "positive, got {}",
            write_interval));
    }
}

void TimeWriter::add(SharedPtr<field::IScalar> f) {
    _registry.add(std::move(f));
}

void TimeWriter::add(SharedPtr<field::IVector> f) {
    _registry.add(std::move(f));
}

void TimeWriter::add(SharedPtr<field::ITensor> f) {
    _registry.add(std::move(f));
}

void TimeWriter::advance(f64 dt) {
    ++_step;
    _time += dt;
}

void TimeWriter::write() {
    if (_step % _interval != 0) {
        return;
    }
    if (_registry.scalars().empty() && _registry.vectors().empty() && _registry.tensors().empty()) {
        return;
    }

    auto file_name = fmt::format("{}_{:04d}.vtu", _prefix, _step);
    _registry.exportAll(file_name, _mode);

    if (_write_pvd) {
        _pvd_entries.emplace_back(_time, file_name);
        flushPvd();
    }
}

void TimeWriter::writeAndAdvance(f64 dt) {
    advance(dt);
    write();
}

void TimeWriter::flushPvd() const {
    auto pvd_name = fmt::format("{}.pvd", _prefix);
    auto f = std::ofstream(pvd_name);

    f.exceptions(std::ios::failbit | std::ios::badbit);
    f << "<?xml version=\"1.0\"?>\n";
    f << "<VTKFile type=\"Collection\" version=\"0.1\">\n";
    f << "  <Collection>\n";

    for (const auto& [t, name] : _pvd_entries) {
        f << fmt::format("    <DataSet timestep=\"{:.6g}\" file=\"{}\"/>\n", t, name);
    }

    f << "  </Collection>\n";
    f << "</VTKFile>\n";
}

} // namespace prism::output
