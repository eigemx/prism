#include "export.h"

#include <algorithm>
#include <array>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <vtu11/vtu11.hpp>

#include "prism/constants.h"
#include "prism/exceptions.h"
#include "prism/log.h"
#include "prism/types.h"

namespace prism::output {

namespace {

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

/** @brief Flattens a scalar field's cell-centred values into a vector. */
auto flattenScalar(const field::IScalar& field, const mesh::PMesh& pmesh) -> std::vector<f64> {
    auto data = std::vector<f64>();
    const auto n_cells = pmesh.cells().size();
    data.reserve(n_cells);
    for (size_t i = 0; i < n_cells; ++i) {
        data.push_back(field.valueAtCell(i));
    }
    return data;
}

/** @brief Flattens a vector field's cell-centred values into an interleaved vector. */
auto flattenVector(const field::IVector& field, const mesh::PMesh& pmesh) -> std::vector<f64> {
    auto data = std::vector<f64>();
    const auto n_cells = pmesh.cells().size();
    data.reserve(n_cells * 3);
    for (size_t i = 0; i < n_cells; ++i) {
        auto v = field.valueAtCell(i);
        data.push_back(v.x());
        data.push_back(v.y());
        data.push_back(v.z());
    }
    return data;
}

/** @brief Flattens a tensor field's cell-centred values into a row-major vector. */
auto flattenTensor(const field::ITensor& field, const mesh::PMesh& pmesh) -> std::vector<f64> {
    const auto n_cells = pmesh.cells().size();
    auto data = std::vector<f64>();
    data.reserve(n_cells * 9);
    for (size_t i = 0; i < n_cells; ++i) {
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

// VTK unstructured cell type ids.
constexpr vtu11::VtkCellType vtk_tetra {10};
constexpr vtu11::VtkCellType vtk_hexahedron {12};
constexpr vtu11::VtkCellType vtk_wedge {13};

struct FieldPointers {
    std::vector<const field::IScalar*> scalars;
    std::vector<const field::IVector*> vectors;
    std::vector<const field::ITensor*> tensors;
};

auto fieldPointers(const FieldRegistry& registry) -> FieldPointers {
    FieldPointers pointers;
    pointers.scalars.reserve(registry.scalars().size());
    pointers.vectors.reserve(registry.vectors().size());
    pointers.tensors.reserve(registry.tensors().size());

    for (const auto& f : registry.scalars()) {
        pointers.scalars.push_back(f.get());
    }
    for (const auto& f : registry.vectors()) {
        pointers.vectors.push_back(f.get());
    }
    for (const auto& f : registry.tensors()) {
        pointers.tensors.push_back(f.get());
    }
    return pointers;
}

} // namespace

namespace detail {

/** @brief Writes VTU files for one mesh, caching the mesh geometry and the vertex-to-cell
 * adjacency used by the cell-to-point interpolation. */
class VtuWriter {
  public:
    explicit VtuWriter(SharedPtr<mesh::PMesh> mesh);

    void write(const std::string& file_name,
               std::span<const field::IScalar*> scalars,
               std::span<const field::IVector*> vectors,
               std::span<const field::ITensor*> tensors,
               const VtuOptions& options) const;

  private:
    template <typename FieldPtr, typename FlattenFn>
    void appendFieldData(FieldPtr field,
                         size_t n_components,
                         FlattenFn flatten,
                         std::vector<vtu11::DataSetInfo>& info,
                         std::vector<vtu11::DataSetData>& data,
                         const VtuOptions& options) const;

    auto cellToPoint(const std::vector<f64>& cell_data, size_t n_components) const
        -> std::vector<f64>;

    SharedPtr<mesh::PMesh> _mesh;
    std::vector<f64> _points;
    std::vector<vtu11::VtkIndexType> _connectivity;
    std::vector<vtu11::VtkIndexType> _offsets;
    std::vector<vtu11::VtkCellType> _types;

    // Flattened vertex-to-cell adjacency: cell ids for vertex i are
    // [_vertex_cell_offsets[i], _vertex_cell_offsets[i + 1]).
    std::vector<size_t> _vertex_cell_offsets;
    std::vector<size_t> _vertex_cells;
};

VtuWriter::VtuWriter(SharedPtr<mesh::PMesh> mesh) : _mesh(std::move(mesh)) {
    const auto& mesh_ref = *_mesh;
    const auto& cells = mesh_ref.cells();
    const auto& vertices = mesh_ref.vertices();

    _points.reserve(vertices.size() * 3);
    for (const auto& vertex : vertices) {
        _points.push_back(vertex.x());
        _points.push_back(vertex.y());
        _points.push_back(vertex.z());
    }

    _offsets.reserve(cells.size());
    _types.reserve(cells.size());
    for (const auto& cell : cells) {
        const auto& vertices_ids = cell.verticesIds();
        _connectivity.insert(_connectivity.end(), vertices_ids.begin(), vertices_ids.end());
        _offsets.push_back(static_cast<vtu11::VtkIndexType>(_connectivity.size()));

        // TODO: Cell type is inferred from the vertex count. This works only while prism
        // accepts hexahedral, tetrahedral and wedge cells. Store an explicit mesh::CellType
        // and map it to VTK types when 2D or polyhedral cell support is added.
        switch (vertices_ids.size()) {
            case 8: _types.push_back(vtk_hexahedron); break;
            case 4: _types.push_back(vtk_tetra); break;
            case 6: _types.push_back(vtk_wedge); break;

            default: throw std::runtime_error("Unsupported cell type for vtu export.");
        }
    }

    _vertex_cell_offsets.assign(vertices.size() + 1, 0);
    for (const auto& cell : cells) {
        for (const auto vertex_id : cell.verticesIds()) {
            ++_vertex_cell_offsets[vertex_id + 1];
        }
    }
    for (size_t i = 0; i < vertices.size(); ++i) {
        _vertex_cell_offsets[i + 1] += _vertex_cell_offsets[i];
    }

    _vertex_cells.resize(_vertex_cell_offsets.back());
    auto cursor = _vertex_cell_offsets;
    for (const auto& cell : cells) {
        for (const auto vertex_id : cell.verticesIds()) {
            _vertex_cells[cursor[vertex_id]++] = cell.id();
        }
    }
}

template <typename FieldPtr, typename FlattenFn>
void VtuWriter::appendFieldData(FieldPtr field,
                                size_t n_components,
                                FlattenFn flatten,
                                std::vector<vtu11::DataSetInfo>& info,
                                std::vector<vtu11::DataSetData>& data,
                                const VtuOptions& options) const {
    if (field->mesh().get() != _mesh.get()) {
        throw error::InvalidMesh(fmt::format(
            "Field '{}' does not share the mesh of the other exported fields.", field->name()));
    }

    auto cell_data = flatten(*field, *_mesh);
    std::vector<f64> point_data;
    if (options.write_point_data) {
        point_data = cellToPoint(cell_data, n_components);
    }

    info.emplace_back(field->name(), vtu11::DataSetType::CellData, n_components);
    data.push_back(std::move(cell_data));

    if (options.write_point_data) {
        info.emplace_back(field->name(), vtu11::DataSetType::PointData, n_components);
        data.push_back(std::move(point_data));
    }
}

auto VtuWriter::cellToPoint(const std::vector<f64>& cell_data, size_t n_components) const
    -> std::vector<f64> {
    const auto& vertices = _mesh->vertices();
    const auto& cells = _mesh->cells();
    const auto n_vertices = vertices.size();
    auto result = std::vector<f64>(n_vertices * n_components, 0.0);

    for (size_t vertex_id = 0; vertex_id < n_vertices; ++vertex_id) {
        const auto begin = _vertex_cell_offsets[vertex_id];
        const auto end = _vertex_cell_offsets[vertex_id + 1];

        if (begin == end) {
            throw error::InvalidMesh(fmt::format(
                "Vertex {} has no adjacent cells; point data cannot be computed.", vertex_id));
        }

        f64 sum_weights = 0.0;
        for (auto k = begin; k < end; ++k) {
            const auto cell_id = _vertex_cells[k];
            const f64 distance = (vertices[vertex_id] - cells[cell_id].center()).norm();
            sum_weights += 1.0 / std::max(distance, EPSILON);
        }

        for (auto k = begin; k < end; ++k) {
            const auto cell_id = _vertex_cells[k];
            const f64 distance = (vertices[vertex_id] - cells[cell_id].center()).norm();
            const f64 weight = (1.0 / std::max(distance, EPSILON)) / sum_weights;

            for (size_t component = 0; component < n_components; ++component) {
                result[(vertex_id * n_components) + component] +=
                    weight * cell_data[(cell_id * n_components) + component];
            }
        }
    }
    return result;
}

void VtuWriter::write(const std::string& file_name,
                      std::span<const field::IScalar*> scalars,
                      std::span<const field::IVector*> vectors,
                      std::span<const field::ITensor*> tensors,
                      const VtuOptions& options) const {
    auto info = std::vector<vtu11::DataSetInfo>();
    auto data = std::vector<vtu11::DataSetData>();
    const size_t sets_per_field = options.write_point_data ? 2 : 1;
    info.reserve((scalars.size() + vectors.size() + tensors.size()) * sets_per_field);
    data.reserve((scalars.size() + vectors.size() + tensors.size()) * sets_per_field);

    for (const auto* field : scalars) {
        appendFieldData(field, 1, flattenScalar, info, data, options);
    }
    for (const auto* field : vectors) {
        appendFieldData(field, 3, flattenVector, info, data, options);
    }
    for (const auto* field : tensors) {
        appendFieldData(field, 9, flattenTensor, info, data, options);
    }

    auto mesh = vtu11::Vtu11UnstructuredMesh {
        .points_ = _points,
        .connectivity_ = _connectivity,
        .offsets_ = _offsets,
        .types_ = _types,
    };
    vtu11::writeVtu(file_name, mesh, info, data, toVtuMode(options.mode));
}

} // namespace detail

void toVtu(const field::IScalar& field, const std::string& file_name, const VtuOptions& options) {
    std::array<const field::IScalar*, 1> scalars {&field};
    toVtu(scalars, {}, {}, file_name, options);
}

void toVtu(const field::IVector& field, const std::string& file_name, const VtuOptions& options) {
    std::array<const field::IVector*, 1> vectors {&field};
    toVtu({}, vectors, {}, file_name, options);
}

void toVtu(const field::ITensor& field, const std::string& file_name, const VtuOptions& options) {
    std::array<const field::ITensor*, 1> tensors {&field};
    toVtu({}, {}, tensors, file_name, options);
}

void toVtu(std::span<const field::IScalar*> scalars,
           std::span<const field::IVector*> vectors,
           std::span<const field::ITensor*> tensors,
           const std::string& file_name,
           const VtuOptions& options) {
    SharedPtr<mesh::PMesh> pmesh;
    if (!scalars.empty()) {
        pmesh = scalars.front()->mesh();
    } else if (!vectors.empty()) {
        pmesh = vectors.front()->mesh();
    } else if (!tensors.empty()) {
        pmesh = tensors.front()->mesh();
    } else {
        log::debug("prism::output::toVtu(): no fields provided; skipping export.");
        return;
    }

    detail::VtuWriter writer(std::move(pmesh));
    writer.write(file_name, scalars, vectors, tensors, options);
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

void FieldRegistry::exportAll(const std::string& file_name, const VtuOptions& options) const {
    auto pointers = fieldPointers(*this);
    toVtu(pointers.scalars, pointers.vectors, pointers.tensors, file_name, options);
}

auto FieldRegistry::empty() const noexcept -> bool {
    return _scalars.empty() && _vectors.empty() && _tensors.empty();
}

auto FieldRegistry::mesh() const -> SharedPtr<mesh::PMesh> {
    if (!_scalars.empty()) {
        return _scalars.front()->mesh();
    }
    if (!_vectors.empty()) {
        return _vectors.front()->mesh();
    }
    if (!_tensors.empty()) {
        return _tensors.front()->mesh();
    }
    throw std::runtime_error("prism::output::FieldRegistry::mesh(): no fields are registered.");
}

VtuTimeWriter::VtuTimeWriter(VtuTimeWriterOptions options) : _options(std::move(options)) {
    if (_options.interval == 0) {
        throw std::invalid_argument(
            fmt::format("prism::output::VtuTimeWriter(): write interval must be positive, got {}",
                        _options.interval));
    }
}

VtuTimeWriter::VtuTimeWriter(VtuTimeWriter&&) noexcept = default;

auto VtuTimeWriter::operator=(VtuTimeWriter&&) noexcept -> VtuTimeWriter& = default;

VtuTimeWriter::~VtuTimeWriter() = default;

void VtuTimeWriter::add(SharedPtr<field::IScalar> f) {
    _registry.add(std::move(f));
}

void VtuTimeWriter::add(SharedPtr<field::IVector> f) {
    _registry.add(std::move(f));
}

void VtuTimeWriter::add(SharedPtr<field::ITensor> f) {
    _registry.add(std::move(f));
}

void VtuTimeWriter::advance(f64 dt) {
    ++_step;
    _time += dt;
}

void VtuTimeWriter::write() {
    if (_step == 0 && !_options.write_initial) {
        return;
    }
    if (_step % _options.interval != 0) {
        return;
    }
    if (_registry.empty()) {
        return;
    }

    if (!_vtu_writer) {
        _vtu_writer = std::make_unique<detail::VtuWriter>(_registry.mesh());
    }

    auto file_name = fmt::format("{}_{:04d}.vtu", _options.prefix, _step);
    auto pointers = fieldPointers(_registry);
    const VtuOptions vtu_options {
        .mode = _options.mode,
        .write_point_data = _options.write_point_data,
    };
    _vtu_writer->write(
        file_name, pointers.scalars, pointers.vectors, pointers.tensors, vtu_options);

    if (_options.write_pvd) {
        _pvd_entries.emplace_back(_time, file_name);
        flushPvd();
    }
}

void VtuTimeWriter::writeAndAdvance(f64 dt) {
    advance(dt);
    write();
}

void VtuTimeWriter::flushPvd() const {
    auto pvd_name = fmt::format("{}.pvd", _options.prefix);
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
