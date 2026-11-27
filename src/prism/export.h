#pragma once

#include <cstdint>
#include <memory>
#include <span>
#include <string>
#include <vector>

#include "prism/field/ifield.h"
#include "prism/types.h"

namespace prism::output {

/** @brief Supported VTU write modes. */
enum class WriteMode : std::uint8_t { Ascii, RawBinaryCompressed };

/** @brief VTU output options. */
struct VtuOptions {
    // VTU encoding.
    WriteMode mode {WriteMode::RawBinaryCompressed};

    // If true, also write IDW-interpolated PointData for each field.
    bool write_point_data {true};
};

/** @brief Time-series writer options. */
struct VtuTimeWriterOptions {
    // Base name for the output files (e.g. "result" produces result_0000.vtu).
    std::string prefix {"result"};

    // Write every N advances. Must be > 0.
    size_t interval {1};

    // If true, generate and update the .pvd collection.
    bool write_pvd {true};

    // VTU encoding.
    WriteMode mode {WriteMode::RawBinaryCompressed};

    // If true, also write IDW-interpolated PointData for each field.
    bool write_point_data {true};

    // If true, write() may emit the initial state at step 0.
    bool write_initial {false};
};

namespace detail {
class VtuWriter;
} // namespace detail

/** @brief Writes a scalar field to a VTU file (CellData and optionally PointData).
 * @param field The scalar field to export.
 * @param file_name The output file name.
 * @param options The VTU output options. */
void toVtu(const field::IScalar& field,
           const std::string& file_name,
           const VtuOptions& options = {});

/** @brief Writes a vector field to a VTU file (CellData and optionally PointData).
 * @param field The vector field to export.
 * @param file_name The output file name.
 * @param options The VTU output options. */
void toVtu(const field::IVector& field,
           const std::string& file_name,
           const VtuOptions& options = {});

/** @brief Writes a tensor field to a VTU file (CellData and optionally PointData).
 * @param field The tensor field to export.
 * @param file_name The output file name.
 * @param options The VTU output options. */
void toVtu(const field::ITensor& field,
           const std::string& file_name,
           const VtuOptions& options = {});

/** @brief Writes multiple fields into a single VTU file.
 *
 * All fields share the same mesh. CellData is always written. PointData is
 * computed via inverse distance weighting (IDW) from cell centres to vertices
 * when options.write_point_data is true.
 *
 * @param scalars Scalar fields to export.
 * @param vectors Vector fields to export.
 * @param tensors Tensor fields to export.
 * @param file_name The output file name.
 * @param options The VTU output options. */
void toVtu(std::span<const field::IScalar*> scalars,
           std::span<const field::IVector*> vectors,
           std::span<const field::ITensor*> tensors,
           const std::string& file_name,
           const VtuOptions& options = {});

/** @brief A registry of fields that can be exported together.
 *
 * Holds shared pointers to scalar, vector, and tensor fields. Fields registered
 * here are kept alive by the registry. Use exportAll() to write all registered
 * fields into a single VTU file. */
class FieldRegistry {
  public:
    /** @brief Register a scalar field. */
    void add(SharedPtr<field::IScalar> f);

    /** @brief Register a vector field. */
    void add(SharedPtr<field::IVector> f);

    /** @brief Register a tensor field. */
    void add(SharedPtr<field::ITensor> f);

    /** @brief Export all registered fields into a single VTU file.
     * @param file_name The output file name.
     * @param options The VTU output options. */
    void exportAll(const std::string& file_name, const VtuOptions& options = {}) const;

    /** @brief Returns true when no fields are registered. */
    auto empty() const noexcept -> bool;

    /** @brief Returns the mesh shared by the registered fields.
     * @throws std::runtime_error If the registry is empty. */
    auto mesh() const -> SharedPtr<mesh::PMesh>;

    auto scalars() const noexcept -> const std::vector<SharedPtr<field::IScalar>>& {
        return _scalars;
    }
    auto vectors() const noexcept -> const std::vector<SharedPtr<field::IVector>>& {
        return _vectors;
    }
    auto tensors() const noexcept -> const std::vector<SharedPtr<field::ITensor>>& {
        return _tensors;
    }

  private:
    std::vector<SharedPtr<field::IScalar>> _scalars;
    std::vector<SharedPtr<field::IVector>> _vectors;
    std::vector<SharedPtr<field::ITensor>> _tensors;
};

/** @brief Automates time-stepped VTU export with interval control and PVD generation.
 *
 * Fields are registered once and written automatically on matching steps.
 * Each write produces a zero-padded VTU file and updates a PVD collection
 * for ParaView time-series support.
 *
 * Call write() once before the first advance()/writeAndAdvance() to emit the
 * initial state when write_initial is enabled; the writer itself does not force
 * an initial write. */
class VtuTimeWriter {
  public:
    /** @brief Construct a time-stepped writer.
     * @param options Time-series and VTU options. */
    explicit VtuTimeWriter(VtuTimeWriterOptions options = {});

    VtuTimeWriter(const VtuTimeWriter&) = delete;
    auto operator=(const VtuTimeWriter&) -> VtuTimeWriter& = delete;
    VtuTimeWriter(VtuTimeWriter&&) noexcept;
    auto operator=(VtuTimeWriter&&) noexcept -> VtuTimeWriter&;
    ~VtuTimeWriter();

    /** @brief Register a scalar field. */
    void add(SharedPtr<field::IScalar> f);

    /** @brief Register a vector field. */
    void add(SharedPtr<field::IVector> f);

    /** @brief Register a tensor field. */
    void add(SharedPtr<field::ITensor> f);

    /** @brief Write the current fields if the step counter is a multiple of the interval. */
    void write();

    /** @brief Advance the step counter and accumulated time.
     * @param dt Time step size (added to accumulated time). */
    void advance(f64 dt = 0.0);

    /** @brief Advance then write. Convenience for solver loops.
     * @param dt Time step size. */
    void writeAndAdvance(f64 dt = 0.0);

    /** @brief The current accumulated time. */
    auto currentTime() const noexcept -> f64 { return _time; }

    /** @brief The current step counter (number of advances). */
    auto currentStep() const noexcept -> size_t { return _step; }

  private:
    void flushPvd() const;

    FieldRegistry _registry;
    size_t _step {0};
    f64 _time {0.0};
    VtuTimeWriterOptions _options;
    std::vector<Pair<f64, std::string>> _pvd_entries;
    std::unique_ptr<detail::VtuWriter> _vtu_writer;
};

} // namespace prism::output
