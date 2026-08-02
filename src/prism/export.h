#pragma once

#include <cstdint>
#include <span>
#include <string>
#include <vector>

#include "prism/field/ifield.h"
#include "prism/types.h"

namespace prism::output {

/** @brief Supported VTU write modes. */
enum class WriteMode : std::uint8_t { Ascii, RawBinaryCompressed };

/** @brief Writes a scalar field to a VTU file (both CellData and PointData).
 * @param field The scalar field to export.
 * @param file_name The output file name.
 * @param mode The VTU write mode. */
void toVtu(const field::IScalar& field,
           const std::string& file_name,
           WriteMode mode = WriteMode::RawBinaryCompressed);

/** @brief Writes a vector field to a VTU file (both CellData and PointData).
 * @param field The vector field to export.
 * @param file_name The output file name.
 * @param mode The VTU write mode. */
void toVtu(const field::IVector& field,
           const std::string& file_name,
           WriteMode mode = WriteMode::RawBinaryCompressed);

/** @brief Writes a tensor field to a VTU file (both CellData and PointData).
 * @param field The tensor field to export.
 * @param file_name The output file name.
 * @param mode The VTU write mode. */
void toVtu(const field::ITensor& field,
           const std::string& file_name,
           WriteMode mode = WriteMode::RawBinaryCompressed);

/** @brief Writes multiple fields into a single VTU file.
 *
 * All fields share the same mesh. Each field is written as both CellData and
 * PointData. PointData is computed via inverse distance weighting (IDW) from
 * cell centres to vertices.
 *
 * @param scalars Scalar fields to export.
 * @param vectors Vector fields to export.
 * @param tensors Tensor fields to export.
 * @param file_name The output file name.
 * @param mode The VTU write mode. */
void toVtu(std::span<const field::IScalar*> scalars,
           std::span<const field::IVector*> vectors,
           std::span<const field::ITensor*> tensors,
           const std::string& file_name,
           WriteMode mode = WriteMode::RawBinaryCompressed);

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
     * @param mode The VTU write mode. */
    void exportAll(const std::string& file_name,
                   WriteMode mode = WriteMode::RawBinaryCompressed) const;

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
 * for ParaView time-series support. */
class TimeWriter {
  public:
    /** @brief Construct a time-stepped writer.
     * @param prefix Base name for output files (e.g., "result" produces result_0000.vtu).
     * @param write_interval Write every N steps. 1 writes every step. Must be > 0.
     * @param mode The VTU write mode.
     * @param write_pvd If true, generate a .pvd file for ParaView. */
    TimeWriter(std::string prefix = "result",
               size_t write_interval = 1,
               WriteMode mode = WriteMode::RawBinaryCompressed,
               bool write_pvd = true);

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
    std::string _prefix;
    size_t _interval;
    WriteMode _mode;
    bool _write_pvd;
    std::vector<Pair<f64, std::string>> _pvd_entries;
};

} // namespace prism::output
