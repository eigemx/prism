#pragma once

#include "history.h"
#include "ifield.h"
#include "prism/gradient/gradient.h"
#include "scalar_boundary.h"
#include "units.h"

/// TODO: we need to make sure that constructors are not leaving _data uninitialized, and we can
/// avoid checks for it later in member functions.
namespace prism::field {

namespace detail {
/// TODO: this function should be moved to a more appropriate place, like a utility file. Same for
/// coordToStr in ifield.h
auto inline coordToIndex(VectorCoord coord) -> std::uint8_t {
    switch (coord) {
        case VectorCoord::X: return 0;
        case VectorCoord::Y: return 1;
        case VectorCoord::Z: return 2;
        default: break;
    }
    throw std::invalid_argument("prism::field::coordToIndex(): Invalid coordinate value");
}
} // namespace detail

class Scalar
    : public IScalar,
      public units::Measurable,
      public prism::boundary::BHManagerProvider<boundary::scalar::IScalarBoundaryHandler> {
  public:
    // Scalar field based on f64 value and an optional coordinate (if vector sub-field)
    Scalar(std::string name, const SharedPtr<mesh::PMesh>& mesh, f64 value);
    Scalar(std::string name, const SharedPtr<mesh::PMesh>& mesh, f64 value, VectorCoord coord);

    // Scalar field based on VectorXd values and an optional coordinate.
    Scalar(std::string name, const SharedPtr<mesh::PMesh>& mesh, VectorXd data);
    Scalar(std::string name,
           const SharedPtr<mesh::PMesh>& mesh,
           VectorXd values,
           VectorCoord coord);

    // Scalar field based on VectorXd values, with field face values and an optional coordinate.
    Scalar(std::string name,
           const SharedPtr<mesh::PMesh>& mesh,
           VectorXd values,
           VectorXd face_values);
    Scalar(std::string name,
           const SharedPtr<mesh::PMesh>& mesh,
           VectorXd values,
           VectorXd face_values,
           VectorCoord coord);

    // default constructors
    Scalar(const Scalar&) = default;
    Scalar(Scalar&&) = default;
    auto operator=(const Scalar&) -> Scalar& = default;
    auto operator=(Scalar&&) -> Scalar& = default;
    ~Scalar() override = default;

    auto values() const -> const VectorXd&;
    auto values() -> VectorXd&;

    auto coord() const noexcept -> Optional<VectorCoord> override;

    auto hasFaceValues() const -> bool override;
    void setFaceValues(VectorXd values);
    void clearFaceValues();

    auto valueAtCell(std::size_t cell_id) const -> f64 override;
    auto valueAtCell(const mesh::Cell& cell) const -> f64 override;
    auto valueAtFace(std::size_t face_id) const -> f64 override;
    auto valueAtFace(const mesh::Face& face) const -> f64 override;

    auto gradAtFace(const mesh::Face& face) -> Vector3d override;
    auto gradAtCell(const mesh::Cell& cell) -> Vector3d override;
    auto gradAtCellStored(const mesh::Cell& cell) const -> Vector3d override;

    /// TODO: overload for the case of zero arguments, to update with current values.
    void update(VectorXd values);
    void updatePrevTimeSteps();

    template <typename Func>
    void updateInteriorFaces(Func func);

    template <typename Func>
    void updateFaces(Func func);

    template <typename Func>
    void updateCells(Func func);

    void setGradScheme(const SharedPtr<gradient::IGradient>& grad_scheme);

    void setHistorySize(std::size_t num_time_steps);
    auto prevValues() const -> Optional<VectorXd>;
    auto prevPrevValues() const -> Optional<VectorXd>;
    auto getHistory(std::size_t index) const -> Optional<VectorXd>;

    auto operator[](std::size_t i) const -> f64;
    auto operator[](std::size_t i) -> f64&;

  protected:
    auto valueAtInteriorFace(const mesh::Face& face) const -> f64;
    auto valueAtBoundaryFace(const mesh::Face& face) const -> f64;

    void setGradScheme();
    virtual void addDefaultBoundaryHandlers();
    virtual void setUnits();

  private:
    Scalar(std::string name,
           const SharedPtr<mesh::PMesh>& mesh,
           VectorXd cell_values,
           VectorXd face_values,
           Optional<VectorCoord> coord);

    void validateCellValues(const VectorXd& values) const;
    void validateFaceValues(const VectorXd& values) const;
    void logConstruction() const;

    VectorXd _cell_values;
    HistoryManager _history_manager;

    /// TODO: _face_data should not include empty faces
    VectorXd _face_values;
    SharedPtr<gradient::IGradient> _grad_scheme = nullptr;

    // This should have a value only when the object is a component of a field::Vector instance.
    Optional<VectorCoord> _coord = NullOption;
};
} // namespace prism::field
