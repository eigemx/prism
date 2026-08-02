#pragma once

#include "history.h"
#include "ifield.h"
#include "prism/gradient/gradient.h"
#include "scalar_boundary.h"
#include "units.h"
#include "utilities.h"

/// TODO: we need to make sure that constructors are not leaving _data uninitialized, and we can
/// avoid checks for it later in member functions.
namespace prism::field {

class Scalar : public IScalar,
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

    auto valueAtCell(size_t cell_id) const -> f64 override;
    auto valueAtCell(const mesh::Cell& cell) const -> f64 override;
    auto valueAtFace(size_t face_id) const -> f64 override;
    auto valueAtFace(const mesh::Face& face) const -> f64 override;

    auto gradAtFace(const mesh::Face& face) const -> Vector3d override;
    auto gradAtCell(const mesh::Cell& cell) const -> Vector3d override;

    void computeAllCellGradients();

    auto eventNo() const noexcept -> size_t override;

    void update(VectorXd values);
    void updatePrevTimeSteps();

    template <FaceMapper Func>
    void mapInteriorFaceValues(Func func);

    template <FaceMapper Func>
    void mapFaceValues(Func func);

    template <CellMapper Func>
    void mapCellValues(Func func);

    auto clone() const -> Scalar;

    void setGradScheme(const SharedPtr<gradient::IGradient>& grad_scheme);

    void setHistorySize(size_t num_time_steps);
    auto prevValues() const -> Optional<Ref<const VectorXd>>;
    auto prevPrevValues() const -> Optional<Ref<const VectorXd>>;
    auto getHistory(size_t index) const -> Optional<Ref<const VectorXd>>;

    auto operator[](size_t i) const -> f64;
    auto operator[](size_t i) -> f64&;

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

    // Gradient scheme that caches per-field cell gradients. mutable: gradient reads are const
    // with respect to the field values, only the scheme's internal cache is updated.
    mutable SharedPtr<gradient::IGradient> _grad_scheme = nullptr;

    // Revision of the cell values; bumped on every mutation. Cached gradients are invalidated
    // when this changes.
    size_t _event_no {0};

    // This should have a value only when the object is a component of a field::Vector instance.
    Optional<VectorCoord> _coord = NullOption;
};

template <FaceMapper Func>
void Scalar::mapInteriorFaceValues(Func func) {
    if (!hasFaceValues()) {
        VectorXd face_values(mesh()->faceCount());

        for (const auto& patch : mesh()->boundaryPatches()) {
            if (patch.isEmpty()) {
                continue;
            }
            for (const auto& face_id : patch.facesIds()) {
                face_values[face_id] = valueAtFace(face_id);
            }
        }
        _face_values = std::move(face_values);
    }

    for (const auto& face : mesh()->interiorFaces()) {
        _face_values[face.id()] = func(face);
    }
}

template <FaceMapper Func>
void Scalar::mapFaceValues(Func func) {
    mapInteriorFaceValues(func);

    for (const auto& patch : mesh()->boundaryPatches()) {
        if (patch.isEmpty()) {
            continue;
        }
        for (const auto& face_id : patch.facesIds()) {
            const auto& face = mesh()->face(face_id);
            _face_values[face_id] = func(face);
        }
    }
}

template <CellMapper Func>
void Scalar::mapCellValues(Func func) {
    for (const auto& cell : mesh()->cells()) {
        _cell_values[cell.id()] = func(cell);
    }
    ++_event_no;
}

// Element-wise field arithmetic

inline auto mul(const SharedPtr<Scalar>& a, const SharedPtr<Scalar>& b) -> VectorXd {
    return a->values().cwiseProduct(b->values());
}

inline auto divide(const SharedPtr<Scalar>& a, const SharedPtr<Scalar>& b) -> VectorXd {
    return a->values().cwiseQuotient(b->values());
}

inline auto add(const SharedPtr<Scalar>& a, const SharedPtr<Scalar>& b) -> VectorXd {
    return a->values() + b->values();
}

inline auto sub(const SharedPtr<Scalar>& a, const SharedPtr<Scalar>& b) -> VectorXd {
    return a->values() - b->values();
}

// Mixed: SharedPtr with VectorXd

inline auto mul(const SharedPtr<Scalar>& a, const VectorXd& b) -> VectorXd {
    return a->values().cwiseProduct(b);
}

inline auto mul(const VectorXd& a, const SharedPtr<Scalar>& b) -> VectorXd {
    return a.cwiseProduct(b->values());
}

inline auto divide(const VectorXd& a, const SharedPtr<Scalar>& b) -> VectorXd {
    return a.cwiseQuotient(b->values());
}

inline auto add(const VectorXd& a, const SharedPtr<Scalar>& b) -> VectorXd {
    return a + b->values();
}

inline auto sub(const VectorXd& a, const SharedPtr<Scalar>& b) -> VectorXd {
    return a - b->values();
}

// Variadic mul: mul(a, b, c, ...) for 3+ fields

template <typename... Rest>
inline auto mul(const SharedPtr<Scalar>& a,
                const SharedPtr<Scalar>& b,
                const SharedPtr<Scalar>& c,
                const Rest&... rest) -> VectorXd {
    return mul(a, mul(b, c, rest...));
}

} // namespace prism::field
