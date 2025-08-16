#pragma once

#include <algorithm>

#include "history.h"
#include "ifield.h"
#include "prism/exceptions.h"
#include "prism/gradient/gradient.h"
#include "prism/gradient/green_gauss.h"
#include "prism/gradient/least_squares.h"
#include "prism/log.h"
#include "prism/mesh/utilities.h"
#include "scalar_boundary.h"
#include "units.h"

/// TODO: we need to make sure that constructors are not leaving _data uninitialized, and we can
/// avoid checks for it later in member functions.
namespace prism::field {

namespace detail {
/// TODO: this function should be moved to a more appropriate place, like a utility file. Same for
/// coordToStr in ifield.h
auto inline coordToIndex(Coord coord) -> std::uint8_t {
    switch (coord) {
        case Coord::X: return 0;
        case Coord::Y: return 1;
        case Coord::Z: return 2;
        default: break;
    }
    throw std::invalid_argument("prism::field::coordToIndex(): Invalid coordinate value");
}
} // namespace detail

// forward declaration for GeneralScalar
template <typename Units, typename BHManagerSetter>
class GeneralScalar;

// forward declaration for GeneralVector
template <typename Component, typename BHManagerSetter>
class GeneralVector;

class UniformScalar : public IScalar {
  public:
    UniformScalar(std::string name, const SharedPtr<mesh::PMesh>& mesh, double value);

    auto valueAtCell(std::size_t cell_id) const -> double override;
    auto valueAtCell(const mesh::Cell& cell) const -> double override;
    auto valueAtFace(std::size_t face_id) const -> double override;
    auto valueAtFace(const mesh::Face& face) const -> double override;

    auto gradAtFace(const mesh::Face& face) -> Vector3d override;
    auto gradAtCell(const mesh::Cell& cell) -> Vector3d override;
    auto gradAtCellStored(const mesh::Cell& cell) const -> Vector3d override;

    template <IScalarBased ScalarField>
    auto operator*(const ScalarField& other) -> ScalarField;

    template <IVectorBased VectorField>
    auto operator*(const VectorField& other) -> VectorField;

  private:
    double _value {0.0};
};

template <typename Units, typename BHManagerSetter>
class GeneralScalar
    : public IScalar,
      public Units,
      public prism::boundary::BHManagerProvider<boundary::scalar::IScalarBoundaryHandler> {
  public:
    GeneralScalar(std::string name, const SharedPtr<mesh::PMesh>& mesh, f64 value);

    GeneralScalar(std::string name, const SharedPtr<mesh::PMesh>& mesh, f64 value, Coord coord);

    GeneralScalar(std::string name, const SharedPtr<mesh::PMesh>& mesh, VectorXd data);

    GeneralScalar(std::string name,
                  const SharedPtr<mesh::PMesh>& mesh,
                  VectorXd data,
                  Coord coord);

    GeneralScalar(std::string name,
                  const SharedPtr<mesh::PMesh>& mesh,
                  VectorXd data,
                  VectorXd face_data);

    GeneralScalar(std::string name,
                  const SharedPtr<mesh::PMesh>& mesh,
                  VectorXd data,
                  VectorXd face_data,
                  Coord coord);

    GeneralScalar() = delete;
    GeneralScalar(const GeneralScalar&) = default;
    GeneralScalar(GeneralScalar&&) = default;
    auto operator=(const GeneralScalar&) -> GeneralScalar& = default;
    auto operator=(GeneralScalar&&) -> GeneralScalar& = default;
    ~GeneralScalar() override = default;

    auto values() const -> const VectorXd&;
    auto values() -> VectorXd&;

    auto coord() const noexcept -> Optional<Coord> override;
    auto hasFaceValues() const -> bool override;
    void setFaceValues(VectorXd values);
    void clearFaceValues();

    auto valueAtCell(std::size_t cell_id) const -> double override;
    auto valueAtCell(const mesh::Cell& cell) const -> double override;
    auto valueAtFace(std::size_t face_id) const -> double override;
    auto valueAtFace(const mesh::Face& face) const -> double override;

    auto gradAtFace(const mesh::Face& face) -> Vector3d override;
    auto gradAtCell(const mesh::Cell& cell) -> Vector3d override;
    auto gradAtCellStored(const mesh::Cell& cell) const -> Vector3d override;

    void update(VectorXd values);

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

    auto clone() const -> GeneralScalar;

    auto operator[](std::size_t i) const -> double;
    auto operator[](std::size_t i) -> double&;

  protected:
    auto valueAtInteriorFace(const mesh::Face& face) const -> double;
    auto valueAtBoundaryFace(const mesh::Face& face) const -> double;

  private:
    void setGradScheme();
    void addDefaultBoundaryHandlers();

    SharedPtr<VectorXd> _cell_values = nullptr;
    SharedPtr<HistoryManager> _history_manager = nullptr;

    /// TODO: _face_data should not include empty faces
    SharedPtr<VectorXd> _face_values = nullptr;
    SharedPtr<gradient::IGradient> _grad_scheme = nullptr;

    const IVector* _parent = nullptr;
    Optional<Coord> _coord = NullOption;
    BHManagerSetter _setter;
};

class ScalarBHManagerSetter {
  public:
    using IScalarBHManager =
        prism::boundary::BoundaryHandlersManager<boundary::scalar::IScalarBoundaryHandler>;

    static void set(IScalarBHManager& manager);
};

using Scalar = GeneralScalar<units::Measurable, ScalarBHManagerSetter>;

template <IScalarBased ScalarField>
auto UniformScalar::operator*(const ScalarField& other) -> ScalarField {
    /// TODO: this is not a good way to compare meshes, we should use a mesh equality operator.
    if (this->mesh()->cellCount() != other.mesh()->cellCount()) {
        throw std::runtime_error(
            fmt::format("UniformScalar::operator*(): cannot multiply scalar field '{}' with "
                        "scalar field '{}', because they are defined on different meshes.",
                        this->name(),
                        other.name()));
    }
    /// TODO: this ignores the units of the scalar field
    return ScalarField(fmt::format("{}{}", this->name(), other.name()),
                       this->mesh(),
                       this->_value * other.values());
}

template <IVectorBased VectorField>
auto UniformScalar::operator*(const VectorField& other) -> VectorField {
    /// TODO: this function ignores the units of the scalar field, which is not ideal, and returns
    /// a vector of components that preserves `other`'s units. We need to fix this.
    if (this->mesh()->cellCount() != other.mesh()->cellCount()) {
        throw std::runtime_error(
            fmt::format("UniformScalar::operator*(): cannot multiply scalar field '{}' with "
                        "vector field '{}', because they are defined on different meshes.",
                        this->name(),
                        other.name()));
    }

    auto output_name = fmt::format("{}{}", this->name(), other.name());
    std::vector<Vector3d> face_values(other.mesh()->faceCount(), Vector3d::Zero());
    VectorXd face_flux_values = VectorXd::Zero(other.mesh()->faceCount());

    for (const auto& face : this->mesh()->interiorFaces()) {
        const Vector3d face_value = other.valueAtFace(face) * this->_value;
        face_values[face.id()] = face_value;
        face_flux_values[face.id()] = face_value.dot(face.areaVector());
        ;
    }
    for (const auto& patch : this->mesh()->boundaryPatches()) {
        if (patch.isEmpty()) {
            continue; // skip empty patches
        }
        for (const auto& face_id : patch.facesIds()) {
            Vector3d result = other.valueAtFace(face_id) * this->_value;
            const auto& face = this->mesh()->face(face_id);
            face_values[face_id] = result;
            face_flux_values[face_id] = result.dot(face.areaVector());
        }
    }

    using Component = typename VectorField::ComponentType;
    std::array<Component, 3> components = {
        Component(output_name + "_x", this->mesh(), this->_value * other.x().values().array()),
        Component(output_name + "_y", this->mesh(), this->_value * other.y().values().array()),
        Component(output_name + "_z", this->mesh(), this->_value * other.z().values().array())};

    auto U =
        VectorField(fmt::format("{}{}", this->name(), other.name()), this->mesh(), components);
    U.setFaceValues(std::move(face_values));
    U.setFaceFluxValues(std::move(face_flux_values));
    return U;
}

template <typename Units, typename BHManagerSetter>
GeneralScalar<Units, BHManagerSetter>::GeneralScalar(std::string name,
                                                     const SharedPtr<mesh::PMesh>& mesh,
                                                     double value)
    : IScalar(std::move(name), mesh),
      _cell_values(std::make_shared<VectorXd>(VectorXd::Ones(mesh->cellCount()) * value)) {
    log::debug("Creating scalar field: '{}' with double value = {}", this->name(), value);
    addDefaultBoundaryHandlers();
    setGradScheme();
}

template <typename Units, typename BHManagerSetter>
GeneralScalar<Units, BHManagerSetter>::GeneralScalar(std::string name,
                                                     const SharedPtr<mesh::PMesh>& mesh,
                                                     double value,
                                                     Coord coord)
    : IScalar(std::move(name), mesh),
      _cell_values(std::make_shared<VectorXd>(VectorXd::Ones(mesh->cellCount()) * value)),
      _coord(coord) {
    log::debug("Creating scalar field: '{}' (as {}-coordinate) with double value = {}",
               this->name(),
               coordToStr(coord),
               value);

    addDefaultBoundaryHandlers();
    setGradScheme();
}

template <typename Units, typename BHManagerSetter>
GeneralScalar<Units, BHManagerSetter>::GeneralScalar(std::string name,
                                                     const SharedPtr<mesh::PMesh>& mesh,
                                                     VectorXd data)
    : IScalar(std::move(name), mesh), _cell_values(std::make_shared<VectorXd>(std::move(data))) {
    if (_cell_values->size() != mesh->cellCount()) {
        throw std::runtime_error(
            fmt::format("field::GeneralScalar() cannot create a scalar field '{}' given a "
                        "vector that has a different size than mesh's cell count.",
                        this->name()));
    }

    log::debug("Creating scalar field: '{}' with a cell vector data of size = {}",
               this->name(),
               _cell_values->size());

    addDefaultBoundaryHandlers();
    setGradScheme();
}

template <typename Units, typename BHManagerSetter>
GeneralScalar<Units, BHManagerSetter>::GeneralScalar(std::string name,
                                                     const SharedPtr<mesh::PMesh>& mesh,
                                                     VectorXd data,
                                                     Coord coord)
    : IScalar(std::move(name), mesh),
      _cell_values(std::make_shared<VectorXd>(std::move(data))),
      _coord(coord) {
    if (_cell_values->size() != mesh->cellCount()) {
        throw std::runtime_error(fmt::format(
            "field::Scalar() cannot create a scalar field '{}' given a vector that has a "
            "different size than mesh's cell count.",
            this->name()));
    }

    log::debug(
        "Creating scalar field: '{}' (as {}-coordinate) with a cell vector data of size = {}",
        this->name(),
        coordToStr(coord),
        _cell_values->size());
    addDefaultBoundaryHandlers();
    setGradScheme();
}

template <typename Units, typename BHManagerSetter>
GeneralScalar<Units, BHManagerSetter>::GeneralScalar(std::string name,
                                                     const SharedPtr<mesh::PMesh>& mesh,
                                                     VectorXd data,
                                                     VectorXd face_data)
    : IScalar(std::move(name), mesh),
      _cell_values(std::make_shared<VectorXd>(std::move(data))),
      _face_values(std::make_shared<VectorXd>(std::move(face_data))) {
    if (_cell_values->size() != mesh->cellCount()) {
        throw std::runtime_error(fmt::format(
            "field::Scalar() cannot create a scalar field '{}' given a vector that has a "
            "different size than mesh's cell count.",
            this->name()));
    }

    if (_face_values->size() != mesh->faceCount()) {
        throw std::runtime_error(
            fmt::format("field::Scalar() cannot create a scalar field '{}' given a face data "
                        "vector that has a different size than mesh's faces count.",
                        this->name()));
    }

    log::debug(
        "Creating scalar field: '{}' with a cell data vector of size = {} and face data "
        "vector of size = {}",
        this->name(),
        _cell_values->size(),
        _face_values->size());

    addDefaultBoundaryHandlers();
    setGradScheme();
}

template <typename Units, typename BHManagerSetter>
GeneralScalar<Units, BHManagerSetter>::GeneralScalar(std::string name,
                                                     const SharedPtr<mesh::PMesh>& mesh,
                                                     VectorXd data,
                                                     VectorXd face_data,
                                                     Coord coord)
    : IScalar(std::move(name), mesh),
      _cell_values(std::make_shared<VectorXd>(std::move(data))),
      _face_values(std::make_shared<VectorXd>(std::move(face_data))),
      _coord(coord) {
    if (_cell_values->size() != mesh->cellCount()) {
        throw std::runtime_error(fmt::format(
            "field::Scalar() cannot create a scalar field '{}' given a vector that has a "
            "different size than mesh's cell count.",
            this->name()));
    }

    if (_face_values->size() != mesh->faceCount()) {
        throw std::runtime_error(
            fmt::format("field::Scalar() cannot create a scalar field '{}' given a face data "
                        "vector that has a different size than mesh's faces count.",
                        this->name()));
    }

    log::debug(
        "Creating scalar field: '{}' (as {}-coordinate) with a cell data vector of size = {} "
        "and "
        "face data vector of size = {}",
        this->name(),
        coordToStr(coord),
        _cell_values->size(),
        _face_values->size());

    addDefaultBoundaryHandlers();
    setGradScheme();
}

/// TODO: delete this in favor of updateFaceValues()
template <typename Units, typename BHManagerSetter>
void GeneralScalar<Units, BHManagerSetter>::setFaceValues(VectorXd values) {
    if (values.size() != mesh()->faceCount()) {
        throw std::runtime_error(
            fmt::format("prism::field::GeneralScalar<Units, BHManagerProvider, "
                        "BHManagerSetter>::setFaceValues(): cannot set face values for scalar "
                        "field {}, to a face data vector having a different size that "
                        "field's faces count.",
                        name()));
    }

    if (hasFaceValues()) {
        log::debug(
            "GeneralScalar<Units, BHManagerSetter>::setFaceValues() was called for field "
            "'{}', which already has face values set.",
            name());
    }

    log::debug("Setting face values for field '{}'", name());
    _face_values = std::make_shared<VectorXd>(std::move(values));
}

template <typename Units, typename BHManagerSetter>
void GeneralScalar<Units, BHManagerSetter>::clearFaceValues() {
    if (_face_values == nullptr) {
        log::warn(
            "GeneralScalar<Units, BHManagerSetter>::clearFaceValues() was called for field "
            "'{}', but the face data is not initialized.",
            name());
        return;
    }
    _face_values = nullptr;
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::values() const -> const VectorXd& {
    if (_cell_values == nullptr) {
        throw std::runtime_error(
            fmt::format("prism::field::GeneralScalar<Units, BHManagerSetter>::values() was "
                        "called for field `{}`, but the data is not initialized.",
                        name()));
    }
    return *_cell_values;
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::values() -> VectorXd& {
    if (_cell_values == nullptr) {
        throw std::runtime_error(
            fmt::format("prism::field::GeneralScalar<Units, BHManagerSetter>::values() was "
                        "called for field `{}`, but the data is not initialized.",
                        name()));
    }
    return *_cell_values;
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::coord() const noexcept -> Optional<Coord> {
    return _coord;
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::hasFaceValues() const -> bool {
    return _face_values != nullptr;
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::valueAtCell(const mesh::Cell& cell) const -> double {
    return valueAtCell(cell.id());
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::valueAtCell(std::size_t cell_id) const -> double {
    assert(_cell_values != nullptr);       // NOLINT
    assert(cell_id < mesh()->cellCount()); // NOLINT
    return (*_cell_values)[cell_id];
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::valueAtFace(std::size_t face_id) const -> double {
    if (hasFaceValues()) {
        // Face data were calculated already, just return the value (as in Rhie-Chow
        // correction).
        return (*_face_values)[face_id];
    }

    const auto& face = mesh()->face(face_id);

    if (face.isInterior()) {
        return valueAtInteriorFace(face);
    }

    return valueAtBoundaryFace(face);
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::valueAtFace(const mesh::Face& face) const -> double {
    return valueAtFace(face.id());
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::valueAtInteriorFace(const mesh::Face& face) const
    -> double {
    assert(face.isInterior()); // NOLINT
    const auto& owner = mesh()->cell(face.owner());
    const auto& neighbor = mesh()->cell(face.neighbor().value());

    const auto gc = mesh::geometricWeight(owner, neighbor, face);
    double val = gc * (*_cell_values)[owner.id()];
    val += (1 - gc) * (*_cell_values)[neighbor.id()];

    return val;
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::valueAtBoundaryFace(const mesh::Face& face) const
    -> double {
    // if this scalar field is a component of a vector field, we need to return the value of the
    // parent field at the face in the coordinate direction.
    if ((_parent != nullptr) && (_coord.has_value())) {
        auto idx = detail::coordToIndex(_coord.value());
        return _parent->valueAtFace(face)[idx];
    }
    const auto& patch = mesh()->boundaryPatch(face);
    const auto& bc = patch.getBoundaryCondition(name());

    auto handler = this->boundaryHandlersManager().getHandler(bc.kindString());

    if (handler == nullptr) {
        throw error::NonImplementedBoundaryCondition(
            fmt::format("prism::field::GeneralScalar<Units, BHManagerProvider, "
                        "BHManagerSetter>::valueAtBoundaryFace() for field `{}`",
                        name()),
            patch.name(),
            bc.kindString());
    }

    return handler->get(*this, face);
}

template <typename Units, typename BHManagerSetter>
void GeneralScalar<Units, BHManagerSetter>::addDefaultBoundaryHandlers() {
    _setter.set(this->boundaryHandlersManager());
}

template <typename Units, typename BHManagerSetter>
void GeneralScalar<Units, BHManagerSetter>::setGradScheme(
    const SharedPtr<gradient::IGradient>& grad_scheme) {
    if (grad_scheme == nullptr) {
        throw std::runtime_error(
            "GeneralScalar::setGradScheme() was given a null gradient scheme pointer");
    }
    _grad_scheme = grad_scheme;
}

template <typename Units, typename BHManagerSetter>
void GeneralScalar<Units, BHManagerSetter>::setGradScheme() {
    // did user specify gradient scheme for the field in `fields.json`?
    auto field_infos = this->mesh()->fieldsInfo();
    auto it = std::find_if(field_infos.begin(), field_infos.end(), [this](const auto& fi) {
        return fi.name() == this->name() && fi.gradScheme().has_value();
    });

    /// TODO: this is buggy, it doesn't find the grad scheme defined in fields.json for the
    // field, also does not consider vector fields.
    if (it == field_infos.end()) {
        log::debug(
            "GeneralScalar::setGradScheme(): couldn't find a specified gradient scheme for "
            "field `{}` in `fields.json`, setting the gradient scheme to least squares.",
            this->name());

        _grad_scheme = std::make_shared<gradient::LeastSquares>(this->mesh());
        return;
    }

    auto grad_scheme_name = it->gradScheme().value();

    if (grad_scheme_name == "green-gauss" || grad_scheme_name == "greenGauss") {
        log::debug(
            "GeneralScalar::setGradScheme(): setting the gradient scheme to Green-Gauss for "
            "field `{}`",
            this->name());
        _grad_scheme = std::make_shared<gradient::GreenGauss>(this->mesh());
        return;
    }

    /// TODO: if scalar field is a component of a vector field, prism should get the grad scheme
    /// type from it's parent field, and we should also implement setGradScheme(scheme) for vector
    /// field type.

    log::debug(
        "GeneralScalar::setGradScheme(): setting the gradient scheme to Least-Squares "
        "for field `{}`",
        this->name());

    _grad_scheme = std::make_shared<gradient::LeastSquares>(this->mesh());
}

template <typename Units, typename BHManagerSetter>
void GeneralScalar<Units, BHManagerSetter>::setHistorySize(std::size_t num_time_steps) {
    // If the requested history size is greater than zero, we need to ensure
    // a history manager is available and configured to the correct size.
    if (num_time_steps > 0) {
        // Check if a history manager is already allocated.
        if (!_history_manager) {
            // If not, create a new shared history manager.
            _history_manager = std::make_shared<HistoryManager>(num_time_steps);
        } else {
            // If a manager already exists, resize it. This preserves existing history.
            _history_manager->resize(num_time_steps);
        }
    } else {
        // If the requested size is 0, we release the history manager.
        _history_manager.reset();
    }
}

template <typename Units, typename BHManagerSetter>
void GeneralScalar<Units, BHManagerSetter>::update(VectorXd values) {
    if (_history_manager) {
        _history_manager->update(*_cell_values);
    }
    *_cell_values = std::move(values);

    if (_face_values != nullptr) {
        clearFaceValues();
    }
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::prevValues() const -> Optional<VectorXd> {
    if (_history_manager) {
        return _history_manager->prevValues();
    }
    return NullOption;
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::prevPrevValues() const -> Optional<VectorXd> {
    if (_history_manager) {
        return _history_manager->prevPrevValues();
    }
    return NullOption;
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::getHistory(std::size_t index) const
    -> Optional<VectorXd> {
    /**
     * @brief Retrieves historical field values at a specific time step.
     *
     * This method provides access to the field's values from previous time steps
     * stored in the history manager. The `index` parameter specifies how many
     * steps back in time to retrieve the data.
     *
     * @param index The zero-based index representing the number of time steps
     *              ago to retrieve the values. An index of 0 refers to the
     *              immediately previous time step (t-1), 1 refers to two time
     *              steps ago (t-2), and so on.
     * @return An `Optional<VectorXd>` containing the historical field values
     *         if available, or `NullOption` if the history manager is not
     *         enabled or the requested index is out of bounds.
     *
     * @note The history manager must be enabled and configured with a sufficient
     *       size using `setHistorySize()` for this method to return meaningful
     *       data. The current field values (t) are *not* part of the history
     *       retrieved by this method; they are accessed via `values()`.
     *
     * @par Examples:
     * @code
     * // Assuming 'T' is a GeneralScalar instance with history enabled
     * // and updated multiple times.
     *
     * // Get values from the immediately previous time step (t-1)
     * Optional<VectorXd> t_minus_1 = T.getHistory(0);
     * if (t_minus_1.has_value()) {
     *     // Use t_minus_1.value()
     * }
     *
     * // Get values from two time steps ago (t-2)
     * Optional<VectorXd> t_minus_2 = T.getHistory(1);
     * if (t_minus_2.has_value()) {
     *     // Use t_minus_2.value()
     * }
     *
     * // Attempt to get values from an index beyond the stored history size
     * Optional<VectorXd> out_of_bounds = T.getHistory(100);
     * if (!out_of_bounds.has_value()) {
     *     // Handle case where history is not available for this index
     * }
     * @endcode
     */

    if (_history_manager) {
        return _history_manager->valuesAt(index);
    }
    return NullOption;
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::gradAtFace(const mesh::Face& face) -> Vector3d {
    return _grad_scheme->gradAtFace(face, *this);
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::gradAtCell(const mesh::Cell& cell) -> Vector3d {
    return _grad_scheme->gradAtCell(cell, *this);
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::gradAtCellStored(const mesh::Cell& cell) const
    -> Vector3d {
    return _grad_scheme->gradAtCellStored(cell, *this);
}

template <typename Units, typename BHManagerSetter>
template <typename Func>
void GeneralScalar<Units, BHManagerSetter>::updateInteriorFaces(Func func) {
    if (!hasFaceValues()) {
        /// TODO: implement initFaceDataVector() to initialize face data vector
        VectorXd face_values(this->mesh()->faceCount());
        face_values.setZero();

        for (const auto& patch : this->mesh()->boundaryPatches()) {
            if (patch.isEmpty()) {
                continue; // skip empty patches
            }
            for (const auto& face_id : patch.facesIds()) {
                face_values[face_id] = valueAtFace(face_id);
            }
        }
        _face_values = std::make_shared<VectorXd>(std::move(face_values));
    }

    for (const auto& face : this->mesh()->interiorFaces()) {
        (*_face_values)[face.id()] = func(face);
    }
}

template <typename Units, typename BHManagerSetter>
template <typename Func>
void GeneralScalar<Units, BHManagerSetter>::updateFaces(Func func) {
    updateInteriorFaces(func);

    for (const auto& patch : this->mesh()->boundaryPatches()) {
        if (patch.isEmpty()) {
            continue; // skip empty patches
        }
        for (const auto& face_id : patch.facesIds()) {
            const auto& face = this->mesh()->face(face_id);
            (*_face_values)[face_id] = func(face);
        }
    }
}

template <typename Units, typename BHManagerSetter>
template <typename Func>
void GeneralScalar<Units, BHManagerSetter>::updateCells(Func func) {
    /// TODO: test this.
    for (const auto& cell : this->mesh()->cells()) {
        (*_cell_values)[cell.id()] = func(cell);
    }
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::clone() const -> GeneralScalar {
    /// NOTE: cloned field is parentless.
    if (this->coord().has_value()) {
        auto clone =
            GeneralScalar(this->name(), this->mesh(), this->values(), this->coord().value());
        if (hasFaceValues()) {
            clone.setFaceValues(*_face_values);
        }
        return clone;
    }
    auto clone = GeneralScalar(this->name(), this->mesh(), this->values());
    if (hasFaceValues()) {
        clone.setFaceValues(*_face_values);
    }
    return clone;
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::operator[](std::size_t i) const -> double {
    if (_cell_values == nullptr) {
        throw std::runtime_error(
            fmt::format("prism::field::GeneralScalar<Units, BHManagerSetter>::operator[]() was "
                        "called for field `{}`, but the data is not initialized.",
                        name()));
    }
    if (i >= _cell_values->size()) {
        throw std::out_of_range(
            fmt::format("prism::field::GeneralScalar<Units, BHManagerSetter>::operator[]() was "
                        "called for field `{}`, but the index {} is out of range (size = {}).",
                        name(),
                        i,
                        _cell_values->size()));
    }
    return (*_cell_values)[i];
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::operator[](std::size_t i) -> double& {
    if (_cell_values == nullptr) {
        throw std::runtime_error(
            fmt::format("prism::field::GeneralScalar<Units, BHManagerSetter>::operator[]() was "
                        "called for field `{}`, but the data is not initialized.",
                        name()));
    }

    if (i >= _cell_values->size()) {
        throw std::out_of_range(
            fmt::format("prism::field::GeneralScalar<Units, BHManagerSetter>::operator[]() was "
                        "called for field `{}`, but the index {} is out of range (size = {}).",
                        name(),
                        i,
                        _cell_values->size()));
    }
    return (*_cell_values)[i];
}

} // namespace prism::field
