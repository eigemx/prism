#include "scalar.h"

#include <algorithm>
#include <cassert>

#include "prism/exceptions.h"
#include "prism/gradient/green_gauss.h"
#include "prism/gradient/least_squares.h"
#include "prism/log.h"
#include "prism/mesh/utilities.h"

/// TODO: make namespaces consistent before function names in log calls.

namespace prism::field {

Scalar::Scalar(std::string name, const SharedPtr<mesh::PMesh>& mesh, f64 value)
    : Scalar(std::move(name), mesh, VectorXd::Ones(mesh->cellCount()) * value, {}, NullOption) {}


Scalar::Scalar(std::string name, const SharedPtr<mesh::PMesh>& mesh, f64 value, VectorCoord coord)
    : Scalar(std::move(name),
             mesh,
             VectorXd::Ones(mesh->cellCount()) * value,
             {},
             Optional<VectorCoord>(coord)) {}


Scalar::Scalar(std::string name, const SharedPtr<mesh::PMesh>& mesh, VectorXd data)
    : Scalar(std::move(name), mesh, std::move(data), {}, NullOption) {}


Scalar::Scalar(std::string name,
               const SharedPtr<mesh::PMesh>& mesh,
               VectorXd values,
               VectorCoord coord)
    : Scalar(std::move(name), mesh, std::move(values), {}, Optional<VectorCoord>(coord)) {}


Scalar::Scalar(std::string name,
               const SharedPtr<mesh::PMesh>& mesh,
               VectorXd values,
               VectorXd face_values)
    : Scalar(std::move(name), mesh, std::move(values), std::move(face_values), NullOption) {}


Scalar::Scalar(std::string name,
               const SharedPtr<mesh::PMesh>& mesh,
               VectorXd values,
               VectorXd face_values,
               VectorCoord coord)
    : Scalar(std::move(name),
             mesh,
             std::move(values),
             std::move(face_values),
             Optional<VectorCoord>(coord)) {}


Scalar::Scalar(std::string name,
               const SharedPtr<mesh::PMesh>& mesh,
               VectorXd cell_values,
               VectorXd face_values,
               Optional<VectorCoord> coord)
    : IScalar(std::move(name), mesh),
      _cell_values(std::move(cell_values)),
      _face_values(std::move(face_values)),
      _history_manager(0),
      _coord(coord) {
    validateCellValues(_cell_values);
    if (_face_values.size() > 0) {
        validateFaceValues(_face_values);
    }
    logConstruction();
    addDefaultBoundaryHandlers();
    setGradScheme();
}


void Scalar::validateCellValues(const VectorXd& values) const {
    if (values.size() != mesh()->cellCount()) {
        throw std::runtime_error(
            fmt::format("field::Scalar() cannot create a scalar field '{}' given a "
                        "vector that has a different size than mesh's cell count.",
                        this->name()));
    }
}


void Scalar::validateFaceValues(const VectorXd& values) const {
    if (values.size() != mesh()->faceCount()) {
        throw std::runtime_error(
            fmt::format("field::Scalar() cannot create a scalar field '{}' given a face data "
                        "vector that has a different size than mesh's faces count.",
                        this->name()));
    }
}


void Scalar::logConstruction() const {
    if (_coord.has_value()) {
        log::debug(
            "Creating scalar field: '{}' (as {}-coordinate) with a cell data vector of "
            "size = {} and face data vector of size = {}",
            this->name(),
            detail::coordToStr(_coord.value()),
            _cell_values.size(),
            _face_values.size());
    } else if (_face_values.size() > 0) {
        log::debug(
            "Creating scalar field: '{}' with a cell data vector of size = {} and face data "
            "vector of size = {}",
            this->name(),
            _cell_values.size(),
            _face_values.size());
    } else {
        log::debug("Creating scalar field: '{}' with a cell vector data of size = {}",
                   this->name(),
                   _cell_values.size());
    }
}

/// TODO: delete this in favor of updateFaceValues()
void Scalar::setFaceValues(VectorXd values) {
    if (values.size() != mesh()->faceCount()) {
        throw std::runtime_error(
            fmt::format("prism::field::Scalar<Units, BHManagerProvider, "
                        "BHManagerSetter>::setFaceValues(): cannot set face values for scalar "
                        "field {}, to a face data vector having a different size that "
                        "field's faces count.",
                        name()));
    }

    if (hasFaceValues()) {
        log::debug(
            "Scalar::setFaceValues() was called for field '{}', which already has face values set.",
            name());
    }

    log::debug("Setting face values for field '{}'", name());
    _face_values = std::move(values);
}

void Scalar::clearFaceValues() {
    if (_face_values.size() == 0) {
        log::warn(
            "Scalar::clearFaceValues() was called for field '{}', but the face data is not "
            "initialized.",
            name());
        return;
    }
    _face_values.resize(0);
}

auto Scalar::values() const -> const VectorXd& {
    assert(_cell_values.size() > 0);
    return _cell_values;
}

auto Scalar::values() -> VectorXd& {
    assert(_cell_values.size() > 0);
    return _cell_values;
}

auto Scalar::coord() const noexcept -> Optional<VectorCoord> {
    return _coord;
}

auto Scalar::clone() const -> Scalar {
    if (hasFaceValues()) {
        return Scalar(name(), mesh(), VectorXd(_cell_values), VectorXd(_face_values));
    }
    return Scalar(name(), mesh(), VectorXd(_cell_values));
}

auto Scalar::hasFaceValues() const -> bool {
    return _face_values.size() > 0;
}

auto Scalar::valueAtCell(const mesh::Cell& cell) const -> f64 {
    return valueAtCell(cell.id());
}

auto Scalar::valueAtCell(std::size_t cell_id) const -> f64 {
    assert(_cell_values.size() > 0);       // NOLINT
    assert(cell_id < mesh()->cellCount()); // NOLINT
    return _cell_values[cell_id];
}

auto Scalar::valueAtFace(std::size_t face_id) const -> f64 {
    if (hasFaceValues()) {
        // Face data were calculated already, just return the value (as in Rhie-Chow
        // correction).
        return _face_values[face_id];
    }

    const auto& face = mesh()->face(face_id);

    if (face.isInterior()) {
        return valueAtInteriorFace(face);
    }

    return valueAtBoundaryFace(face);
}

auto Scalar::valueAtFace(const mesh::Face& face) const -> f64 {
    return valueAtFace(face.id());
}

auto Scalar::valueAtInteriorFace(const mesh::Face& face) const -> f64 {
    assert(face.isInterior()); // NOLINT
    const auto& owner = mesh()->cell(face.owner());
    const auto& neighbor = mesh()->cell(face.neighbor().value());

    const auto gc = mesh::geometricWeight(owner, neighbor, face);
    f64 val = gc * _cell_values[owner.id()];
    val += (1 - gc) * _cell_values[neighbor.id()];

    return val;
}

auto Scalar::valueAtBoundaryFace(const mesh::Face& face) const -> f64 {
    const auto& patch = mesh()->boundaryPatch(face);
    const auto& bc = patch.getBoundaryCondition(name());

    auto handler = this->boundaryHandlersManager().getHandler(bc.kindString());

    if (handler == nullptr) {
        throw error::NonImplementedBoundaryCondition(
            fmt::format("prism::field::Scalar<Units, BHManagerProvider, "
                        "BHManagerSetter>::valueAtBoundaryFace() for field `{}`",
                        name()),
            patch.name(),
            bc.kindString());
    }

    return handler->get(*this, face);
}

void Scalar::addDefaultBoundaryHandlers() {
    log::debug(
        "prism::field::Scalar::set(): adding default boundary handlers for a "
        "scalar field instance");
    this->boundaryHandlersManager().addHandler<field::boundary::scalar::Fixed<Scalar>>();
    this->boundaryHandlersManager().addHandler<field::boundary::scalar::NoSlip<Scalar>>();
    this->boundaryHandlersManager().addHandler<field::boundary::scalar::Symmetry<Scalar>>();
    this->boundaryHandlersManager().addHandler<field::boundary::scalar::Outlet<Scalar>>();
    this->boundaryHandlersManager().addHandler<field::boundary::scalar::ZeroGradient<Scalar>>();
}

void Scalar::setUnits() {}

void Scalar::setGradScheme(const SharedPtr<gradient::IGradient>& grad_scheme) {
    if (grad_scheme == nullptr) {
        throw std::runtime_error(
            "Scalar::setGradScheme() was given a null gradient scheme pointer");
    }
    _grad_scheme = grad_scheme;
}

void Scalar::setGradScheme() {
    // Vector field components resolve the gradient scheme from their parent field, so a
    // `gradScheme` declared for e.g. `U` in `fields.json` applies to `U_x`, `U_y` and `U_z`.
    auto lookup_name = name();
    if (_coord.has_value()) {
        // Component names are `{parent}_x`, `{parent}_y` or `{parent}_z` (enforced by
        // GeneralVector), so stripping the trailing two characters yields the parent name.
        lookup_name = name().substr(0, name().size() - 2);
    }

    // did user specify gradient scheme for the field in `fields.json`?
    auto field_infos = this->mesh()->fieldsInfo();
    auto it = std::find_if(field_infos.begin(), field_infos.end(), [&lookup_name](const auto& fi) {
        return fi.name() == lookup_name && fi.gradScheme().has_value();
    });

    if (it == field_infos.end()) {
        log::debug(
            "Scalar::setGradScheme(): couldn't find a specified gradient scheme for "
            "field `{}` in `fields.json`, setting the gradient scheme to least squares.",
            this->name());

        _grad_scheme = makeShared<gradient::LeastSquares>(this->mesh());
        return;
    }

    auto grad_scheme_name = it->gradScheme().value();

    if (grad_scheme_name == "green-gauss" || grad_scheme_name == "greenGauss") {
        log::debug(
            "Scalar::setGradScheme(): setting the gradient scheme to Green-Gauss for "
            "field `{}`",
            this->name());
        _grad_scheme = makeShared<gradient::GreenGauss>(this->mesh());
        return;
    }


    log::debug(
        "Scalar::setGradScheme(): setting the gradient scheme to Least-Squares "
        "for field `{}`",
        this->name());

    _grad_scheme = makeShared<gradient::LeastSquares>(this->mesh());
}

void Scalar::setHistorySize(std::size_t num_time_steps) {
    _history_manager.resize(num_time_steps);
}


void Scalar::update(VectorXd values) {
    updatePrevTimeSteps();
    _cell_values = std::move(values);

    // clear face values if they exist
    if (hasFaceValues()) {
        clearFaceValues();
    }
}

void Scalar::updatePrevTimeSteps() {
    _history_manager.update(_cell_values);
}

auto Scalar::prevValues() const -> Optional<Ref<const VectorXd>> {
    return _history_manager.prevValues();
}


auto Scalar::prevPrevValues() const -> Optional<Ref<const VectorXd>> {
    return _history_manager.prevPrevValues();
}


auto Scalar::getHistory(std::size_t index) const -> Optional<Ref<const VectorXd>> {
    return _history_manager.valuesAt(index);
}


auto Scalar::gradAtFace(const mesh::Face& face) -> Vector3d {
    return _grad_scheme->gradAtFace(face, *this);
}


auto Scalar::gradAtCell(const mesh::Cell& cell) -> Vector3d {
    return _grad_scheme->gradAtCell(cell, *this);
}


auto Scalar::gradAtCellStored(const mesh::Cell& cell) const -> Vector3d {
    return _grad_scheme->gradAtCellStored(cell, *this);
}

void Scalar::computeAllCellGradients() {
    for (const auto& cell : this->mesh()->cells()) {
        this->gradAtCell(cell);
    }
}


auto Scalar::operator[](std::size_t i) const -> f64 {
    assert(_cell_values.size() > 0);
    assert(i < _cell_values.size());
    return _cell_values[i];
}


auto Scalar::operator[](std::size_t i) -> f64& {
    assert(_cell_values.size() > 0);
    assert(i < _cell_values.size());
    return _cell_values[i];
}

} // namespace prism::field
