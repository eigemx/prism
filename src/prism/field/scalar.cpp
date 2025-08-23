#include "scalar.h"

#include <algorithm>

#include "prism/exceptions.h"
#include "prism/gradient/green_gauss.h"
#include "prism/gradient/least_squares.h"
#include "prism/log.h"
#include "prism/mesh/utilities.h"

/// TODO: make namespaces consistent before function names in log calls.

namespace prism::field {
Scalar::Scalar(std::string name, const SharedPtr<mesh::PMesh>& mesh, f64 value)
    : IScalar(std::move(name), mesh),
      _history_manager(0),
      _cell_values(VectorXd::Ones(mesh->cellCount()) * value) {
    log::debug("Creating scalar field: '{}' with f64 value = {}", this->name(), value);
    addDefaultBoundaryHandlers();
    setGradScheme();
}


Scalar::Scalar(std::string name, const SharedPtr<mesh::PMesh>& mesh, f64 value, Coord coord)
    : IScalar(std::move(name), mesh),
      _cell_values(VectorXd::Ones(mesh->cellCount()) * value),
      _history_manager(0),
      _coord(coord) {
    log::debug("Creating scalar field: '{}' (as {}-coordinate) with f64 value = {}",
               this->name(),
               coordToStr(coord),
               value);

    addDefaultBoundaryHandlers();
    setGradScheme();
}


Scalar::Scalar(std::string name, const SharedPtr<mesh::PMesh>& mesh, VectorXd data)
    : IScalar(std::move(name), mesh), _cell_values(std::move(data)), _history_manager(0) {
    if (_cell_values.size() != mesh->cellCount()) {
        throw std::runtime_error(
            fmt::format("field::Scalar() cannot create a scalar field '{}' given a "
                        "vector that has a different size than mesh's cell count.",
                        this->name()));
    }

    log::debug("Creating scalar field: '{}' with a cell vector data of size = {}",
               this->name(),
               _cell_values.size());

    addDefaultBoundaryHandlers();
    setGradScheme();
}


Scalar::Scalar(std::string name, const SharedPtr<mesh::PMesh>& mesh, VectorXd data, Coord coord)
    : IScalar(std::move(name), mesh),
      _cell_values(std::move(data)),
      _history_manager(0),
      _coord(coord) {
    if (_cell_values.size() != mesh->cellCount()) {
        throw std::runtime_error(fmt::format(
            "field::Scalar() cannot create a scalar field '{}' given a vector that has a "
            "different size than mesh's cell count.",
            this->name()));
    }

    log::debug(
        "Creating scalar field: '{}' (as {}-coordinate) with a cell vector data of size = {}",
        this->name(),
        coordToStr(coord),
        _cell_values.size());
    addDefaultBoundaryHandlers();
    setGradScheme();
}


Scalar::Scalar(std::string name,
               const SharedPtr<mesh::PMesh>& mesh,
               VectorXd data,
               VectorXd face_data)
    : IScalar(std::move(name), mesh),
      _cell_values(std::move(data)),
      _face_values(std::move(face_data)),
      _history_manager(0) {
    if (_cell_values.size() != mesh->cellCount()) {
        throw std::runtime_error(fmt::format(
            "field::Scalar() cannot create a scalar field '{}' given a vector that has a "
            "different size than mesh's cell count.",
            this->name()));
    }

    if (_face_values.size() != mesh->faceCount()) {
        throw std::runtime_error(
            fmt::format("field::Scalar() cannot create a scalar field '{}' given a face data "
                        "vector that has a different size than mesh's faces count.",
                        this->name()));
    }

    log::debug(
        "Creating scalar field: '{}' with a cell data vector of size = {} and face data "
        "vector of size = {}",
        this->name(),
        _cell_values.size(),
        _face_values.size());

    addDefaultBoundaryHandlers();
    setGradScheme();
}


Scalar::Scalar(std::string name,
               const SharedPtr<mesh::PMesh>& mesh,
               VectorXd data,
               VectorXd face_data,
               Coord coord)
    : IScalar(std::move(name), mesh),
      _cell_values(std::move(data)),
      _face_values(std::move(face_data)),
      _history_manager(0),
      _coord(coord) {
    if (_cell_values.size() != mesh->cellCount()) {
        throw std::runtime_error(fmt::format(
            "field::Scalar() cannot create a scalar field '{}' given a vector that has a "
            "different size than mesh's cell count.",
            this->name()));
    }

    if (_face_values.size() != mesh->faceCount()) {
        throw std::runtime_error(
            fmt::format("field::Scalar() cannot create a scalar field '{}' given a face data "
                        "vector that has a different size than mesh's faces count.",
                        this->name()));
    }

    log::debug(
        "Creating scalar field: '{}' (as {}-coordinate) with a cell data vector of size = {} "
        "and face data vector of size = {}",
        this->name(),
        coordToStr(coord),
        _cell_values.size(),
        _face_values.size());

    addDefaultBoundaryHandlers();
    setGradScheme();
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
            "Scalar::setFaceValues() was called for field "
            "'{}', which already has face values set.",
            name());
    }

    log::debug("Setting face values for field '{}'", name());
    _face_values = std::move(values);
}

void Scalar::clearFaceValues() {
    if (_face_values.size() > 0) {
        log::warn(
            "Scalar::clearFaceValues() was called for field '{}', but the face data is not "
            "initialized.",
            name());
        return;
    }
    _face_values.resize(0);
}

auto Scalar::values() const -> const VectorXd& {
    if (_cell_values.size() == 0) {
        throw std::runtime_error(
            fmt::format("prism::field::Scalar::values() was called for field `{}`, but the data "
                        "is not initialized.",
                        name()));
    }
    return _cell_values;
}

auto Scalar::values() -> VectorXd& {
    if (_cell_values.size() == 0) {
        throw std::runtime_error(
            fmt::format("prism::field::Scalar::values() was "
                        "called for field `{}`, but the data is not initialized.",
                        name()));
    }
    return _cell_values;
}

auto Scalar::coord() const noexcept -> Optional<Coord> {
    return _coord;
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
    // if this scalar field is a component of a vector field, we need to return the value of the
    // parent field at the face in the coordinate direction.
    /// TODO: implement this using make_shared_from_this and std::weak_ptr
    if (_coord.has_value()) {
        auto idx = detail::coordToIndex(_coord.value());
        // return _parent->valueAtFace(face)[idx];
    }
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
    // did user specify gradient scheme for the field in `fields.json`?
    auto field_infos = this->mesh()->fieldsInfo();
    auto it = std::find_if(field_infos.begin(), field_infos.end(), [this](const auto& fi) {
        return fi.name() == this->name() && fi.gradScheme().has_value();
    });

    /// TODO: this is buggy, it doesn't find the grad scheme defined in fields.json for the
    // field, also does not consider vector fields.
    if (it == field_infos.end()) {
        log::debug(
            "Scalar::setGradScheme(): couldn't find a specified gradient scheme for "
            "field `{}` in `fields.json`, setting the gradient scheme to least squares.",
            this->name());

        _grad_scheme = std::make_shared<gradient::LeastSquares>(this->mesh());
        return;
    }

    auto grad_scheme_name = it->gradScheme().value();

    if (grad_scheme_name == "green-gauss" || grad_scheme_name == "greenGauss") {
        log::debug(
            "Scalar::setGradScheme(): setting the gradient scheme to Green-Gauss for "
            "field `{}`",
            this->name());
        _grad_scheme = std::make_shared<gradient::GreenGauss>(this->mesh());
        return;
    }

    /// TODO: if scalar field is a component of a vector field, prism should get the grad scheme
    /// type from it's parent field, and we should also implement setGradScheme(scheme) for vector
    /// field type.

    log::debug(
        "Scalar::setGradScheme(): setting the gradient scheme to Least-Squares "
        "for field `{}`",
        this->name());

    _grad_scheme = std::make_shared<gradient::LeastSquares>(this->mesh());
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

auto Scalar::prevValues() const -> Optional<VectorXd> {
    /// TODO: implement this
    return NullOption;
}


auto Scalar::prevPrevValues() const -> Optional<VectorXd> {
    /// TODO: implement this
    return NullOption;
}


auto Scalar::getHistory(std::size_t index) const -> Optional<VectorXd> {
    /// TODO: implement this
    return NullOption;
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


template <typename Func>
void Scalar::updateInteriorFaces(Func func) {
    if (!hasFaceValues()) {
        /// TODO: implement initFaceDataVector() to initialize face data vector
        VectorXd face_values(mesh()->faceCount());
        face_values.setZero();

        for (const auto& patch : mesh()->boundaryPatches()) {
            if (patch.isEmpty()) {
                continue; // skip empty patches
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


template <typename Func>
void Scalar::updateFaces(Func func) {
    updateInteriorFaces(func);

    for (const auto& patch : mesh()->boundaryPatches()) {
        if (patch.isEmpty()) {
            continue; // skip empty patches
        }
        for (const auto& face_id : patch.facesIds()) {
            const auto& face = mesh()->face(face_id);
            _face_values[face_id] = func(face);
        }
    }
}


template <typename Func>
void Scalar::updateCells(Func func) {
    /// TODO: test this.
    for (const auto& cell : mesh()->cells()) {
        _cell_values[cell.id()] = func(cell);
    }
}

auto Scalar::operator[](std::size_t i) const -> f64 {
    if (_cell_values.size() == 0) {
        throw std::runtime_error(
            fmt::format("prism::field::Scalar::operator[]() was called for field `{}`, but the "
                        "data is not initialized.",
                        name()));
    }
    if (i >= _cell_values.size()) {
        throw std::out_of_range(
            fmt::format("prism::field::Scalar::operator[]() was called for field `{}`, but the "
                        "index {} is out of range (size = {}).",
                        name(),
                        i,
                        _cell_values.size()));
    }
    return _cell_values[i];
}


auto Scalar::operator[](std::size_t i) -> f64& {
    if (_cell_values.size() == 0) {
        throw std::runtime_error(
            fmt::format("prism::field::Scalar::operator[]() was called for field `{}`, but the "
                        "data is not initialized.",
                        name()));
    }

    if (i >= _cell_values.size()) {
        throw std::out_of_range(
            fmt::format("prism::field::Scalar::operator[]() was called for field `{}`, but the "
                        "index {} is out of range (size = {}).",
                        name(),
                        i,
                        _cell_values.size()));
    }
    return _cell_values[i];
}

} // namespace prism::field
