#pragma once

#include <algorithm>

#include "ifield.h"
#include "prism/exceptions.h"
#include "prism/gradient/gradient.h"
#include "prism/log.h"
#include "prism/mesh/utilities.h"
#include "scalar_boundary.h"
#include "units.h"

namespace prism::field {

class UniformScalar : public IScalar {
  public:
    UniformScalar(std::string name, const mesh::PMesh& mesh, double value);

    auto valueAtCell(std::size_t cell_id) const -> double override;
    auto valueAtCell(const mesh::Cell& cell) const -> double override;
    auto valueAtFace(std::size_t face_id) const -> double override;
    auto valueAtFace(const mesh::Face& face) const -> double override;

    auto gradAtFace(const mesh::Face& face) const -> Vector3d override;
    auto gradAtCell(const mesh::Cell& cell) const -> Vector3d override;
    auto gradAtCellStored(const mesh::Cell& cell) const -> Vector3d override;

  private:
    double _value {0.0};
    SharedPtr<gradient::IGradient> _grad_scheme = nullptr;
};

template <typename Units, typename BHManagerSetter>
class GeneralScalar
    : public IScalar,
      public Units,
      public prism::boundary::BHManagersProvider<boundary::IScalarBoundaryHandler> {
  public:
    // Uniform double value constructors
    GeneralScalar(std::string name,
                  const mesh::PMesh& mesh,
                  double value,
                  IVector* parent = nullptr);
    GeneralScalar(std::string name,
                  const mesh::PMesh& mesh,
                  double value,
                  Coord coord,
                  IVector* parent = nullptr);

    // VectorXd cell values constructors
    GeneralScalar(std::string name,
                  const mesh::PMesh& mesh,
                  VectorXd data,
                  IVector* parent = nullptr);
    GeneralScalar(std::string name,
                  const mesh::PMesh& mesh,
                  VectorXd data,
                  Coord coord,
                  IVector* parent = nullptr);

    // VectorXd cell & face values constructors
    GeneralScalar(std::string name,
                  const mesh::PMesh& mesh,
                  VectorXd data,
                  VectorXd face_data,
                  IVector* parent = nullptr);
    GeneralScalar(std::string name,
                  const mesh::PMesh& mesh,
                  VectorXd data,
                  VectorXd face_data,
                  Coord coord,
                  IVector* parent = nullptr);

    GeneralScalar() = delete;
    GeneralScalar(const GeneralScalar&) = default;
    GeneralScalar(GeneralScalar&&) = default;
    auto operator=(const GeneralScalar&) -> GeneralScalar& = default;
    auto operator=(GeneralScalar&&) -> GeneralScalar& = default;
    ~GeneralScalar() override = default;

    // TODO: check that _data is not null before returning, and maybe wrap it in an optional type
    auto inline values() const -> const VectorXd& { return *_data; }
    auto inline values() -> VectorXd& { return *_data; }

    auto inline coord() const noexcept -> std::optional<Coord> override { return _coord; }
    auto inline hasFaceValues() const -> bool override { return _face_data != nullptr; }
    void setFaceValues(VectorXd values);

    auto valueAtCell(std::size_t cell_id) const -> double override;
    auto valueAtCell(const mesh::Cell& cell) const -> double override;
    auto valueAtFace(std::size_t face_id) const -> double override;
    auto valueAtFace(const mesh::Face& face) const -> double override;

    auto parent() -> IVector*;
    void setParent(IVector* parent);

    auto gradAtFace(const mesh::Face& face) const -> Vector3d override;
    auto gradAtCell(const mesh::Cell& cell) const -> Vector3d override;
    auto gradAtCellStored(const mesh::Cell& cell) const -> Vector3d override;

    void setGradScheme(const SharedPtr<gradient::IGradient>& grad_scheme);

    auto inline operator[](std::size_t i) const -> double { return (*_data)[i]; }
    auto inline operator[](std::size_t i) -> double& { return (*_data)[i]; }

  protected:
    auto valueAtInteriorFace(const mesh::Face& face) const -> double;
    auto valueAtBoundaryFace(const mesh::Face& face) const -> double;

  private:
    void setGradScheme();
    void addDefaultBoundaryHandlers();

    SharedPtr<VectorXd> _data = nullptr;

    // TODO: _face_data should not include empty faces
    SharedPtr<VectorXd> _face_data = nullptr;
    SharedPtr<gradient::IGradient> _grad_scheme;

    IVector* _parent = nullptr;
    std::optional<Coord> _coord = std::nullopt;
    BHManagerSetter _setter;
};

class ScalarBHManagerSetter {
  public:
    using IScalarBHManager =
        prism::boundary::BoundaryHandlersManager<boundary::IScalarBoundaryHandler>;

    static void set(IScalarBHManager& manager);
};

using Scalar = GeneralScalar<units::Measurable, ScalarBHManagerSetter>;

template <typename Units, typename BHManagerSetter>
GeneralScalar<Units, BHManagerSetter>::GeneralScalar(std::string name,
                                                     const mesh::PMesh& mesh,
                                                     double value,
                                                     IVector* parent)
    : IScalar(std::move(name), mesh),
      _data(std::make_shared<VectorXd>(VectorXd::Ones(mesh.cellCount()) * value)),
      _parent(parent) {
    log::debug("Creating scalar field: '{}' with double value = {}", this->name(), value);
    addDefaultBoundaryHandlers();
    setGradScheme();
}

template <typename Units, typename BHManagerSetter>
GeneralScalar<Units, BHManagerSetter>::GeneralScalar(std::string name,
                                                     const mesh::PMesh& mesh,
                                                     double value,
                                                     Coord coord,
                                                     IVector* parent)
    : IScalar(std::move(name), mesh),
      _data(std::make_shared<VectorXd>(VectorXd::Ones(mesh.cellCount()) * value)),
      _coord(coord),
      _parent(parent) {
    log::debug("Creating scalar field: '{}' (as {}-coordinate) with double value = {}",
               this->name(),
               coordToStr(coord),
               value);
    addDefaultBoundaryHandlers();
    setGradScheme();
}

template <typename Units, typename BHManagerSetter>
GeneralScalar<Units, BHManagerSetter>::GeneralScalar(std::string name,
                                                     const mesh::PMesh& mesh,
                                                     VectorXd data,
                                                     IVector* parent)
    : IScalar(std::move(name), mesh),
      _data(std::make_shared<VectorXd>(std::move(data))),
      _parent(parent) {
    if (_data->size() != mesh.cellCount()) {
        throw std::runtime_error(fmt::format(
            "field::Scalar() cannot create a scalar field '{}' given a vector that has a "
            "different size than mesh's cell count.",
            this->name()));
    }

    log::debug("Creating scalar field: '{}' with a cell vector data of size = {}",
               this->name(),
               _data->size());
    addDefaultBoundaryHandlers();
    setGradScheme();
}

template <typename Units, typename BHManagerSetter>
GeneralScalar<Units, BHManagerSetter>::GeneralScalar(std::string name,
                                                     const mesh::PMesh& mesh,
                                                     VectorXd data,
                                                     Coord coord,
                                                     IVector* parent)
    : IScalar(std::move(name), mesh),
      _data(std::make_shared<VectorXd>(std::move(data))),
      _coord(coord),
      _parent(parent) {
    if (_data->size() != mesh.cellCount()) {
        throw std::runtime_error(fmt::format(
            "field::Scalar() cannot create a scalar field '{}' given a vector that has a "
            "different size than mesh's cell count.",
            this->name()));
    }

    log::debug(
        "Creating scalar field: '{}' (as {}-coordinate) with a cell vector data of size = {}",
        this->name(),
        coordToStr(coord),
        _data->size());
    addDefaultBoundaryHandlers();
    setGradScheme();
}

template <typename Units, typename BHManagerSetter>
GeneralScalar<Units, BHManagerSetter>::GeneralScalar(std::string name,
                                                     const mesh::PMesh& mesh,
                                                     VectorXd data,
                                                     VectorXd face_data,
                                                     IVector* parent)
    : IScalar(std::move(name), mesh),
      _data(std::make_shared<VectorXd>(std::move(data))),
      _face_data(std::make_shared<VectorXd>(std::move(face_data))),
      _parent(parent) {
    if (_data->size() != mesh.cellCount()) {
        throw std::runtime_error(fmt::format(
            "field::Scalar() cannot create a scalar field '{}' given a vector that has a "
            "different size than mesh's cell count.",
            this->name()));
    }

    if (_face_data->size() != mesh.faceCount()) {
        throw std::runtime_error(
            fmt::format("field::Scalar() cannot create a scalar field '{}' given a face data "
                        "vector that has a different size than mesh's faces count.",
                        this->name()));
    }

    log::debug(
        "Creating scalar field: '{}' with a cell data vector of size = {} and face data vector "
        "of size = {}",
        this->name(),
        _data->size(),
        _face_data->size());

    addDefaultBoundaryHandlers();
    setGradScheme();
}

template <typename Units, typename BHManagerSetter>
GeneralScalar<Units, BHManagerSetter>::GeneralScalar(std::string name,
                                                     const mesh::PMesh& mesh,
                                                     VectorXd data,
                                                     VectorXd face_data,
                                                     Coord coord,
                                                     IVector* parent)
    : IScalar(std::move(name), mesh),
      _data(std::make_shared<VectorXd>(std::move(data))),
      _face_data(std::make_shared<VectorXd>(std::move(face_data))),
      _coord(coord),
      _parent(parent) {
    if (_data->size() != mesh.cellCount()) {
        throw std::runtime_error(fmt::format(
            "field::Scalar() cannot create a scalar field '{}' given a vector that has a "
            "different size than mesh's cell count.",
            this->name()));
    }

    if (_face_data->size() != mesh.faceCount()) {
        throw std::runtime_error(
            fmt::format("field::Scalar() cannot create a scalar field '{}' given a face data "
                        "vector that has a different size than mesh's faces count.",
                        this->name()));
    }

    log::debug(
        "Creating scalar field: '{}' (as {}-coordinate) with a cell data vector of size = {} and "
        "face data vector "
        "of size = {}",
        this->name(),
        coordToStr(coord),
        _data->size(),
        _face_data->size());

    addDefaultBoundaryHandlers();
    setGradScheme();
}

template <typename Units, typename BHManagerSetter>
void GeneralScalar<Units, BHManagerSetter>::setFaceValues(VectorXd values) {
    if (values.size() != mesh().faceCount()) {
        throw std::runtime_error(fmt::format(
            "prism::field::GeneralScalar<Units, BHManagerProvider, "
            "BHManagerSetter>::setFaceValues(): cannot set face values for scalar field {}, to a "
            "face data vector having a different size that field's faces count.",
            name()));
    }

    if (hasFaceValues()) {
        log::debug("Setting new face values to scalar field '{}', discarding old face values.",
                   name());
    }

    _face_data = std::make_shared<VectorXd>(std::move(values));
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::valueAtCell(const mesh::Cell& cell) const -> double {
    return valueAtCell(cell.id());
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::valueAtCell(std::size_t cell_id) const -> double {
    assert(_data != nullptr);             // NOLINT
    assert(cell_id < mesh().cellCount()); // NOLINT
    return (*_data)[cell_id];
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::valueAtFace(std::size_t face_id) const -> double {
    if (hasFaceValues()) {
        // Face data were calculataed for us, just return the value (as in Rhie-Chow corrected
        // face values).
        return (*_face_data)[face_id];
    }

    // We need to interpolate the value of the field at the face
    const auto& face = mesh().face(face_id);

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
    const auto& owner = mesh().cell(face.owner());
    const auto& neighbor = mesh().cell(face.neighbor().value());

    const auto gc = mesh::geometricWeight(owner, neighbor, face);
    double val = gc * (*_data)[owner.id()];
    val += (1 - gc) * (*_data)[neighbor.id()];

    return val;
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::valueAtBoundaryFace(const mesh::Face& face) const
    -> double {
    const auto& patch = mesh().boundaryPatch(face);
    const auto& bc = patch.getBoundaryCondition(name());

    auto handler = this->boundaryHandlersManager().getHandler(bc.kindString());

    if (handler == nullptr) {
        throw error::NonImplementedBoundaryCondition(
            fmt::format("prism::field::GeneralScalar<Units, BHManagerProvider, "
                        "BHManagerSetter>::valueAtBoundaryFace() {}",
                        name()),
            patch.name(),
            bc.kindString());
    }

    return handler->get(*this, face);
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::parent() -> IVector* {
    return _parent;
}

template <typename Units, typename BHManagerSetter>
void GeneralScalar<Units, BHManagerSetter>::setParent(IVector* parent) {
    _parent = parent;
    // TODO: check parent and component names consistency
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
    auto field_infos = this->mesh().fieldsInfo();
    auto it = std::find_if(field_infos.begin(), field_infos.end(), [this](const auto& fi) {
        return fi.name() == this->name() && fi.gradScheme().has_value();
    });

    // TODO: this is buggy, it doesn't find the grad scheme defined in fields.json for the field,
    // also does not consider vector fields.
    if (it == field_infos.end()) {
        log::debug(
            "GeneralScalar::setGradScheme(): couldn't find a specified gradient scheme for field "
            "`{}` in `fields.json`, setting the gradient scheme to least squares.",
            this->name());

        _grad_scheme = std::make_shared<gradient::LeastSquares>(this);
        return;
    }

    auto grad_scheme_name = it->gradScheme().value();

    if (grad_scheme_name == "green-gauss" || grad_scheme_name == "greenGauss") {
        _grad_scheme = std::make_shared<gradient::GreenGauss>(this);
        return;
    }

    if (grad_scheme_name == "least-squares" || grad_scheme_name == "leastSquares") {
        _grad_scheme = std::make_shared<gradient::LeastSquares>(this);
        return;
    }
    _grad_scheme = std::make_shared<gradient::LeastSquares>(this);
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::gradAtFace(const mesh::Face& face) const -> Vector3d {
    return _grad_scheme->gradAtFace(face);
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::gradAtCell(const mesh::Cell& cell) const -> Vector3d {
    return _grad_scheme->gradAtCell(cell);
}

template <typename Units, typename BHManagerSetter>
auto GeneralScalar<Units, BHManagerSetter>::gradAtCellStored(const mesh::Cell& cell) const
    -> Vector3d {
    return _grad_scheme->gradAtCellStored(cell);
}

} // namespace prism::field
