
#include "field.h"

#include <fmt/core.h>

#include <stdexcept>

#include "exceptions.h"
#include "prism/mesh/boundary.h"
#include "prism/mesh/cell.h"
#include "prism/mesh/face.h"
#include "prism/mesh/pmesh.h"
#include "prism/mesh/utilities.h"

namespace prism {

void inline check_field_name(const std::string& name) {
    if (name.empty()) {
        throw std::runtime_error("Cannot create a Field with an empty name.");
    }
}

void inline check_mesh(const mesh::PMesh& mesh) {
    if (mesh.cells().empty() || mesh.faces().empty() || mesh.boundary_patches().empty()) {
        throw std::runtime_error("Cannot create a field over an empty mesh.");
    }
}

ScalarField::ScalarField(std::string name, const mesh::PMesh& mesh, double value)
    : _mesh(&mesh),
      _name(std::move(name)),
      _data(std::make_shared<VectorXd>(VectorXd::Ones(mesh.n_cells()) * value)) {
    check_field_name(_name);
    check_mesh(mesh);
}

ScalarField::ScalarField(std::string name, const mesh::PMesh& mesh, VectorXd data)
    : _mesh(&mesh), _name(std::move(name)), _data(std::make_shared<VectorXd>(std::move(data))) {
    check_field_name(_name);
    check_mesh(mesh);
}

ScalarField::ScalarField(std::string name,
                         const mesh::PMesh& mesh,
                         VectorXd data,
                         VectorXd face_data)
    : _mesh(&mesh),
      _name(std::move(name)),
      _data(std::make_shared<VectorXd>(std::move(data))),
      _face_data(std::make_shared<VectorXd>(std::move(face_data))) {
    check_field_name(_name);
    check_mesh(mesh);
}

auto ScalarField::clone() const -> ScalarField {
    if (has_face_data()) {
        return {_name, *_mesh, *_data, *_face_data};
    }
    return {_name, *_mesh, *_data};
}

auto ScalarField::map(CellMapper* mapper) -> ScalarField& {
    for (std::size_t i = 0; i < _mesh->n_cells(); ++i) {
        data()[i] = mapper(_mesh->cell(i));
    }
    return *this;
}

auto ScalarField::map(CoordinatesMapper* mapper) -> ScalarField& {
    const auto n_cells = _mesh->n_cells();
    for (std::size_t i = 0; i < n_cells; ++i) {
        const auto& cell = _mesh->cell(i);
        const auto& center = cell.center();
        data()[i] = mapper(center.x(), center.y(), center.z());
    }
    return *this;
}

auto ScalarField::value_at_cell(std::size_t cell_id) const -> double {
    return (*_data)[cell_id];
}

auto ScalarField::value_at_cell(const mesh::Cell& cell) const -> double {
    return value_at_cell(cell.id());
}

auto ScalarField::value_at_face(std::size_t face_id) const -> double {
    if (has_face_data()) {
        // Face data were calculataed for us before calling the constructor,
        // just return the value
        return (*_face_data)[face_id];
    }

    // We need to interpolate the value of the field at the face
    const auto& face = _mesh->face(face_id);

    if (face.is_interior()) {
        return value_at_interior_face(face);
    }

    return value_at_boundary_face(face);
}

auto ScalarField::value_at_face(const mesh::Face& face) const -> double {
    return value_at_face(face.id());
}

auto ScalarField::value_at_interior_face(const mesh::Face& face) const -> double {
    const auto& owner = _mesh->cell(face.owner());
    const auto& neighbor = _mesh->cell(face.neighbor().value());
    const auto gc = mesh::geo_weight(owner, neighbor, face);

    double val = gc * (*_data)[owner.id()];
    val += (1 - gc) * (*_data)[neighbor.id()];

    return val;
}

auto ScalarField::value_at_boundary_face(const mesh::Face& face) const -> double {
    const auto& patch = _mesh->boundary_patch(face);
    const auto& bc = patch.get_bc(_name);

    switch (bc.bc_type()) {
        case mesh::BoundaryConditionType::Fixed:
        case mesh::BoundaryConditionType::Inlet: {
            return patch.get_scalar_bc(_name);
        }

        // TODO: We return the field value of an empty face the same value as its owner cell.
        // This makes many schemes and gradient methods to work without a special check for an
        // empty face (like LeastSquares), check if this assumption is correct.
        case mesh::BoundaryConditionType::Empty:
        case mesh::BoundaryConditionType::Symmetry:
        case mesh::BoundaryConditionType::Outlet: {
            return (*_data)[face.owner()];
        }

        case mesh::BoundaryConditionType::FixedGradient: {
            Vector3d grad_at_boundary = patch.get_vector_bc(_name);
            const auto& owner = _mesh->cell(face.owner());
            Vector3d e = face.center() - owner.center();
            double d_Cf = e.norm();
            e = e / e.norm();
            grad_at_boundary = grad_at_boundary * d_Cf;

            return grad_at_boundary.dot(e) + value_at_cell(owner);
        }

        default: {
            throw error::NonImplementedBoundaryCondition(
                fmt::format("ScalarField({})::value_at_boundary_face()", _name),
                patch.name(),
                bc.bc_type_str());
        }
    }
}

VectorField::VectorField(std::string name, const mesh::PMesh& mesh, double value)
    : _mesh(&mesh),
      _name(std::move(name)),
      _x(_name + "_x", mesh, value),
      _y(_name + "_y", mesh, value),
      _z(_name + "_z", mesh, value) {
    check_field_name(_name);
    check_mesh(mesh);
}

VectorField::VectorField(std::string name, const mesh::PMesh& mesh, const Vector3d& data)
    : _mesh(&mesh),
      _name(std::move(name)),
      _x(_name + "_x", mesh, data[0]),
      _y(_name + "_y", mesh, data[1]),
      _z(_name + "_z", mesh, data[2]) {
    check_field_name(_name);
    check_mesh(mesh);
}

VectorField::VectorField(std::string name,
                         const mesh::PMesh& mesh,
                         const std::array<ScalarField, 3>& fields)
    : _mesh(&mesh), _name(std::move(name)), _x(fields[0]), _y(fields[1]), _z(fields[2]) {
    check_field_name(_name);
    check_mesh(mesh);

    // check mesh consistency
    for (const auto& field : fields) {
        if (&mesh != &field.mesh()) {
            throw std::runtime_error(
                fmt::format("VectorField constructor was given a ScalarField component with name "
                            "`{}` that is defined over a different mesh",
                            field.name()));
        }
    }

    // check sub-fields naming consistency
    if ((_x.name() != (_name + "_x")) || (_y.name() != (_name + "_y")) ||
        (_z.name() != (_name + "_z"))) {
        throw std::runtime_error(fmt::format(
            "All VectorField component names should end with '_x', '_y' or '_z'. VectorField "
            "constructor for `{}` vector field was given the following ScalarFields names: '{}', "
            "'{}', '{}",
            _name,
            _x.name(),
            _y.name(),
            _z.name()));
    }
}

auto VectorField::value_at_cell(std::size_t cell_id) const -> Vector3d {
    return operator[](cell_id);
}

auto VectorField::value_at_cell(const mesh::Cell& cell) const -> Vector3d {
    return value_at_cell(cell.id());
}

auto VectorField::value_at_face(std::size_t face_id) const -> Vector3d {
    return {_x.value_at_face(face_id), _y.value_at_face(face_id), _z.value_at_face(face_id)};
}

auto VectorField::value_at_face(const mesh::Face& face) const -> Vector3d {
    return value_at_face(face.id());
}

auto VectorField::has_face_data() const -> bool {
    return _x.has_face_data() && _y.has_face_data() && _z.has_face_data();
}

} // namespace prism