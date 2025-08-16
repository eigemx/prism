#include <string>

#include "prism/boundary.h"
#include "prism/field/ifield.h"
#include "prism/mesh/face.h"
#include "prism/types.h"

namespace prism::gradient::boundary {
class IGradSchemeBoundaryHandler : public prism::boundary::IBoundaryHandler {
  public:
    virtual auto name() const -> std::string = 0;
    virtual auto get(field::IScalar& field, const mesh::Face& face) -> Vector3d = 0;
};

class Symmetry : public IGradSchemeBoundaryHandler {
  public:
    auto name() const -> std::string override { return "symmetry"; }
    auto get(field::IScalar& field, const mesh::Face& face) -> Vector3d override;
};

class Outlet : public IGradSchemeBoundaryHandler {
  public:
    auto name() const -> std::string override { return "outlet"; }
    auto get(field::IScalar& field, const mesh::Face& face) -> Vector3d override;
};

class Fixed : public IGradSchemeBoundaryHandler {
  public:
    auto name() const -> std::string override { return "fixed"; }
    auto get(field::IScalar& field, const mesh::Face& face) -> Vector3d override;
};

class NoSlip : public IGradSchemeBoundaryHandler {
  public:
    auto name() const -> std::string override { return "no-slip"; }
    auto get(field::IScalar& field, const mesh::Face& face) -> Vector3d override;
};

class ZeroGradient : public IGradSchemeBoundaryHandler {
  public:
    auto name() const -> std::string override { return "zero-gradient"; }
    auto get(field::IScalar& field, const mesh::Face& face) -> Vector3d override;
};
} // namespace prism::gradient::boundary
