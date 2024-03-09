#pragma once

#include <string>

#include "prism/field.h"
#include "prism/mesh/boundary.h"

namespace prism::boundary {

template <typename Scheme, typename Field>
class BoundaryCondition {
  public:
    BoundaryCondition() = default;
    BoundaryCondition(BoundaryCondition&) = default;
    BoundaryCondition(BoundaryCondition&&) noexcept = default;
    auto operator=(const BoundaryCondition&) -> BoundaryCondition& = default;
    auto operator=(BoundaryCondition&&) noexcept -> BoundaryCondition& = default;
    virtual ~BoundaryCondition() = default;

    virtual void apply(Scheme& scheme,
                       const Field& field,
                       const mesh::BoundaryPatch& patch) const = 0;

    virtual auto name() const -> std::string = 0;
};

template <typename Scheme, typename Field>
class Fixed : public BoundaryCondition<Scheme, Field> {
  public:
    void apply(Scheme& scheme,
               const Field& field,
               const mesh::BoundaryPatch& patch) const override = 0;
    auto inline name() const -> std::string override { return "fixed"; }
};

template <typename Scheme, typename Field>
class VelocityInlet : public BoundaryCondition<Scheme, Field> {
  public:
    void apply(Scheme& scheme,
               const Field& field,
               const mesh::BoundaryPatch& patch) const override = 0;
    auto inline name() const -> std::string override { return "velocity-intlet"; }
};

template <typename Scheme, typename Field>
class Outlet : public BoundaryCondition<Scheme, Field> {
  public:
    void apply(Scheme& scheme,
               const Field& field,
               const mesh::BoundaryPatch& patch) const override = 0;
    auto inline name() const -> std::string override { return "outlet"; }
};

template <typename Scheme, typename Field>
class Symmetry : public BoundaryCondition<Scheme, Field> {
  public:
    void apply(Scheme& scheme,
               const Field& field,
               const mesh::BoundaryPatch& patch) const override = 0;
    auto inline name() const -> std::string override { return "symmetry"; }
};

template <typename Scheme, typename Field>
class FixedGradient : public BoundaryCondition<Scheme, Field> {
  public:
    void apply(Scheme& scheme,
               const Field& field,
               const mesh::BoundaryPatch& patch) const override = 0;
    auto inline name() const -> std::string override { return "fixed-gradient"; }
};

template <typename Scheme, typename Field>
class SlipWall : public BoundaryCondition<Scheme, Field> {
  public:
    void apply(Scheme& scheme,
               const Field& field,
               const mesh::BoundaryPatch& patch) const override = 0;
    auto inline name() const -> std::string override { return "slip-wall"; }
};

template <typename Scheme, typename Field>
class NonSlipWall : public BoundaryCondition<Scheme, Field> {
  public:
    void apply(Scheme& scheme,
               const Field& field,
               const mesh::BoundaryPatch& patch) const override = 0;
    auto inline name() const -> std::string override { return "nonslip-wall"; }
};

} // namespace prism::boundary