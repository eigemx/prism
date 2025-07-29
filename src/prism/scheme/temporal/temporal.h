#include "prism/scheme/scheme.h"

namespace prism::scheme::temporal {

template <typename Field>
class ITemporal : public scheme::IFullScheme<Field> {
  public:
    ITemporal(std::size_t n_cells) : scheme::IFullScheme<Field>(n_cells) {}
    virtual void apply() = 0;

  private:
    void applyInterior(const mesh::Face& face) override {}
    void applyBoundary() override {}
};

} // namespace prism::scheme::temporal
