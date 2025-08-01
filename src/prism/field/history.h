#pragma once

#include <vector>

#include "prism/types.h"

namespace prism::field {

class HistoryManager {
  public:
    HistoryManager(int num_time_steps);

    auto prevValues() const -> Optional<VectorXd>;
    auto prevPrevValues() const -> Optional<VectorXd>;
    auto valuesAt(std::size_t n) const -> Optional<VectorXd>;
    void update(const VectorXd& current_values);
    void update(VectorXd&& current_values);
    void resize(std::size_t new_size);

  private:
    std::size_t _max_steps {0};
    std::vector<VectorXd> _history;
};

} // namespace prism::field
