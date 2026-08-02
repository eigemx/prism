#pragma once

#include <vector>

#include "prism/types.h"

namespace prism::field {

class HistoryManager {
  public:
    /**
     * @brief Constructs a HistoryManager with the specified maximum number of time steps.
     * @param num_time_steps The maximum number of time steps to store in history.
     */
    HistoryManager(size_t num_time_steps);

    /**
     * @brief Returns a reference to the values from the previous time step.
     * @return Optional containing the previous values, or NullOption if history is empty.
     */
    auto prevValues() const -> Optional<Ref<const VectorXd>>;

    /**
     * @brief Returns a reference to the values from two time steps ago.
     * @return Optional containing the values from two steps ago, or NullOption if history size
     * < 2.
     */
    auto prevPrevValues() const -> Optional<Ref<const VectorXd>>;

    /**
     * @brief Returns a reference to the values at a specific time step index.
     * @param n The index of the time step (0 = most recent).
     * @return Optional containing the values at index n, or NullOption if out of bounds.
     */
    auto valuesAt(size_t n) const -> Optional<Ref<const VectorXd>>;

    /**
     * @brief Updates the history with new values from the current time step.
     * @param current_values The values to add to the history.
     */
    void update(const VectorXd& current_values);

    /**
     * @brief Updates the history with new values from the current time step (move version).
     * @param current_values The values to add to the history (moved).
     */
    void update(VectorXd&& current_values);

    /**
     * @brief Resizes the history to a new maximum number of time steps.
     * @param new_size The new maximum number of time steps to store.
     */
    void resize(size_t new_size);

  private:
    size_t _max_steps {0};
    std::vector<VectorXd> _history;
};

} // namespace prism::field
