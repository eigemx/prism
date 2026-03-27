#pragma once

#include <vector>

#include "prism/types.h"

namespace prism::field {

/// TODO: the prevValues getters are not efficient, because they return copies of the data.
/// We should return references instead. std::option does not allow returning references. We
/// need to find an alternative. Consider std::reference_wrapper

class HistoryManager {
  public:
    /**
     * @brief Constructs a HistoryManager with the specified maximum number of time steps.
     * @param num_time_steps The maximum number of time steps to store in history.
     */
    HistoryManager(std::size_t num_time_steps);

    /**
     * @brief Returns the values from the previous time step.
     * @return Optional containing the previous values, or NullOption if history is empty.
     */
    auto prevValues() const -> Optional<VectorXd>;

    /**
     * @brief Returns the values from two time steps ago.
     * @return Optional containing the values from two steps ago, or NullOption if history size
     * < 2.
     */
    auto prevPrevValues() const -> Optional<VectorXd>;

    /**
     * @brief Returns the values at a specific time step index.
     * @param n The index of the time step (0 = most recent).
     * @return Optional containing the values at index n, or NullOption if out of bounds.
     */
    auto valuesAt(std::size_t n) const -> Optional<VectorXd>;

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
    void resize(std::size_t new_size);

  private:
    std::size_t _max_steps {0};
    std::vector<VectorXd> _history;
};

} // namespace prism::field
