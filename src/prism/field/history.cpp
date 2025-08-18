#include "history.h"

namespace prism::field {

HistoryManager::HistoryManager(int num_time_steps) : _max_steps(num_time_steps) {
    _history.reserve(_max_steps);
}

auto HistoryManager::prevValues() const -> Optional<VectorXd> {
    if (_history.empty()) {
        return NullOption;
    }
    return _history.front();
}

auto HistoryManager::prevPrevValues() const -> Optional<VectorXd> {
    if (_history.size() < 2) {
        return NullOption;
    }
    return _history[1];
}

auto HistoryManager::valuesAt(std::size_t n) const -> Optional<VectorXd> {
    if (n >= _history.size()) {
        return NullOption;
    }
    return _history[n];
}

void HistoryManager::update(const VectorXd& current_values) {
    if (_max_steps == 0) {
        return;
    }
    if (_history.size() == _max_steps) {
        _history.pop_back();
    }
    _history.insert(_history.begin(), current_values);
}

void HistoryManager::update(VectorXd&& current_values) {
    if (_max_steps == 0) {
        return;
    }
    if (_history.size() == _max_steps) {
        _history.pop_back();
    }
    _history.insert(_history.begin(), std::move(current_values));
}

void HistoryManager::resize(std::size_t new_size) {
    if (new_size == _max_steps) {
        return;
    }

    // If the new size is smaller, truncate the history
    if (new_size < _history.size()) {
        _history.resize(new_size);
    }

    // Update the max steps and reserve capacity
    _max_steps = new_size;
    _history.reserve(_max_steps);
}

} // namespace prism::field
