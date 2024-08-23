#pragma once

#include <spdlog/spdlog.h>

namespace prism::log {
enum class Level {
    Debug,
    Info,
    Warn,
    Error,
    Critical,
};

void inline setLevel(Level level) {
    switch (level) {
        case Level::Debug:
            spdlog::set_level(spdlog::level::debug);
            break;
        case Level::Info:
            spdlog::set_level(spdlog::level::info);
            break;
        case Level::Warn:
            spdlog::set_level(spdlog::level::warn);
            break;
        case Level::Error:
            spdlog::set_level(spdlog::level::err);
            break;
        case Level::Critical:
            spdlog::set_level(spdlog::level::critical);
            break;
    }
}

template <typename... Args>
void inline debug(spdlog::format_string_t<Args...> fmt, Args... args) {
    spdlog::debug(fmt, std::forward<Args>(args)...);
}

template <typename... Args>
void inline info(spdlog::format_string_t<Args...> fmt, Args... args) {
    spdlog::info(fmt, std::forward<Args>(args)...);
}

template <typename... Args>
void inline warn(spdlog::format_string_t<Args...> fmt, Args... args) {
    spdlog::warn(fmt, std::forward<Args>(args)...);
}

template <typename... Args>
void inline error(spdlog::format_string_t<Args...> fmt, Args... args) {
    spdlog::error(fmt, std::forward<Args>(args)...);
}

template <typename... Args>
void inline critical(spdlog::format_string_t<Args...> fmt, Args... args) {
    spdlog::critical(fmt, std::forward<Args>(args)...);
}
} // namespace prism::log