add_executable(duct2d main.cpp)
target_link_libraries(duct2d prism unvpp fmt::fmt-header-only nlohmann_json::nlohmann_json spdlog::spdlog)