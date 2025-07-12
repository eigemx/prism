add_executable(SIMPLESolver main.cpp)
target_link_libraries(SIMPLESolver prism fmt::fmt-header-only)

add_executable(divTest test_div.cpp)
target_link_libraries(divTest prism)