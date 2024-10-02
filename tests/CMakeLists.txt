find_package(Catch2 CONFIG REQUIRED)

add_executable(tests test_advection_1d.cpp test_poisson.cpp)
target_link_libraries(tests PRIVATE prism Catch2::Catch2WithMain)

add_executable(test_backward_step test_backward_step.cpp)
target_link_libraries(test_backward_step PRIVATE prism fmt::fmt-header-only nlohmann_json::nlohmann_json)

include(CTest)
include(Catch)
catch_discover_tests(tests)
