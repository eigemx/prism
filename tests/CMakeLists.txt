find_package(Catch2 CONFIG REQUIRED)

#add_library(foam_json_parser foam_json_parser.cpp)
add_executable(foam_json_parser foam_json_parser.cpp)
target_link_libraries(foam_json_parser 
                            PRIVATE nlohmann_json::nlohmann_json prism fmt::fmt-header-only)

# Maybe this should be add_test?
add_executable(tests test_advection_1d.cpp test_poisson.cpp)
target_link_libraries(tests PRIVATE prism Catch2::Catch2WithMain)

add_executable(test_backward_step test_backward_step.cpp)
target_link_libraries(test_backward_step PRIVATE prism nlohmann_json::nlohmann_json fmt::fmt-header-only)

include(CTest)
include(Catch)
catch_discover_tests(tests)
