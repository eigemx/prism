find_package(Catch2 CONFIG REQUIRED)

#add_library(foam_json_parser foam_json_parser.cpp)
add_executable(foam_json_parser foam_json_parser.cpp)
target_link_libraries(foam_json_parser PRIVATE nlohmann_json::nlohmann_json prism)

add_executable(tests test_advection_1d.cpp test_poisson.cpp)
target_link_libraries(tests PRIVATE prism Catch2::Catch2WithMain)

include(CTest)
include(Catch)
catch_discover_tests(tests)
