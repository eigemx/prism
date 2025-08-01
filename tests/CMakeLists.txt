find_package(Catch2 CONFIG REQUIRED)

add_executable(tests 
                    test_advection_1d.cpp 
                    test_poisson.cpp 
                    test_unv_converter.cpp
                    test_scalar_history.cpp
                    test_implicit_source.cpp
               )
target_link_libraries(tests PRIVATE prism Catch2::Catch2WithMain)

include(CTest)
include(Catch)
catch_discover_tests(tests)
