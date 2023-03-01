add_executable(pmesh main.cpp export_vtu.h export_vtu.cpp)

target_link_libraries(pmesh prism-core fmt unvpp vtu11)