add_library(fem_core STATIC
    mesh.cpp
    solver.cpp
    material.cpp
    benchmark.cpp
    vtk_export.cpp
)

target_include_directories(fem_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(fem_core PUBLIC Eigen3::Eigen)

# OpenMP support
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(fem_core PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(fem_core PUBLIC USE_OPENMP)
endif()