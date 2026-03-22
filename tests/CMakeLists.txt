add_executable(fem_tests
    test_mesh.cpp
    test_solver.cpp
)

target_link_libraries(fem_tests PRIVATE
    fem_core
    GTest::gtest
    GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(fem_tests)