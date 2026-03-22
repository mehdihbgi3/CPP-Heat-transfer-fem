#include <gtest/gtest.h>
#include "mesh.hpp"
#include "solver.hpp"
#include "benchmark.hpp"
#include <cmath>

TEST(SolverTest, LinearTemperatureGradient) {
    fem::Mesh mesh;
    mesh.generate_rectangle(1.0, 1.0, 10, 10);
    mesh.set_left_edge_temperature(100.0);
    mesh.set_right_edge_temperature(0.0);

    fem::Solver solver(mesh);
    solver.solve();

    const auto& T = solver.temperatures();
    const auto& nodes = mesh.nodes();

    for (size_t i = 0; i < mesh.num_nodes(); ++i) {
        double expected = 100.0 * (1.0 - nodes[i].x);
        EXPECT_NEAR(T(i), expected, 1.0);
    }
}

TEST(SolverTest, UniformTemperature) {
    fem::Mesh mesh;
    mesh.generate_rectangle(1.0, 1.0, 5, 5);

    for (size_t i = 0; i < mesh.num_nodes(); ++i) {
        mesh.set_fixed_temperature(i, 50.0);
    }

    fem::Solver solver(mesh);
    solver.solve();

    const auto& T = solver.temperatures();
    for (int i = 0; i < T.size(); ++i) {
        EXPECT_NEAR(T(i), 50.0, 1e-6);
    }
}

TEST(SolverTest, ConjugateGradientSolver) {
    fem::Mesh mesh;
    mesh.generate_rectangle(1.0, 1.0, 10, 10);
    mesh.set_left_edge_temperature(100.0);
    mesh.set_right_edge_temperature(0.0);

    fem::Solver solver(mesh);
    solver.set_solver_type(fem::SolverType::ConjugateGradient);
    solver.solve();

    EXPECT_GT(solver.max_temperature(), 90.0);
    EXPECT_LT(solver.min_temperature(), 10.0);
}

TEST(SolverTest, BiCGSTABSolver) {
    fem::Mesh mesh;
    mesh.generate_rectangle(1.0, 1.0, 10, 10);
    mesh.set_left_edge_temperature(100.0);
    mesh.set_right_edge_temperature(0.0);

    fem::Solver solver(mesh);
    solver.set_solver_type(fem::SolverType::BiCGSTAB);
    solver.solve();

    EXPECT_GT(solver.max_temperature(), 90.0);
    EXPECT_LT(solver.min_temperature(), 10.0);
}

TEST(SolverTest, DifferentMaterials) {
    fem::Mesh mesh1;
    mesh1.generate_rectangle(1.0, 1.0, 10, 10);
    mesh1.set_left_edge_temperature(100.0);
    mesh1.set_right_edge_temperature(0.0);

    fem::Solver solver1(mesh1);
    solver1.set_material(0, fem::Material::aluminum());
    solver1.solve();

    fem::Mesh mesh2;
    mesh2.generate_rectangle(1.0, 1.0, 10, 10);
    mesh2.set_left_edge_temperature(100.0);
    mesh2.set_right_edge_temperature(0.0);

    fem::Solver solver2(mesh2);
    solver2.set_material(0, fem::Material::steel());
    solver2.solve();

    for (int i = 0; i < solver1.temperatures().size(); ++i) {
        EXPECT_NEAR(solver1.temperatures()(i), solver2.temperatures()(i), 1e-6);
    }
}

TEST(SolverTest, SolverStats) {
    fem::Mesh mesh;
    mesh.generate_rectangle(1.0, 1.0, 20, 20);
    mesh.set_left_edge_temperature(100.0);
    mesh.set_right_edge_temperature(0.0);

    fem::Solver solver(mesh);
    solver.solve();

    const auto& stats = solver.stats();
    EXPECT_GT(stats.assembly_time_ms, 0.0);
    EXPECT_GT(stats.solve_time_ms, 0.0);
    EXPECT_GT(stats.total_time_ms, 0.0);
}

TEST(TransientSolverTest, ConvergesToSteadyState) {
    fem::Mesh mesh;
    mesh.generate_rectangle(1.0, 1.0, 10, 10);
    mesh.set_left_edge_temperature(100.0);
    mesh.set_right_edge_temperature(0.0);

    fem::Solver steady_solver(mesh);
    steady_solver.solve();

    fem::TransientSolver transient_solver(mesh);
    transient_solver.set_time_step(0.01);
    transient_solver.set_end_time(5.0);
    transient_solver.set_initial_temperature(50.0);
    transient_solver.solve();

    const auto& T_steady = steady_solver.temperatures();
    const auto& T_transient = transient_solver.temperatures();

    for (int i = 0; i < T_steady.size(); ++i) {
        EXPECT_NEAR(T_transient(i), T_steady(i), 5.0);
    }
}

TEST(TransientSolverTest, InitialCondition) {
    fem::Mesh mesh;
    mesh.generate_rectangle(1.0, 1.0, 5, 5);

    fem::TransientSolver solver(mesh);
    solver.set_initial_temperature(75.0);

    const auto& T = solver.temperatures();
    for (int i = 0; i < T.size(); ++i) {
        EXPECT_DOUBLE_EQ(T(i), 75.0);
    }
}

TEST(TransientSolverTest, ThetaMethods) {
    std::vector<double> thetas = { 0.0, 0.5, 1.0 };

    for (double theta : thetas) {
        fem::Mesh mesh;
        mesh.generate_rectangle(1.0, 1.0, 5, 5);
        mesh.set_left_edge_temperature(100.0);
        mesh.set_right_edge_temperature(0.0);

        fem::TransientSolver solver(mesh);
        solver.set_theta(theta);
        solver.set_time_step(0.001);
        solver.set_end_time(0.01);
        solver.set_initial_temperature(50.0);

        EXPECT_NO_THROW(solver.solve());
    }
}

TEST(BenchmarkTest, ScalingTest) {
    std::vector<size_t> resolutions = { 5, 10 };
    auto results = fem::Benchmark::run_scaling_test(resolutions, false);

    EXPECT_EQ(results.size(), 2);

    EXPECT_GT(results[1].total_time_ms, results[0].total_time_ms);
}

TEST(BenchmarkTest, SolverComparison) {
    auto results = fem::Benchmark::run_solver_comparison(10);

    EXPECT_EQ(results.size(), 3);
}