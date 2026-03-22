#include "benchmark.hpp"
#include <iostream>
#include <iomanip>
#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace fem {

    std::vector<BenchmarkResult> Benchmark::run_scaling_test(
        const std::vector<size_t>& resolutions,
        bool test_parallel)
    {
        std::vector<BenchmarkResult> results;

        for (size_t res : resolutions) {
            {
                Mesh mesh;
                mesh.generate_rectangle(1.0, 1.0, res, res);
                mesh.set_left_edge_temperature(100.0);
                mesh.set_right_edge_temperature(0.0);

                Solver solver(mesh);
                solver.enable_openmp(false);
                solver.solve();

                BenchmarkResult result;
                result.name = "Serial";
                result.num_nodes = mesh.num_nodes();
                result.num_elements = mesh.num_elements();
                result.assembly_time_ms = solver.stats().assembly_time_ms;
                result.solve_time_ms = solver.stats().solve_time_ms;
                result.total_time_ms = solver.stats().total_time_ms;
                result.iterations = solver.stats().iterations;
                result.memory_mb = estimate_memory_mb(mesh.num_nodes(), mesh.num_elements());
                result.parallel = false;
                result.num_threads = 1;

                results.push_back(result);
            }

            if (test_parallel) {
#ifdef _OPENMP
                Mesh mesh;
                mesh.generate_rectangle(1.0, 1.0, res, res);
                mesh.set_left_edge_temperature(100.0);
                mesh.set_right_edge_temperature(0.0);

                Solver solver(mesh);
                solver.enable_openmp(true);
                solver.solve();

                BenchmarkResult result;
                result.name = "Parallel";
                result.num_nodes = mesh.num_nodes();
                result.num_elements = mesh.num_elements();
                result.assembly_time_ms = solver.stats().assembly_time_ms;
                result.solve_time_ms = solver.stats().solve_time_ms;
                result.total_time_ms = solver.stats().total_time_ms;
                result.iterations = solver.stats().iterations;
                result.memory_mb = estimate_memory_mb(mesh.num_nodes(), mesh.num_elements());
                result.parallel = true;
                result.num_threads = omp_get_max_threads();

                results.push_back(result);
#endif
            }
        }

        return results;
    }

    std::vector<BenchmarkResult> Benchmark::run_solver_comparison(size_t resolution) {
        std::vector<BenchmarkResult> results;

        std::vector<std::pair<std::string, SolverType>> solvers = {
            {"SparseLU", SolverType::SparseLU},
            {"ConjugateGradient", SolverType::ConjugateGradient},
            {"BiCGSTAB", SolverType::BiCGSTAB}
        };

        for (const auto& [name, type] : solvers) {
            Mesh mesh;
            mesh.generate_rectangle(1.0, 1.0, resolution, resolution);
            mesh.set_left_edge_temperature(100.0);
            mesh.set_right_edge_temperature(0.0);

            Solver solver(mesh);
            solver.set_solver_type(type);
            solver.solve();

            BenchmarkResult result;
            result.name = name;
            result.num_nodes = mesh.num_nodes();
            result.num_elements = mesh.num_elements();
            result.assembly_time_ms = solver.stats().assembly_time_ms;
            result.solve_time_ms = solver.stats().solve_time_ms;
            result.total_time_ms = solver.stats().total_time_ms;
            result.iterations = solver.stats().iterations;
            result.memory_mb = estimate_memory_mb(mesh.num_nodes(), mesh.num_elements());
            result.parallel = false;
            result.num_threads = 1;

            results.push_back(result);
        }

        return results;
    }

    std::vector<BenchmarkResult> Benchmark::run_transient_benchmark(
        size_t resolution,
        double end_time,
        const std::vector<double>& time_steps)
    {
        std::vector<BenchmarkResult> results;

        for (double dt : time_steps) {
            Mesh mesh;
            mesh.generate_rectangle(1.0, 1.0, resolution, resolution);
            mesh.set_left_edge_temperature(100.0);
            mesh.set_right_edge_temperature(0.0);

            TransientSolver solver(mesh);
            solver.set_time_step(dt);
            solver.set_end_time(end_time);
            solver.set_initial_temperature(20.0);
            solver.solve();

            BenchmarkResult result;
            std::ostringstream oss;
            oss << "dt=" << dt;
            result.name = oss.str();
            result.num_nodes = mesh.num_nodes();
            result.num_elements = mesh.num_elements();
            result.assembly_time_ms = solver.stats().assembly_time_ms;
            result.solve_time_ms = solver.stats().solve_time_ms;
            result.total_time_ms = solver.stats().total_time_ms;
            result.iterations = solver.current_step();
            result.memory_mb = estimate_memory_mb(mesh.num_nodes(), mesh.num_elements());
            result.parallel = false;
            result.num_threads = 1;

            results.push_back(result);
        }

        return results;
    }

    void Benchmark::export_to_csv(
        const std::vector<BenchmarkResult>& results,
        const std::string& filename)
    {
        std::ofstream file(filename);
        if (!file.is_open()) return;

        file << "Name,Nodes,Elements,AssemblyTime_ms,SolveTime_ms,TotalTime_ms,"
            << "Iterations,Memory_MB,Parallel,NumThreads\n";

        for (const auto& r : results) {
            file << r.name << ","
                << r.num_nodes << ","
                << r.num_elements << ","
                << r.assembly_time_ms << ","
                << r.solve_time_ms << ","
                << r.total_time_ms << ","
                << r.iterations << ","
                << r.memory_mb << ","
                << (r.parallel ? "true" : "false") << ","
                << r.num_threads << "\n";
        }
    }

    void Benchmark::print_results(const std::vector<BenchmarkResult>& results) {
        std::cout << std::left
            << std::setw(20) << "Name"
            << std::setw(10) << "Nodes"
            << std::setw(12) << "Elements"
            << std::setw(15) << "Assembly(ms)"
            << std::setw(12) << "Solve(ms)"
            << std::setw(12) << "Total(ms)"
            << std::setw(10) << "Iters"
            << std::setw(10) << "Mem(MB)"
            << "\n";

        std::cout << std::string(101, '-') << "\n";

        for (const auto& r : results) {
            std::cout << std::left
                << std::setw(20) << r.name
                << std::setw(10) << r.num_nodes
                << std::setw(12) << r.num_elements
                << std::setw(15) << std::fixed << std::setprecision(2) << r.assembly_time_ms
                << std::setw(12) << r.solve_time_ms
                << std::setw(12) << r.total_time_ms
                << std::setw(10) << r.iterations
                << std::setw(10) << r.memory_mb
                << "\n";
        }
    }

    double Benchmark::compute_speedup(
        const BenchmarkResult& serial,
        const BenchmarkResult& parallel)
    {
        return serial.total_time_ms / parallel.total_time_ms;
    }

    double Benchmark::estimate_memory_mb(size_t num_nodes, size_t num_elements) {
        double node_mem = num_nodes * 24.0;
        double elem_mem = num_elements * 32.0;
        double matrix_mem = num_nodes * 10.0 * 8.0;
        double vector_mem = num_nodes * 8.0 * 3.0;

        return (node_mem + elem_mem + matrix_mem + vector_mem) / (1024.0 * 1024.0);
    }

}