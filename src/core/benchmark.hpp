#pragma once

#include "mesh.hpp"
#include "solver.hpp"
#include <vector>
#include <string>
#include <chrono>
#include <functional>
#include <fstream>

namespace fem {

    struct BenchmarkResult {
        std::string name;
        size_t num_nodes;
        size_t num_elements;
        double assembly_time_ms;
        double solve_time_ms;
        double total_time_ms;
        size_t iterations;
        double memory_mb;
        bool parallel;
        int num_threads;
    };

    class Benchmark {
    public:
        static std::vector<BenchmarkResult> run_scaling_test(
            const std::vector<size_t>& resolutions,
            bool test_parallel = true);

        static std::vector<BenchmarkResult> run_solver_comparison(
            size_t resolution);

        static std::vector<BenchmarkResult> run_transient_benchmark(
            size_t resolution,
            double end_time,
            const std::vector<double>& time_steps);

        static void export_to_csv(
            const std::vector<BenchmarkResult>& results,
            const std::string& filename);

        static void print_results(const std::vector<BenchmarkResult>& results);

        static double compute_speedup(
            const BenchmarkResult& serial,
            const BenchmarkResult& parallel);

        static double estimate_memory_mb(size_t num_nodes, size_t num_elements);
    };

}