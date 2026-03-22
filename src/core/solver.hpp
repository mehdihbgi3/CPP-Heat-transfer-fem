#pragma once

#include "mesh.hpp"
#include "material.hpp"
#include <Eigen/Sparse>
#include <vector>
#include <chrono>
#include <functional>

namespace fem {

    struct SolverStats {
        double assembly_time_ms = 0.0;
        double solve_time_ms = 0.0;
        double total_time_ms = 0.0;
        size_t iterations = 0;
        double residual = 0.0;
    };

    enum class SolverType {
        SparseLU,
        ConjugateGradient,
        BiCGSTAB
    };

    class Solver {
    public:
        explicit Solver(Mesh& mesh);

        void set_material(size_t material_id, const Material& material);
        void set_solver_type(SolverType type) { solver_type_ = type; }
        void set_tolerance(double tol) { tolerance_ = tol; }
        void set_max_iterations(int max_iter) { max_iterations_ = max_iter; }
        void enable_openmp(bool enable) { use_openmp_ = enable; }

        void solve();
        void solve_nonlinear(int max_newton_iterations = 10, double newton_tol = 1e-6);

        [[nodiscard]] const Eigen::VectorXd& temperatures() const { return temperatures_; }
        [[nodiscard]] double min_temperature() const { return temperatures_.minCoeff(); }
        [[nodiscard]] double max_temperature() const { return temperatures_.maxCoeff(); }
        [[nodiscard]] const SolverStats& stats() const { return stats_; }

        [[nodiscard]] const Eigen::SparseMatrix<double>& stiffness_matrix() const { return K_; }
        [[nodiscard]] const Eigen::VectorXd& rhs_vector() const { return rhs_; }

        void assemble_stiffness_matrix();
        void assemble_stiffness_matrix_parallel();
        void assemble_rhs();
        void apply_boundary_conditions();

    protected:
        Eigen::Matrix3d element_stiffness(size_t elem_idx) const;
        Eigen::Matrix3d element_stiffness(size_t elem_idx, const Eigen::VectorXd& T) const;
        Eigen::Matrix3d element_mass(size_t elem_idx) const;
        Eigen::Vector3d element_load(size_t elem_idx) const;

        Mesh& mesh_;
        std::vector<Material> materials_;
        Eigen::SparseMatrix<double> K_;
        Eigen::VectorXd rhs_;
        Eigen::VectorXd temperatures_;

        SolverType solver_type_ = SolverType::SparseLU;
        double tolerance_ = 1e-10;
        int max_iterations_ = 1000;
        bool use_openmp_ = false;

        SolverStats stats_;

        friend class TransientSolver;
    };

    class TransientSolver : public Solver {
    public:
        explicit TransientSolver(Mesh& mesh);

        void set_time_step(double dt) { dt_ = dt; }
        void set_end_time(double t_end) { t_end_ = t_end; }
        void set_theta(double theta) { theta_ = theta; }
        void set_initial_temperature(double T0);
        void set_initial_temperature(const Eigen::VectorXd& T0);
        void set_initial_temperature(std::function<double(double, double)> T0_func);

        void step();
        void solve();

        [[nodiscard]] double current_time() const { return current_time_; }
        [[nodiscard]] size_t current_step() const { return current_step_; }
        [[nodiscard]] bool is_finished() const { return current_time_ >= t_end_; }
        [[nodiscard]] const std::vector<Eigen::VectorXd>& temperature_history() const { return history_; }

        using StepCallback = std::function<void(size_t step, double time, const Eigen::VectorXd& T)>;
        void set_step_callback(StepCallback callback) { step_callback_ = callback; }

    private:
        void assemble_mass_matrix();
        void assemble_mass_matrix_parallel();

        Eigen::SparseMatrix<double> M_;
        double dt_ = 0.01;
        double t_end_ = 1.0;
        double theta_ = 0.5;
        double current_time_ = 0.0;
        size_t current_step_ = 0;

        std::vector<Eigen::VectorXd> history_;
        StepCallback step_callback_ = nullptr;
    };

}