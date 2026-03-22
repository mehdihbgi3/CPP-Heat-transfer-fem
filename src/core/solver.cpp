#include "solver.hpp"
#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>
#include <stdexcept>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace fem {

    Solver::Solver(Mesh& mesh)
        : mesh_(mesh)
        , materials_(1, Material::aluminum())
        , K_(mesh.num_nodes(), mesh.num_nodes())
        , rhs_(Eigen::VectorXd::Zero(mesh.num_nodes()))
        , temperatures_(Eigen::VectorXd::Zero(mesh.num_nodes()))
    {
    }

    void Solver::set_material(size_t material_id, const Material& material) {
        if (material_id >= materials_.size()) {
            materials_.resize(material_id + 1, Material::aluminum());
        }
        materials_[material_id] = material;
    }

    Eigen::Matrix3d Solver::element_stiffness(size_t elem_idx) const {
        const auto& elem = mesh_.elements_[elem_idx];
        auto coords = mesh_.element_coords(elem_idx);

        double x1 = coords(0, 0), y1 = coords(0, 1);
        double x2 = coords(1, 0), y2 = coords(1, 1);
        double x3 = coords(2, 0), y3 = coords(2, 1);

        double area = mesh_.element_area(elem_idx);

        Eigen::Matrix<double, 3, 2> B;
        B(0, 0) = (y2 - y3) / (2.0 * area);
        B(0, 1) = (x3 - x2) / (2.0 * area);
        B(1, 0) = (y3 - y1) / (2.0 * area);
        B(1, 1) = (x1 - x3) / (2.0 * area);
        B(2, 0) = (y1 - y2) / (2.0 * area);
        B(2, 1) = (x2 - x1) / (2.0 * area);

        double k = materials_[elem.material_id].thermal_conductivity;

        return k * area * B * B.transpose();
    }

    Eigen::Matrix3d Solver::element_stiffness(size_t elem_idx, const Eigen::VectorXd& T) const {
        const auto& elem = mesh_.elements_[elem_idx];
        auto coords = mesh_.element_coords(elem_idx);

        double x1 = coords(0, 0), y1 = coords(0, 1);
        double x2 = coords(1, 0), y2 = coords(1, 1);
        double x3 = coords(2, 0), y3 = coords(2, 1);

        double area = mesh_.element_area(elem_idx);

        Eigen::Matrix<double, 3, 2> B;
        B(0, 0) = (y2 - y3) / (2.0 * area);
        B(0, 1) = (x3 - x2) / (2.0 * area);
        B(1, 0) = (y3 - y1) / (2.0 * area);
        B(1, 1) = (x1 - x3) / (2.0 * area);
        B(2, 0) = (y1 - y2) / (2.0 * area);
        B(2, 1) = (x2 - x1) / (2.0 * area);

        double T_avg = 0.0;
        for (int i = 0; i < 3; ++i) {
            T_avg += T(elem.node_indices[i]);
        }
        T_avg /= 3.0;

        double k = materials_[elem.material_id].get_conductivity(T_avg);

        return k * area * B * B.transpose();
    }

    Eigen::Matrix3d Solver::element_mass(size_t elem_idx) const {
        const auto& elem = mesh_.elements_[elem_idx];
        double area = mesh_.element_area(elem_idx);
        double rho = materials_[elem.material_id].density;
        double cp = materials_[elem.material_id].specific_heat;

        Eigen::Matrix3d M;
        M << 2, 1, 1,
            1, 2, 1,
            1, 1, 2;

        return (rho * cp * area / 12.0) * M;
    }

    Eigen::Vector3d Solver::element_load(size_t elem_idx) const {
        double area = mesh_.element_area(elem_idx);

        double Q = 0.0;
        for (const auto& hs : mesh_.heat_sources_) {
            if (hs.element_index == elem_idx) {
                Q = hs.value;
                break;
            }
        }

        return Eigen::Vector3d::Constant(Q * area / 3.0);
    }

    void Solver::assemble_stiffness_matrix() {
        auto start = std::chrono::high_resolution_clock::now();

        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(mesh_.num_elements() * 9);

        for (size_t e = 0; e < mesh_.num_elements(); ++e) {
            Eigen::Matrix3d Ke = element_stiffness(e);
            const auto& elem = mesh_.elements_[e];

            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    triplets.emplace_back(
                        static_cast<int>(elem.node_indices[i]),
                        static_cast<int>(elem.node_indices[j]),
                        Ke(i, j)
                    );
                }
            }
        }

        K_.setFromTriplets(triplets.begin(), triplets.end());

        auto end = std::chrono::high_resolution_clock::now();
        stats_.assembly_time_ms = std::chrono::duration<double, std::milli>(end - start).count();
    }

    void Solver::assemble_stiffness_matrix_parallel() {
#ifdef _OPENMP
        auto start = std::chrono::high_resolution_clock::now();

        int num_elements = static_cast<int>(mesh_.num_elements());
        int num_threads = omp_get_max_threads();

        std::vector<std::vector<Eigen::Triplet<double>>> thread_triplets(num_threads);

#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            thread_triplets[tid].reserve(static_cast<size_t>(num_elements) * 9 / num_threads + 100);

#pragma omp for schedule(static)
            for (int e = 0; e < num_elements; ++e) {
                Eigen::Matrix3d Ke = element_stiffness(static_cast<size_t>(e));
                const auto& elem = mesh_.elements_[static_cast<size_t>(e)];

                for (int i = 0; i < 3; ++i) {
                    for (int j = 0; j < 3; ++j) {
                        thread_triplets[tid].emplace_back(
                            static_cast<int>(elem.node_indices[i]),
                            static_cast<int>(elem.node_indices[j]),
                            Ke(i, j)
                        );
                    }
                }
            }
        }

        std::vector<Eigen::Triplet<double>> all_triplets;
        size_t total_size = 0;
        for (const auto& tt : thread_triplets) {
            total_size += tt.size();
        }
        all_triplets.reserve(total_size);

        for (const auto& tt : thread_triplets) {
            all_triplets.insert(all_triplets.end(), tt.begin(), tt.end());
        }

        K_.setFromTriplets(all_triplets.begin(), all_triplets.end());

        auto end = std::chrono::high_resolution_clock::now();
        stats_.assembly_time_ms = std::chrono::duration<double, std::milli>(end - start).count();
#else
        assemble_stiffness_matrix();
#endif
    }

    void Solver::assemble_rhs() {
        rhs_.setZero();

        for (size_t e = 0; e < mesh_.num_elements(); ++e) {
            Eigen::Vector3d fe = element_load(e);
            const auto& elem = mesh_.elements_[e];

            for (int i = 0; i < 3; ++i) {
                rhs_(elem.node_indices[i]) += fe(i);
            }
        }

        for (const auto& bc : mesh_.bcs_) {
            if (bc.type == BCType::Neumann) {
                rhs_(bc.node_index) += bc.value;
            }
            else if (bc.type == BCType::Convective) {
                rhs_(bc.node_index) += bc.h_coefficient * bc.T_ambient;
            }
        }
    }

    void Solver::apply_boundary_conditions() {
        constexpr double penalty = 1e30;

        for (size_t i = 0; i < mesh_.num_nodes(); ++i) {
            if (mesh_.nodes_[i].fixed_temperature.has_value()) {
                double T_fixed = mesh_.nodes_[i].fixed_temperature.value();
                K_.coeffRef(i, i) += penalty;
                rhs_(i) += penalty * T_fixed;
            }
        }

        for (const auto& bc : mesh_.bcs_) {
            if (bc.type == BCType::Convective) {
                K_.coeffRef(bc.node_index, bc.node_index) += bc.h_coefficient;
            }
        }
    }

    void Solver::solve() {
        auto total_start = std::chrono::high_resolution_clock::now();

        if (use_openmp_) {
            assemble_stiffness_matrix_parallel();
        }
        else {
            assemble_stiffness_matrix();
        }

        assemble_rhs();
        apply_boundary_conditions();
        K_.makeCompressed();

        auto solve_start = std::chrono::high_resolution_clock::now();

        switch (solver_type_) {
        case SolverType::SparseLU: {
            Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
            solver.compute(K_);
            if (solver.info() != Eigen::Success) {
                throw std::runtime_error("SparseLU decomposition failed");
            }
            temperatures_ = solver.solve(rhs_);
            stats_.iterations = 1;
            break;
        }

        case SolverType::ConjugateGradient: {
            Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> solver;
            solver.setTolerance(tolerance_);
            solver.setMaxIterations(max_iterations_);
            solver.compute(K_);
            temperatures_ = solver.solve(rhs_);
            stats_.iterations = solver.iterations();
            stats_.residual = solver.error();
            break;
        }

        case SolverType::BiCGSTAB: {
            Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
            solver.setTolerance(tolerance_);
            solver.setMaxIterations(max_iterations_);
            solver.compute(K_);
            temperatures_ = solver.solve(rhs_);
            stats_.iterations = solver.iterations();
            stats_.residual = solver.error();
            break;
        }
        }

        auto solve_end = std::chrono::high_resolution_clock::now();
        stats_.solve_time_ms = std::chrono::duration<double, std::milli>(solve_end - solve_start).count();

        auto total_end = std::chrono::high_resolution_clock::now();
        stats_.total_time_ms = std::chrono::duration<double, std::milli>(total_end - total_start).count();
    }

    void Solver::solve_nonlinear(int max_newton_iterations, double newton_tol) {
        temperatures_.setZero();

        for (int iter = 0; iter < max_newton_iterations; ++iter) {
            std::vector<Eigen::Triplet<double>> triplets;
            triplets.reserve(mesh_.num_elements() * 9);

            for (size_t e = 0; e < mesh_.num_elements(); ++e) {
                Eigen::Matrix3d Ke = element_stiffness(e, temperatures_);
                const auto& elem = mesh_.elements_[e];

                for (int i = 0; i < 3; ++i) {
                    for (int j = 0; j < 3; ++j) {
                        triplets.emplace_back(
                            static_cast<int>(elem.node_indices[i]),
                            static_cast<int>(elem.node_indices[j]),
                            Ke(i, j)
                        );
                    }
                }
            }

            K_.setFromTriplets(triplets.begin(), triplets.end());
            assemble_rhs();
            apply_boundary_conditions();
            K_.makeCompressed();

            Eigen::VectorXd residual = rhs_ - K_ * temperatures_;
            double res_norm = residual.norm();

            if (res_norm < newton_tol) {
                stats_.iterations = iter + 1;
                stats_.residual = res_norm;
                return;
            }

            Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
            solver.compute(K_);
            Eigen::VectorXd delta_T = solver.solve(residual);

            temperatures_ += delta_T;
        }

        stats_.iterations = max_newton_iterations;
    }

    TransientSolver::TransientSolver(Mesh& mesh)
        : Solver(mesh)
        , M_(mesh.num_nodes(), mesh.num_nodes())
    {
    }

    void TransientSolver::set_initial_temperature(double T0) {
        temperatures_.setConstant(mesh_.num_nodes(), T0);
    }

    void TransientSolver::set_initial_temperature(const Eigen::VectorXd& T0) {
        if (T0.size() != static_cast<int>(mesh_.num_nodes())) {
            throw std::invalid_argument("Initial temperature vector size mismatch");
        }
        temperatures_ = T0;
    }

    void TransientSolver::set_initial_temperature(std::function<double(double, double)> T0_func) {
        temperatures_.resize(mesh_.num_nodes());
        for (size_t i = 0; i < mesh_.num_nodes(); ++i) {
            temperatures_(i) = T0_func(mesh_.nodes_[i].x, mesh_.nodes_[i].y);
        }
    }

    void TransientSolver::assemble_mass_matrix() {
        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(mesh_.num_elements() * 9);

        for (size_t e = 0; e < mesh_.num_elements(); ++e) {
            Eigen::Matrix3d Me = element_mass(e);
            const auto& elem = mesh_.elements_[e];

            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    triplets.emplace_back(
                        static_cast<int>(elem.node_indices[i]),
                        static_cast<int>(elem.node_indices[j]),
                        Me(i, j)
                    );
                }
            }
        }

        M_.setFromTriplets(triplets.begin(), triplets.end());
    }

    void TransientSolver::assemble_mass_matrix_parallel() {
#ifdef _OPENMP
        int num_elements = static_cast<int>(mesh_.num_elements());
        int num_threads = omp_get_max_threads();

        std::vector<std::vector<Eigen::Triplet<double>>> thread_triplets(num_threads);

#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            thread_triplets[tid].reserve(static_cast<size_t>(num_elements) * 9 / num_threads + 100);

#pragma omp for schedule(static)
            for (int e = 0; e < num_elements; ++e) {
                Eigen::Matrix3d Me = element_mass(static_cast<size_t>(e));
                const auto& elem = mesh_.elements_[static_cast<size_t>(e)];

                for (int i = 0; i < 3; ++i) {
                    for (int j = 0; j < 3; ++j) {
                        thread_triplets[tid].emplace_back(
                            static_cast<int>(elem.node_indices[i]),
                            static_cast<int>(elem.node_indices[j]),
                            Me(i, j)
                        );
                    }
                }
            }
        }

        std::vector<Eigen::Triplet<double>> all_triplets;
        size_t total_size = 0;
        for (const auto& tt : thread_triplets) {
            total_size += tt.size();
        }
        all_triplets.reserve(total_size);

        for (const auto& tt : thread_triplets) {
            all_triplets.insert(all_triplets.end(), tt.begin(), tt.end());
        }

        M_.setFromTriplets(all_triplets.begin(), all_triplets.end());
#else
        assemble_mass_matrix();
#endif
    }

    void TransientSolver::step() {
        auto start = std::chrono::high_resolution_clock::now();

        if (current_step_ == 0) {
            if (use_openmp_) {
                assemble_stiffness_matrix_parallel();
                assemble_mass_matrix_parallel();
            }
            else {
                assemble_stiffness_matrix();
                assemble_mass_matrix();
            }
        }

        Eigen::SparseMatrix<double> A = M_ + theta_ * dt_ * K_;
        Eigen::SparseMatrix<double> B = M_ - (1.0 - theta_) * dt_ * K_;

        assemble_rhs();

        Eigen::VectorXd rhs_transient = B * temperatures_ + dt_ * rhs_;

        constexpr double penalty = 1e30;
        for (size_t i = 0; i < mesh_.num_nodes(); ++i) {
            if (mesh_.nodes_[i].fixed_temperature.has_value()) {
                double T_fixed = mesh_.nodes_[i].fixed_temperature.value();
                A.coeffRef(i, i) += penalty;
                rhs_transient(i) = penalty * T_fixed;
            }
        }

        A.makeCompressed();

        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(A);
        temperatures_ = solver.solve(rhs_transient);

        current_time_ += dt_;
        current_step_++;

        history_.push_back(temperatures_);

        if (step_callback_) {
            step_callback_(current_step_, current_time_, temperatures_);
        }

        auto end = std::chrono::high_resolution_clock::now();
        stats_.total_time_ms += std::chrono::duration<double, std::milli>(end - start).count();
    }

    void TransientSolver::solve() {
        history_.clear();
        history_.push_back(temperatures_);

        while (current_time_ < t_end_) {
            step();
        }
    }

}