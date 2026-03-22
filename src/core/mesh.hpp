#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <array>
#include <memory>
#include <optional>
#include <string>
#include <functional>

namespace fem {

    struct Node {
        double x, y;
        std::optional<double> fixed_temperature;
    };

    struct Element {
        std::array<size_t, 3> node_indices;
        size_t material_id;
    };

    enum class BCType {
        Dirichlet,
        Neumann,
        Convective
    };

    struct BoundaryCondition {
        BCType type;
        size_t node_index;
        double value;
        double h_coefficient;
        double T_ambient;
    };

    struct HeatSource {
        size_t element_index;
        double value;
    };

    class Mesh {
    public:
        Mesh() = default;

        void generate_rectangle(double width, double height, size_t nx, size_t ny);
        void generate_circle(double radius, size_t n_radial, size_t n_angular);
        void generate_l_shape(double size, size_t n);

        bool load_from_json(const std::string& filename);
        bool save_to_json(const std::string& filename) const;

        [[nodiscard]] size_t num_nodes() const { return nodes_.size(); }
        [[nodiscard]] size_t num_elements() const { return elements_.size(); }
        [[nodiscard]] const std::vector<Node>& nodes() const { return nodes_; }
        [[nodiscard]] const std::vector<Element>& elements() const { return elements_; }
        [[nodiscard]] const std::vector<BoundaryCondition>& boundary_conditions() const { return bcs_; }
        [[nodiscard]] const std::vector<HeatSource>& heat_sources() const { return heat_sources_; }

        void set_fixed_temperature(size_t node_index, double temperature);
        void set_heat_flux(size_t node_index, double flux);
        void set_convective(size_t node_index, double h, double T_ambient);
        void add_heat_source(size_t element_index, double value);

        void set_left_edge_temperature(double temperature);
        void set_right_edge_temperature(double temperature);
        void set_top_edge_temperature(double temperature);
        void set_bottom_edge_temperature(double temperature);
        void set_left_edge_flux(double flux);
        void set_right_edge_flux(double flux);
        void set_left_edge_convective(double h, double T_ambient);
        void set_right_edge_convective(double h, double T_ambient);

        void clear_boundary_conditions();
        void clear_heat_sources();

        [[nodiscard]] double element_area(size_t elem_idx) const;
        [[nodiscard]] Eigen::Matrix<double, 3, 2> element_coords(size_t elem_idx) const;
        [[nodiscard]] std::vector<size_t> get_left_edge_nodes() const;
        [[nodiscard]] std::vector<size_t> get_right_edge_nodes() const;
        [[nodiscard]] std::vector<size_t> get_top_edge_nodes() const;
        [[nodiscard]] std::vector<size_t> get_bottom_edge_nodes() const;

        [[nodiscard]] double min_x() const;
        [[nodiscard]] double max_x() const;
        [[nodiscard]] double min_y() const;
        [[nodiscard]] double max_y() const;

    private:
        std::vector<Node> nodes_;
        std::vector<Element> elements_;
        std::vector<BoundaryCondition> bcs_;
        std::vector<HeatSource> heat_sources_;

        size_t nx_ = 0, ny_ = 0;

        friend class Solver;
        friend class TransientSolver;
    };

}