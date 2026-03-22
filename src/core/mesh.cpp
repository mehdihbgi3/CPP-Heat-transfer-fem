#include "mesh.hpp"
#include <stdexcept>

#define _USE_MATH_DEFINES
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <fstream>
#include <algorithm>
#include <limits>

namespace fem {

    void Mesh::generate_rectangle(double width, double height, size_t nx, size_t ny) {
        nodes_.clear();
        elements_.clear();
        bcs_.clear();
        heat_sources_.clear();

        nx_ = nx;
        ny_ = ny;

        double dx = width / static_cast<double>(nx);
        double dy = height / static_cast<double>(ny);

        for (size_t j = 0; j <= ny; ++j) {
            for (size_t i = 0; i <= nx; ++i) {
                nodes_.push_back({
                    static_cast<double>(i) * dx,
                    static_cast<double>(j) * dy,
                    std::nullopt
                    });
            }
        }

        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {
                size_t n0 = j * (nx + 1) + i;
                size_t n1 = n0 + 1;
                size_t n2 = n0 + (nx + 1);
                size_t n3 = n2 + 1;

                elements_.push_back({ {n0, n1, n2}, 0 });
                elements_.push_back({ {n1, n3, n2}, 0 });
            }
        }
    }

    void Mesh::generate_circle(double radius, size_t n_radial, size_t n_angular) {
        nodes_.clear();
        elements_.clear();
        bcs_.clear();
        heat_sources_.clear();

        nodes_.push_back({ 0.0, 0.0, std::nullopt });

        for (size_t r = 1; r <= n_radial; ++r) {
            double rad = radius * static_cast<double>(r) / static_cast<double>(n_radial);
            for (size_t a = 0; a < n_angular; ++a) {
                double angle = 2.0 * M_PI * static_cast<double>(a) / static_cast<double>(n_angular);
                nodes_.push_back({ rad * std::cos(angle), rad * std::sin(angle), std::nullopt });
            }
        }

        for (size_t a = 0; a < n_angular; ++a) {
            size_t next = (a + 1) % n_angular;
            elements_.push_back({ {0, 1 + a, 1 + next}, 0 });
        }

        for (size_t r = 1; r < n_radial; ++r) {
            size_t inner_start = 1 + (r - 1) * n_angular;
            size_t outer_start = 1 + r * n_angular;

            for (size_t a = 0; a < n_angular; ++a) {
                size_t next = (a + 1) % n_angular;

                size_t i0 = inner_start + a;
                size_t i1 = inner_start + next;
                size_t o0 = outer_start + a;
                size_t o1 = outer_start + next;

                elements_.push_back({ {i0, o0, i1}, 0 });
                elements_.push_back({ {i1, o0, o1}, 0 });
            }
        }
    }

    void Mesh::generate_l_shape(double size, size_t n) {
        nodes_.clear();
        elements_.clear();
        bcs_.clear();
        heat_sources_.clear();

        double h = size / static_cast<double>(n);

        for (size_t j = n; j <= 2 * n; ++j) {
            for (size_t i = 0; i <= 2 * n; ++i) {
                nodes_.push_back({
                    static_cast<double>(i) * h,
                    static_cast<double>(j) * h,
                    std::nullopt
                    });
            }
        }

        for (size_t j = 0; j < n; ++j) {
            for (size_t i = n; i <= 2 * n; ++i) {
                nodes_.push_back({
                    static_cast<double>(i) * h,
                    static_cast<double>(j) * h,
                    std::nullopt
                    });
            }
        }

        size_t cols = 2 * n + 1;
        for (size_t j = 0; j < n; ++j) {
            for (size_t i = 0; i < 2 * n; ++i) {
                size_t n0 = j * cols + i;
                size_t n1 = n0 + 1;
                size_t n2 = n0 + cols;
                size_t n3 = n2 + 1;

                elements_.push_back({ {n0, n1, n2}, 0 });
                elements_.push_back({ {n1, n3, n2}, 0 });
            }
        }

        size_t upper_nodes = (n + 1) * cols;
        size_t lower_cols = n + 1;

        for (size_t j = 0; j < n; ++j) {
            for (size_t i = 0; i < n; ++i) {
                size_t n0, n1, n2, n3;

                if (j == n - 1) {
                    n2 = n * cols + (n + i);
                    n3 = n * cols + (n + i + 1);
                }
                else {
                    n2 = upper_nodes + (j + 1) * lower_cols + i;
                    n3 = n2 + 1;
                }

                n0 = upper_nodes + j * lower_cols + i;
                n1 = n0 + 1;

                elements_.push_back({ {n0, n1, n2}, 0 });
                elements_.push_back({ {n1, n3, n2}, 0 });
            }
        }
    }

    void Mesh::set_fixed_temperature(size_t node_index, double temperature) {
        if (node_index >= nodes_.size()) {
            throw std::out_of_range("Node index out of range");
        }
        nodes_[node_index].fixed_temperature = temperature;
        bcs_.push_back({ BCType::Dirichlet, node_index, temperature, 0.0, 0.0 });
    }

    void Mesh::set_heat_flux(size_t node_index, double flux) {
        if (node_index >= nodes_.size()) {
            throw std::out_of_range("Node index out of range");
        }
        bcs_.push_back({ BCType::Neumann, node_index, flux, 0.0, 0.0 });
    }

    void Mesh::set_convective(size_t node_index, double h, double T_ambient) {
        if (node_index >= nodes_.size()) {
            throw std::out_of_range("Node index out of range");
        }
        bcs_.push_back({ BCType::Convective, node_index, 0.0, h, T_ambient });
    }

    void Mesh::add_heat_source(size_t element_index, double value) {
        if (element_index >= elements_.size()) {
            throw std::out_of_range("Element index out of range");
        }
        heat_sources_.push_back({ element_index, value });
    }

    std::vector<size_t> Mesh::get_left_edge_nodes() const {
        std::vector<size_t> nodes;
        double min_x_val = min_x();
        for (size_t i = 0; i < nodes_.size(); ++i) {
            if (std::abs(nodes_[i].x - min_x_val) < 1e-10) {
                nodes.push_back(i);
            }
        }
        return nodes;
    }

    std::vector<size_t> Mesh::get_right_edge_nodes() const {
        std::vector<size_t> nodes;
        double max_x_val = max_x();
        for (size_t i = 0; i < nodes_.size(); ++i) {
            if (std::abs(nodes_[i].x - max_x_val) < 1e-10) {
                nodes.push_back(i);
            }
        }
        return nodes;
    }

    std::vector<size_t> Mesh::get_top_edge_nodes() const {
        std::vector<size_t> nodes;
        double max_y_val = max_y();
        for (size_t i = 0; i < nodes_.size(); ++i) {
            if (std::abs(nodes_[i].y - max_y_val) < 1e-10) {
                nodes.push_back(i);
            }
        }
        return nodes;
    }

    std::vector<size_t> Mesh::get_bottom_edge_nodes() const {
        std::vector<size_t> nodes;
        double min_y_val = min_y();
        for (size_t i = 0; i < nodes_.size(); ++i) {
            if (std::abs(nodes_[i].y - min_y_val) < 1e-10) {
                nodes.push_back(i);
            }
        }
        return nodes;
    }

    void Mesh::set_left_edge_temperature(double temperature) {
        for (size_t idx : get_left_edge_nodes()) {
            set_fixed_temperature(idx, temperature);
        }
    }

    void Mesh::set_right_edge_temperature(double temperature) {
        for (size_t idx : get_right_edge_nodes()) {
            set_fixed_temperature(idx, temperature);
        }
    }

    void Mesh::set_top_edge_temperature(double temperature) {
        for (size_t idx : get_top_edge_nodes()) {
            set_fixed_temperature(idx, temperature);
        }
    }

    void Mesh::set_bottom_edge_temperature(double temperature) {
        for (size_t idx : get_bottom_edge_nodes()) {
            set_fixed_temperature(idx, temperature);
        }
    }

    void Mesh::set_left_edge_flux(double flux) {
        for (size_t idx : get_left_edge_nodes()) {
            set_heat_flux(idx, flux);
        }
    }

    void Mesh::set_right_edge_flux(double flux) {
        for (size_t idx : get_right_edge_nodes()) {
            set_heat_flux(idx, flux);
        }
    }

    void Mesh::set_left_edge_convective(double h, double T_ambient) {
        for (size_t idx : get_left_edge_nodes()) {
            set_convective(idx, h, T_ambient);
        }
    }

    void Mesh::set_right_edge_convective(double h, double T_ambient) {
        for (size_t idx : get_right_edge_nodes()) {
            set_convective(idx, h, T_ambient);
        }
    }

    void Mesh::clear_boundary_conditions() {
        bcs_.clear();
        for (auto& node : nodes_) {
            node.fixed_temperature = std::nullopt;
        }
    }

    void Mesh::clear_heat_sources() {
        heat_sources_.clear();
    }

    double Mesh::element_area(size_t elem_idx) const {
        auto coords = element_coords(elem_idx);

        double x1 = coords(0, 0), y1 = coords(0, 1);
        double x2 = coords(1, 0), y2 = coords(1, 1);
        double x3 = coords(2, 0), y3 = coords(2, 1);

        return 0.5 * std::abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));
    }

    Eigen::Matrix<double, 3, 2> Mesh::element_coords(size_t elem_idx) const {
        const auto& elem = elements_[elem_idx];
        Eigen::Matrix<double, 3, 2> coords;

        for (int i = 0; i < 3; ++i) {
            coords(i, 0) = nodes_[elem.node_indices[i]].x;
            coords(i, 1) = nodes_[elem.node_indices[i]].y;
        }

        return coords;
    }

    double Mesh::min_x() const {
        double val = std::numeric_limits<double>::max();
        for (const auto& node : nodes_) {
            val = std::min(val, node.x);
        }
        return val;
    }

    double Mesh::max_x() const {
        double val = std::numeric_limits<double>::lowest();
        for (const auto& node : nodes_) {
            val = std::max(val, node.x);
        }
        return val;
    }

    double Mesh::min_y() const {
        double val = std::numeric_limits<double>::max();
        for (const auto& node : nodes_) {
            val = std::min(val, node.y);
        }
        return val;
    }

    double Mesh::max_y() const {
        double val = std::numeric_limits<double>::lowest();
        for (const auto& node : nodes_) {
            val = std::max(val, node.y);
        }
        return val;
    }

    bool Mesh::load_from_json(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) return false;

        return false;
    }

    bool Mesh::save_to_json(const std::string& filename) const {
        std::ofstream file(filename);
        if (!file.is_open()) return false;

        file << "{\n";
        file << "  \"nodes\": [\n";
        for (size_t i = 0; i < nodes_.size(); ++i) {
            file << "    {\"x\": " << nodes_[i].x << ", \"y\": " << nodes_[i].y << "}";
            if (i < nodes_.size() - 1) file << ",";
            file << "\n";
        }
        file << "  ],\n";
        file << "  \"elements\": [\n";
        for (size_t i = 0; i < elements_.size(); ++i) {
            file << "    [" << elements_[i].node_indices[0] << ", "
                << elements_[i].node_indices[1] << ", "
                << elements_[i].node_indices[2] << "]";
            if (i < elements_.size() - 1) file << ",";
            file << "\n";
        }
        file << "  ]\n";
        file << "}\n";

        return true;
    }

}