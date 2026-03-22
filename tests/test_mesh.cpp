#include <gtest/gtest.h>
#include "mesh.hpp"
#include <cmath>

TEST(MeshTest, GenerateRectangle) {
    fem::Mesh mesh;
    mesh.generate_rectangle(1.0, 1.0, 2, 2);

    EXPECT_EQ(mesh.num_nodes(), 9);
    EXPECT_EQ(mesh.num_elements(), 8);
}

TEST(MeshTest, GenerateRectangleLarge) {
    fem::Mesh mesh;
    mesh.generate_rectangle(2.0, 1.0, 100, 50);

    EXPECT_EQ(mesh.num_nodes(), 101 * 51);
    EXPECT_EQ(mesh.num_elements(), 100 * 50 * 2);
}

TEST(MeshTest, GenerateCircle) {
    fem::Mesh mesh;
    mesh.generate_circle(1.0, 5, 12);

    EXPECT_EQ(mesh.num_nodes(), 1 + 5 * 12);
    EXPECT_GT(mesh.num_elements(), 0);
}

TEST(MeshTest, GenerateLShape) {
    fem::Mesh mesh;
    mesh.generate_l_shape(1.0, 5);

    EXPECT_GT(mesh.num_nodes(), 0);
    EXPECT_GT(mesh.num_elements(), 0);
}

TEST(MeshTest, ElementArea) {
    fem::Mesh mesh;
    mesh.generate_rectangle(1.0, 1.0, 1, 1);

    EXPECT_NEAR(mesh.element_area(0), 0.5, 1e-10);
    EXPECT_NEAR(mesh.element_area(1), 0.5, 1e-10);
}

TEST(MeshTest, ElementAreaNonUniform) {
    fem::Mesh mesh;
    mesh.generate_rectangle(2.0, 1.0, 2, 1);

    for (size_t i = 0; i < mesh.num_elements(); ++i) {
        EXPECT_NEAR(mesh.element_area(i), 0.5, 1e-10);
    }
}

TEST(MeshTest, BoundaryConditions) {
    fem::Mesh mesh;
    mesh.generate_rectangle(1.0, 1.0, 2, 2);

    mesh.set_left_edge_temperature(100.0);

    const auto& nodes = mesh.nodes();
    EXPECT_TRUE(nodes[0].fixed_temperature.has_value());
    EXPECT_DOUBLE_EQ(nodes[0].fixed_temperature.value(), 100.0);
}

TEST(MeshTest, GetEdgeNodes) {
    fem::Mesh mesh;
    mesh.generate_rectangle(1.0, 1.0, 3, 3);

    auto left = mesh.get_left_edge_nodes();
    auto right = mesh.get_right_edge_nodes();
    auto top = mesh.get_top_edge_nodes();
    auto bottom = mesh.get_bottom_edge_nodes();

    EXPECT_EQ(left.size(), 4);
    EXPECT_EQ(right.size(), 4);
    EXPECT_EQ(top.size(), 4);
    EXPECT_EQ(bottom.size(), 4);
}

TEST(MeshTest, MeshBounds) {
    fem::Mesh mesh;
    mesh.generate_rectangle(2.0, 3.0, 10, 10);

    EXPECT_NEAR(mesh.min_x(), 0.0, 1e-10);
    EXPECT_NEAR(mesh.max_x(), 2.0, 1e-10);
    EXPECT_NEAR(mesh.min_y(), 0.0, 1e-10);
    EXPECT_NEAR(mesh.max_y(), 3.0, 1e-10);
}

TEST(MeshTest, ClearBoundaryConditions) {
    fem::Mesh mesh;
    mesh.generate_rectangle(1.0, 1.0, 2, 2);
    mesh.set_left_edge_temperature(100.0);

    EXPECT_GT(mesh.boundary_conditions().size(), 0);

    mesh.clear_boundary_conditions();

    EXPECT_EQ(mesh.boundary_conditions().size(), 0);
}

TEST(MeshTest, HeatSource) {
    fem::Mesh mesh;
    mesh.generate_rectangle(1.0, 1.0, 2, 2);

    mesh.add_heat_source(0, 1000.0);

    EXPECT_EQ(mesh.heat_sources().size(), 1);
    EXPECT_EQ(mesh.heat_sources()[0].element_index, 0);
    EXPECT_DOUBLE_EQ(mesh.heat_sources()[0].value, 1000.0);
}