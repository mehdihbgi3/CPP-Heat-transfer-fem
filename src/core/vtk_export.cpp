#include "vtk_export.hpp"
#include <fstream>
#include <iomanip>
#include <sstream>

namespace fem {

    bool VTKExporter::export_vtk(
        const Mesh& mesh,
        const Eigen::VectorXd& temperatures,
        const std::string& filename)
    {
        std::ofstream file(filename);
        if (!file.is_open()) return false;

        file << "# vtk DataFile Version 3.0\n";
        file << "FEM Heat Transfer Solution\n";
        file << "ASCII\n";
        file << "DATASET UNSTRUCTURED_GRID\n";

        file << "POINTS " << mesh.num_nodes() << " double\n";
        for (const auto& node : mesh.nodes()) {
            file << std::fixed << std::setprecision(6)
                << node.x << " " << node.y << " 0.0\n";
        }

        size_t num_cells = mesh.num_elements();
        size_t cell_list_size = num_cells * 4;

        file << "CELLS " << num_cells << " " << cell_list_size << "\n";
        for (const auto& elem : mesh.elements()) {
            file << "3 " << elem.node_indices[0] << " "
                << elem.node_indices[1] << " "
                << elem.node_indices[2] << "\n";
        }

        file << "CELL_TYPES " << num_cells << "\n";
        for (size_t i = 0; i < num_cells; ++i) {
            file << "5\n";
        }

        file << "POINT_DATA " << mesh.num_nodes() << "\n";
        file << "SCALARS Temperature double 1\n";
        file << "LOOKUP_TABLE default\n";
        for (int i = 0; i < temperatures.size(); ++i) {
            file << std::fixed << std::setprecision(6) << temperatures(i) << "\n";
        }

        return true;
    }

    bool VTKExporter::export_vtk_series(
        const Mesh& mesh,
        const std::vector<Eigen::VectorXd>& temperature_history,
        const std::string& base_filename,
        double dt)
    {
        for (size_t i = 0; i < temperature_history.size(); ++i) {
            std::ostringstream filename;
            filename << base_filename << "_" << std::setfill('0') << std::setw(5) << i << ".vtk";

            if (!export_vtk(mesh, temperature_history[i], filename.str())) {
                return false;
            }
        }

        std::string pvd_filename = base_filename + ".pvd";
        std::ofstream pvd(pvd_filename);
        if (!pvd.is_open()) return false;

        pvd << "<?xml version=\"1.0\"?>\n";
        pvd << "<VTKFile type=\"Collection\" version=\"0.1\">\n";
        pvd << "  <Collection>\n";

        for (size_t i = 0; i < temperature_history.size(); ++i) {
            std::ostringstream filename;
            filename << base_filename << "_" << std::setfill('0') << std::setw(5) << i << ".vtk";

            double time = i * dt;
            pvd << "    <DataSet timestep=\"" << time << "\" file=\"" << filename.str() << "\"/>\n";
        }

        pvd << "  </Collection>\n";
        pvd << "</VTKFile>\n";

        return true;
    }

    bool VTKExporter::export_mesh_vtk(
        const Mesh& mesh,
        const std::string& filename)
    {
        Eigen::VectorXd zeros = Eigen::VectorXd::Zero(mesh.num_nodes());
        return export_vtk(mesh, zeros, filename);
    }

}