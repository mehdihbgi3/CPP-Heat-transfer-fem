#pragma once

#include "mesh.hpp"
#include <Eigen/Dense>
#include <string>
#include <vector>

namespace fem {

    class VTKExporter {
    public:
        static bool export_vtk(
            const Mesh& mesh,
            const Eigen::VectorXd& temperatures,
            const std::string& filename);

        static bool export_vtk_series(
            const Mesh& mesh,
            const std::vector<Eigen::VectorXd>& temperature_history,
            const std::string& base_filename,
            double dt);

        static bool export_mesh_vtk(
            const Mesh& mesh,
            const std::string& filename);
    };

}