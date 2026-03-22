#pragma once

#include <functional>

namespace fem {

    struct Material {
        double thermal_conductivity;
        double density;
        double specific_heat;

        double diffusivity() const {
            return thermal_conductivity / (density * specific_heat);
        }

        std::function<double(double)> conductivity_function = nullptr;

        double get_conductivity(double temperature) const {
            if (conductivity_function) {
                return conductivity_function(temperature);
            }
            return thermal_conductivity;
        }

        static Material aluminum() {
            return { 237.0, 2700.0, 900.0 };
        }

        static Material steel() {
            return { 50.0, 7850.0, 500.0 };
        }

        static Material copper() {
            return { 400.0, 8960.0, 385.0 };
        }

        static Material aluminum_nonlinear() {
            Material mat = { 237.0, 2700.0, 900.0 };
            mat.conductivity_function = [](double T) {
                return 237.0 * (1.0 + 0.0005 * (T - 20.0));
                };
            return mat;
        }
    };

}