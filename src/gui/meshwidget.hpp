#pragma once

#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include "mesh.hpp"
#include <Eigen/Dense>
#include <memory>

namespace fem {
    enum class ColorMap {
        BlueRed,
        Viridis,
        Plasma,
        Grayscale
    };
}

class MeshWidget : public QOpenGLWidget, protected QOpenGLFunctions {
    Q_OBJECT

public:
    explicit MeshWidget(QWidget* parent = nullptr);

    void setMesh(std::shared_ptr<fem::Mesh> mesh);
    void setTemperatures(const Eigen::VectorXd& temperatures);
    void setTemperatureRange(double min_temp, double max_temp);
    void setAutoScale(bool enable) { auto_scale_ = enable; update(); }
    void setShowMesh(bool show) { show_mesh_ = show; update(); }
    void setShowNodes(bool show) { show_nodes_ = show; update(); }
    void setColorMap(fem::ColorMap cmap) { colormap_ = cmap; update(); }

    [[nodiscard]] double minTemperature() const { return min_temp_; }
    [[nodiscard]] double maxTemperature() const { return max_temp_; }

signals:
    void nodeClicked(size_t nodeIndex, double x, double y);

protected:
    void initializeGL() override;
    void resizeGL(int w, int h) override;
    void paintGL() override;
    void mousePressEvent(QMouseEvent* event) override;

private:
    QColor temperatureToColor(double temp) const;
    QColor viridis(double t) const;
    QColor plasma(double t) const;
    QPointF worldToScreen(double x, double y) const;

    std::shared_ptr<fem::Mesh> mesh_;
    Eigen::VectorXd temperatures_;
    double min_temp_ = 0.0;
    double max_temp_ = 100.0;
    bool auto_scale_ = true;
    bool show_mesh_ = true;
    bool show_nodes_ = false;
    fem::ColorMap colormap_ = fem::ColorMap::BlueRed;

    double margin_ = 50.0;
    double scale_ = 1.0;
    double offset_x_ = 0.0;
    double offset_y_ = 0.0;
};