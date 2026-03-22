#include "meshwidget.hpp"
#include <QPainter>
#include <QMouseEvent>
#include <cmath>
#include <algorithm>

MeshWidget::MeshWidget(QWidget* parent)
    : QOpenGLWidget(parent)
{
    setMinimumSize(400, 400);
    setMouseTracking(true);
}

void MeshWidget::setMesh(std::shared_ptr<fem::Mesh> mesh) {
    mesh_ = std::move(mesh);
    temperatures_.resize(mesh_->num_nodes());
    temperatures_.setZero();
    update();
}

void MeshWidget::setTemperatures(const Eigen::VectorXd& temperatures) {
    temperatures_ = temperatures;

    if (auto_scale_ && temperatures_.size() > 0) {
        min_temp_ = temperatures_.minCoeff();
        max_temp_ = temperatures_.maxCoeff();

        if (std::abs(max_temp_ - min_temp_) < 1e-10) {
            min_temp_ -= 1.0;
            max_temp_ += 1.0;
        }
    }

    update();
}

void MeshWidget::setTemperatureRange(double min_temp, double max_temp) {
    min_temp_ = min_temp;
    max_temp_ = max_temp;
    auto_scale_ = false;
    update();
}

void MeshWidget::initializeGL() {
    initializeOpenGLFunctions();
    glClearColor(0.15f, 0.15f, 0.15f, 1.0f);
}

void MeshWidget::resizeGL(int w, int h) {
    glViewport(0, 0, w, h);
}

QColor MeshWidget::temperatureToColor(double temp) const {
    double t = (temp - min_temp_) / (max_temp_ - min_temp_ + 1e-10);
    t = std::clamp(t, 0.0, 1.0);

    switch (colormap_) {
    case fem::ColorMap::Viridis:
        return viridis(t);
    case fem::ColorMap::Plasma:
        return plasma(t);
    case fem::ColorMap::Grayscale: {
        int v = static_cast<int>(255 * t);
        return QColor(v, v, v);
    }
    case fem::ColorMap::BlueRed:
    default: {
        int r = static_cast<int>(255 * t);
        int b = static_cast<int>(255 * (1.0 - t));
        int g = static_cast<int>(255 * (1.0 - std::abs(2.0 * t - 1.0)));
        return QColor(r, g, b);
    }
    }
}

QColor MeshWidget::viridis(double t) const {
    double r = 0.267004 + t * (0.993248 - 0.267004);
    double g = 0.004874 + t * (0.906157 - 0.004874);
    double b = 0.329415 + t * (0.143936 - 0.329415);

    if (t < 0.25) {
        r = 0.267 + t * 4 * (0.282 - 0.267);
        g = 0.005 + t * 4 * (0.140 - 0.005);
        b = 0.329 + t * 4 * (0.457 - 0.329);
    }
    else if (t < 0.5) {
        r = 0.282 + (t - 0.25) * 4 * (0.127 - 0.282);
        g = 0.140 + (t - 0.25) * 4 * (0.566 - 0.140);
        b = 0.457 + (t - 0.25) * 4 * (0.550 - 0.457);
    }
    else if (t < 0.75) {
        r = 0.127 + (t - 0.5) * 4 * (0.741 - 0.127);
        g = 0.566 + (t - 0.5) * 4 * (0.873 - 0.566);
        b = 0.550 + (t - 0.5) * 4 * (0.150 - 0.550);
    }
    else {
        r = 0.741 + (t - 0.75) * 4 * (0.993 - 0.741);
        g = 0.873 + (t - 0.75) * 4 * (0.906 - 0.873);
        b = 0.150 + (t - 0.75) * 4 * (0.144 - 0.150);
    }

    return QColor(
        static_cast<int>(255 * std::clamp(r, 0.0, 1.0)),
        static_cast<int>(255 * std::clamp(g, 0.0, 1.0)),
        static_cast<int>(255 * std::clamp(b, 0.0, 1.0))
    );
}

QColor MeshWidget::plasma(double t) const {
    double r = 0.050 + t * (0.940 - 0.050);
    double g = 0.030 + t * (0.975 - 0.030);
    double b = 0.528 + t * (0.131 - 0.528);

    if (t < 0.5) {
        r = 0.050 + t * 2 * (0.798 - 0.050);
        g = 0.030 + t * 2 * (0.280 - 0.030);
        b = 0.528 + t * 2 * (0.470 - 0.528);
    }
    else {
        r = 0.798 + (t - 0.5) * 2 * (0.940 - 0.798);
        g = 0.280 + (t - 0.5) * 2 * (0.975 - 0.280);
        b = 0.470 + (t - 0.5) * 2 * (0.131 - 0.470);
    }

    return QColor(
        static_cast<int>(255 * std::clamp(r, 0.0, 1.0)),
        static_cast<int>(255 * std::clamp(g, 0.0, 1.0)),
        static_cast<int>(255 * std::clamp(b, 0.0, 1.0))
    );
}

QPointF MeshWidget::worldToScreen(double x, double y) const {
    return QPointF(
        margin_ + (x - offset_x_) * scale_,
        height() - margin_ - (y - offset_y_) * scale_
    );
}

void MeshWidget::paintGL() {
    glClear(GL_COLOR_BUFFER_BIT);

    if (!mesh_ || mesh_->num_nodes() == 0) {
        return;
    }

    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing);

    double mesh_width = mesh_->max_x() - mesh_->min_x();
    double mesh_height = mesh_->max_y() - mesh_->min_y();

    double available_width = width() - 2 * margin_ - 80;
    double available_height = height() - 2 * margin_;

    double scale_x = available_width / mesh_width;
    double scale_y = available_height / mesh_height;
    scale_ = std::min(scale_x, scale_y);

    offset_x_ = mesh_->min_x();
    offset_y_ = mesh_->min_y();

    for (size_t e = 0; e < mesh_->num_elements(); ++e) {
        const auto& elem = mesh_->elements()[e];

        double avg_temp = 0;
        QPolygonF triangle;

        for (int i = 0; i < 3; ++i) {
            size_t idx = elem.node_indices[i];
            const auto& node = mesh_->nodes()[idx];
            triangle << worldToScreen(node.x, node.y);

            if (idx < static_cast<size_t>(temperatures_.size())) {
                avg_temp += temperatures_(idx);
            }
        }
        avg_temp /= 3.0;

        painter.setBrush(temperatureToColor(avg_temp));

        if (show_mesh_) {
            painter.setPen(QPen(QColor(30, 30, 30), 0.5));
        }
        else {
            painter.setPen(Qt::NoPen);
        }

        painter.drawPolygon(triangle);
    }

    if (show_nodes_) {
        painter.setPen(Qt::NoPen);
        painter.setBrush(Qt::white);

        for (size_t i = 0; i < mesh_->num_nodes(); ++i) {
            const auto& node = mesh_->nodes()[i];
            QPointF pos = worldToScreen(node.x, node.y);
            painter.drawEllipse(pos, 2, 2);
        }
    }

    int bar_width = 20;
    int bar_height = height() - 100;
    int bar_x = width() - 50;
    int bar_y = 50;

    for (int i = 0; i < bar_height; ++i) {
        double t = 1.0 - static_cast<double>(i) / bar_height;
        double temp = min_temp_ + t * (max_temp_ - min_temp_);
        painter.setPen(temperatureToColor(temp));
        painter.drawLine(bar_x, bar_y + i, bar_x + bar_width, bar_y + i);
    }

    painter.setPen(QPen(Qt::white, 1));
    painter.setBrush(Qt::NoBrush);
    painter.drawRect(bar_x, bar_y, bar_width, bar_height);

    painter.setPen(Qt::white);
    QFont font = painter.font();
    font.setPointSize(9);
    painter.setFont(font);

    painter.drawText(bar_x - 5, bar_y - 8, QString::number(max_temp_, 'f', 1) + " °C");
    painter.drawText(bar_x - 5, bar_y + bar_height + 15, QString::number(min_temp_, 'f', 1) + " °C");

    double mid_temp = (min_temp_ + max_temp_) / 2.0;
    painter.drawText(bar_x - 5, bar_y + bar_height / 2 + 4, QString::number(mid_temp, 'f', 1));
}

void MeshWidget::mousePressEvent(QMouseEvent* event) {
    if (!mesh_) return;

    double min_dist = std::numeric_limits<double>::max();
    size_t closest_node = 0;

    for (size_t i = 0; i < mesh_->num_nodes(); ++i) {
        const auto& node = mesh_->nodes()[i];
        QPointF screen_pos = worldToScreen(node.x, node.y);

        double dx = event->pos().x() - screen_pos.x();
        double dy = event->pos().y() - screen_pos.y();
        double dist = dx * dx + dy * dy;

        if (dist < min_dist) {
            min_dist = dist;
            closest_node = i;
        }
    }

    if (min_dist < 100) {
        const auto& node = mesh_->nodes()[closest_node];
        emit nodeClicked(closest_node, node.x, node.y);
    }

    QOpenGLWidget::mousePressEvent(event);
}