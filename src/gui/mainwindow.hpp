#pragma once

#include <QMainWindow>
#include <QSlider>
#include <QLabel>
#include <QPushButton>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QComboBox>
#include <QCheckBox>
#include <QTimer>
#include <QProgressBar>
#include <QTabWidget>
#include <QTextEdit>
#include "meshwidget.hpp"
#include "mesh.hpp"
#include "solver.hpp"
#include "benchmark.hpp"
#include <memory>

class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    MainWindow(QWidget* parent = nullptr);

private slots:
    void onSolve();
    void onMeshResolutionChanged(int value);
    void onLeftTempChanged(int value);
    void onRightTempChanged(int value);
    void onMaterialChanged(int index);
    void onGeometryChanged(int index);
    void onSolverTypeChanged(int index);

    void onStartTransient();
    void onStopTransient();
    void onResetTransient();
    void onAnimationStep();

    void onColorMapChanged(int index);
    void onShowMeshToggled(bool checked);
    void onShowNodesToggled(bool checked);

    void onRunBenchmarks();
    void onExportResults();

    void onBCTypeChanged(int index);
    void onNodeClicked(size_t nodeIndex, double x, double y);

private:
    void setupUI();
    void setupSteadyStateTab();
    void setupTransientTab();
    void setupVisualizationTab();
    void setupBenchmarkTab();
    void createMesh();
    void updateStatusBar();
    void appendLog(const QString& message);

    MeshWidget* mesh_widget_;
    std::shared_ptr<fem::Mesh> mesh_;
    std::unique_ptr<fem::Solver> solver_;
    std::unique_ptr<fem::TransientSolver> transient_solver_;

    QTabWidget* tab_widget_;

    QSpinBox* resolution_spinbox_;
    QComboBox* geometry_combo_;
    QComboBox* material_combo_;
    QComboBox* solver_type_combo_;
    QSlider* left_temp_slider_;
    QSlider* right_temp_slider_;
    QLabel* left_temp_label_;
    QLabel* right_temp_label_;
    QPushButton* solve_button_;
    QCheckBox* parallel_checkbox_;

    QDoubleSpinBox* dt_spinbox_;
    QDoubleSpinBox* end_time_spinbox_;
    QDoubleSpinBox* initial_temp_spinbox_;
    QComboBox* theta_combo_;
    QPushButton* start_transient_button_;
    QPushButton* stop_transient_button_;
    QPushButton* reset_transient_button_;
    QSlider* animation_speed_slider_;
    QLabel* time_label_;
    QProgressBar* progress_bar_;
    QTimer* animation_timer_;

    QComboBox* colormap_combo_;
    QCheckBox* show_mesh_checkbox_;
    QCheckBox* show_nodes_checkbox_;
    QCheckBox* auto_scale_checkbox_;

    QPushButton* run_benchmark_button_;
    QPushButton* export_results_button_;
    QTextEdit* benchmark_log_;
    std::vector<fem::BenchmarkResult> benchmark_results_;

    QComboBox* bc_type_combo_;
    QDoubleSpinBox* bc_value_spinbox_;
    QDoubleSpinBox* bc_h_spinbox_;
    QDoubleSpinBox* bc_ambient_spinbox_;

    QLabel* status_label_;

    bool transient_running_ = false;
    size_t current_animation_step_ = 0;
    std::vector<Eigen::VectorXd> temperature_history_;
};