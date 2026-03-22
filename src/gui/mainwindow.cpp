#include "mainwindow.hpp"
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QGroupBox>
#include <QFormLayout>
#include <QStatusBar>
#include <QFileDialog>
#include <QMessageBox>
#include <QScrollArea>
#include <sstream>
#include <iomanip>

MainWindow::MainWindow(QWidget* parent)
    : QMainWindow(parent)
    , animation_timer_(new QTimer(this))
{
    setupUI();
    createMesh();
    onSolve();

    connect(animation_timer_, &QTimer::timeout, this, &MainWindow::onAnimationStep);
}

void MainWindow::setupUI() {
    setWindowTitle("FEM Heat Transfer Visualizer");
    setMinimumSize(1200, 800);

    auto* central_widget = new QWidget(this);
    setCentralWidget(central_widget);

    auto* main_layout = new QHBoxLayout(central_widget);

    mesh_widget_ = new MeshWidget(this);
    main_layout->addWidget(mesh_widget_, 1);

    connect(mesh_widget_, &MeshWidget::nodeClicked, this, &MainWindow::onNodeClicked);

    auto* control_panel = new QWidget(this);
    control_panel->setFixedWidth(320);
    auto* control_layout = new QVBoxLayout(control_panel);

    tab_widget_ = new QTabWidget(this);

    setupSteadyStateTab();
    setupTransientTab();
    setupVisualizationTab();
    setupBenchmarkTab();

    control_layout->addWidget(tab_widget_);
    main_layout->addWidget(control_panel);

    status_label_ = new QLabel(this);
    statusBar()->addPermanentWidget(status_label_);
}

void MainWindow::setupSteadyStateTab() {
    auto* tab = new QWidget();
    auto* layout = new QVBoxLayout(tab);

    auto* mesh_group = new QGroupBox("Mesh Settings", this);
    auto* mesh_layout = new QFormLayout(mesh_group);

    resolution_spinbox_ = new QSpinBox(this);
    resolution_spinbox_->setRange(5, 200);
    resolution_spinbox_->setValue(20);
    mesh_layout->addRow("Resolution:", resolution_spinbox_);

    geometry_combo_ = new QComboBox(this);
    geometry_combo_->addItem("Rectangle");
    geometry_combo_->addItem("Circle");
    geometry_combo_->addItem("L-Shape");
    mesh_layout->addRow("Geometry:", geometry_combo_);

    layout->addWidget(mesh_group);

    auto* material_group = new QGroupBox("Material", this);
    auto* material_layout = new QFormLayout(material_group);

    material_combo_ = new QComboBox(this);
    material_combo_->addItem("Aluminum (k=237)");
    material_combo_->addItem("Steel (k=50)");
    material_combo_->addItem("Copper (k=400)");
    material_combo_->addItem("Aluminum (Nonlinear)");
    material_layout->addRow("Type:", material_combo_);

    layout->addWidget(material_group);

    auto* solver_group = new QGroupBox("Solver", this);
    auto* solver_layout = new QFormLayout(solver_group);

    solver_type_combo_ = new QComboBox(this);
    solver_type_combo_->addItem("SparseLU (Direct)");
    solver_type_combo_->addItem("Conjugate Gradient");
    solver_type_combo_->addItem("BiCGSTAB");
    solver_layout->addRow("Type:", solver_type_combo_);

    parallel_checkbox_ = new QCheckBox("Enable OpenMP", this);
    solver_layout->addRow(parallel_checkbox_);

    layout->addWidget(solver_group);

    auto* bc_group = new QGroupBox("Boundary Conditions", this);
    auto* bc_layout = new QVBoxLayout(bc_group);

    left_temp_label_ = new QLabel("Left Edge: 100 °C", this);
    left_temp_slider_ = new QSlider(Qt::Horizontal, this);
    left_temp_slider_->setRange(0, 500);
    left_temp_slider_->setValue(100);

    right_temp_label_ = new QLabel("Right Edge: 0 °C", this);
    right_temp_slider_ = new QSlider(Qt::Horizontal, this);
    right_temp_slider_->setRange(0, 500);
    right_temp_slider_->setValue(0);

    bc_layout->addWidget(left_temp_label_);
    bc_layout->addWidget(left_temp_slider_);
    bc_layout->addWidget(right_temp_label_);
    bc_layout->addWidget(right_temp_slider_);

    layout->addWidget(bc_group);

    solve_button_ = new QPushButton("Solve", this);
    solve_button_->setStyleSheet("QPushButton { padding: 10px; font-weight: bold; background-color: #4CAF50; color: white; }");
    layout->addWidget(solve_button_);

    layout->addStretch();

    connect(solve_button_, &QPushButton::clicked, this, &MainWindow::onSolve);
    connect(resolution_spinbox_, QOverload<int>::of(&QSpinBox::valueChanged), this, &MainWindow::onMeshResolutionChanged);
    connect(geometry_combo_, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &MainWindow::onGeometryChanged);
    connect(material_combo_, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &MainWindow::onMaterialChanged);
    connect(solver_type_combo_, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &MainWindow::onSolverTypeChanged);
    connect(left_temp_slider_, &QSlider::valueChanged, this, &MainWindow::onLeftTempChanged);
    connect(right_temp_slider_, &QSlider::valueChanged, this, &MainWindow::onRightTempChanged);

    tab_widget_->addTab(tab, "Steady State");
}

void MainWindow::setupTransientTab() {
    auto* tab = new QWidget();
    auto* layout = new QVBoxLayout(tab);

    auto* time_group = new QGroupBox("Time Settings", this);
    auto* time_layout = new QFormLayout(time_group);

    dt_spinbox_ = new QDoubleSpinBox(this);
    dt_spinbox_->setRange(0.0001, 1.0);
    dt_spinbox_->setDecimals(4);
    dt_spinbox_->setValue(0.01);
    dt_spinbox_->setSingleStep(0.001);
    time_layout->addRow("Time Step (s):", dt_spinbox_);

    end_time_spinbox_ = new QDoubleSpinBox(this);
    end_time_spinbox_->setRange(0.1, 100.0);
    end_time_spinbox_->setValue(1.0);
    time_layout->addRow("End Time (s):", end_time_spinbox_);

    theta_combo_ = new QComboBox(this);
    theta_combo_->addItem("Explicit (θ=0)");
    theta_combo_->addItem("Crank-Nicolson (θ=0.5)");
    theta_combo_->addItem("Implicit (θ=1)");
    theta_combo_->setCurrentIndex(1);
    time_layout->addRow("Method:", theta_combo_);

    layout->addWidget(time_group);

    auto* ic_group = new QGroupBox("Initial Condition", this);
    auto* ic_layout = new QFormLayout(ic_group);

    initial_temp_spinbox_ = new QDoubleSpinBox(this);
    initial_temp_spinbox_->setRange(-100, 1000);
    initial_temp_spinbox_->setValue(20.0);
    ic_layout->addRow("Initial Temp (°C):", initial_temp_spinbox_);

    layout->addWidget(ic_group);

    auto* anim_group = new QGroupBox("Animation", this);
    auto* anim_layout = new QVBoxLayout(anim_group);

    auto* button_layout = new QHBoxLayout();
    start_transient_button_ = new QPushButton("▶ Start", this);
    stop_transient_button_ = new QPushButton("⏹ Stop", this);
    reset_transient_button_ = new QPushButton("↺ Reset", this);
    stop_transient_button_->setEnabled(false);

    button_layout->addWidget(start_transient_button_);
    button_layout->addWidget(stop_transient_button_);
    button_layout->addWidget(reset_transient_button_);
    anim_layout->addLayout(button_layout);

    auto* speed_layout = new QHBoxLayout();
    speed_layout->addWidget(new QLabel("Speed:", this));
    animation_speed_slider_ = new QSlider(Qt::Horizontal, this);
    animation_speed_slider_->setRange(1, 100);
    animation_speed_slider_->setValue(50);
    speed_layout->addWidget(animation_speed_slider_);
    anim_layout->addLayout(speed_layout);

    time_label_ = new QLabel("Time: 0.000 s", this);
    anim_layout->addWidget(time_label_);

    progress_bar_ = new QProgressBar(this);
    progress_bar_->setRange(0, 100);
    progress_bar_->setValue(0);
    anim_layout->addWidget(progress_bar_);

    layout->addWidget(anim_group);

    layout->addStretch();

    connect(start_transient_button_, &QPushButton::clicked, this, &MainWindow::onStartTransient);
    connect(stop_transient_button_, &QPushButton::clicked, this, &MainWindow::onStopTransient);
    connect(reset_transient_button_, &QPushButton::clicked, this, &MainWindow::onResetTransient);

    tab_widget_->addTab(tab, "Transient");
}

void MainWindow::setupVisualizationTab() {
    auto* tab = new QWidget();
    auto* layout = new QVBoxLayout(tab);

    auto* vis_group = new QGroupBox("Display Options", this);
    auto* vis_layout = new QFormLayout(vis_group);

    colormap_combo_ = new QComboBox(this);
    colormap_combo_->addItem("Blue-Red");
    colormap_combo_->addItem("Viridis");
    colormap_combo_->addItem("Plasma");
    colormap_combo_->addItem("Grayscale");
    vis_layout->addRow("Colormap:", colormap_combo_);

    show_mesh_checkbox_ = new QCheckBox("Show Mesh Lines", this);
    show_mesh_checkbox_->setChecked(true);
    vis_layout->addRow(show_mesh_checkbox_);

    show_nodes_checkbox_ = new QCheckBox("Show Nodes", this);
    vis_layout->addRow(show_nodes_checkbox_);

    auto_scale_checkbox_ = new QCheckBox("Auto-scale Colors", this);
    auto_scale_checkbox_->setChecked(true);
    vis_layout->addRow(auto_scale_checkbox_);

    layout->addWidget(vis_group);
    layout->addStretch();

    connect(colormap_combo_, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &MainWindow::onColorMapChanged);
    connect(show_mesh_checkbox_, &QCheckBox::toggled, this, &MainWindow::onShowMeshToggled);
    connect(show_nodes_checkbox_, &QCheckBox::toggled, this, &MainWindow::onShowNodesToggled);
    connect(auto_scale_checkbox_, &QCheckBox::toggled, mesh_widget_, &MeshWidget::setAutoScale);

    tab_widget_->addTab(tab, "Visualization");
}

void MainWindow::setupBenchmarkTab() {
    auto* tab = new QWidget();
    auto* layout = new QVBoxLayout(tab);

    auto* bench_group = new QGroupBox("Benchmarks", this);
    auto* bench_layout = new QVBoxLayout(bench_group);

    run_benchmark_button_ = new QPushButton("Run Scaling Test", this);
    export_results_button_ = new QPushButton("Export Results (CSV)", this);
    export_results_button_->setEnabled(false);

    bench_layout->addWidget(run_benchmark_button_);
    bench_layout->addWidget(export_results_button_);

    layout->addWidget(bench_group);

    benchmark_log_ = new QTextEdit(this);
    benchmark_log_->setReadOnly(true);
    benchmark_log_->setFont(QFont("Consolas", 9));
    layout->addWidget(benchmark_log_);

    connect(run_benchmark_button_, &QPushButton::clicked, this, &MainWindow::onRunBenchmarks);
    connect(export_results_button_, &QPushButton::clicked, this, &MainWindow::onExportResults);

    tab_widget_->addTab(tab, "Benchmarks");
}

void MainWindow::createMesh() {
    mesh_ = std::make_shared<fem::Mesh>();
    int res = resolution_spinbox_->value();

    switch (geometry_combo_->currentIndex()) {
    case 0:
        mesh_->generate_rectangle(1.0, 1.0, res, res);
        break;
    case 1:
        mesh_->generate_circle(0.5, res / 2, res);
        break;
    case 2:
        mesh_->generate_l_shape(1.0, res / 2);
        break;
    }

    mesh_->set_left_edge_temperature(left_temp_slider_->value());
    mesh_->set_right_edge_temperature(right_temp_slider_->value());
    mesh_widget_->setMesh(mesh_);

    solver_ = std::make_unique<fem::Solver>(*mesh_);
}

void MainWindow::onSolve() {
    createMesh();

    switch (material_combo_->currentIndex()) {
    case 0: solver_->set_material(0, fem::Material::aluminum()); break;
    case 1: solver_->set_material(0, fem::Material::steel()); break;
    case 2: solver_->set_material(0, fem::Material::copper()); break;
    case 3: solver_->set_material(0, fem::Material::aluminum_nonlinear()); break;
    }

    switch (solver_type_combo_->currentIndex()) {
    case 0: solver_->set_solver_type(fem::SolverType::SparseLU); break;
    case 1: solver_->set_solver_type(fem::SolverType::ConjugateGradient); break;
    case 2: solver_->set_solver_type(fem::SolverType::BiCGSTAB); break;
    }

    solver_->enable_openmp(parallel_checkbox_->isChecked());

    if (material_combo_->currentIndex() == 3) {
        solver_->solve_nonlinear();
    }
    else {
        solver_->solve();
    }

    mesh_widget_->setTemperatures(solver_->temperatures());
    updateStatusBar();
}

void MainWindow::onMeshResolutionChanged(int) {
    onSolve();
}

void MainWindow::onGeometryChanged(int) {
    onSolve();
}

void MainWindow::onLeftTempChanged(int value) {
    left_temp_label_->setText(QString("Left Edge: %1 °C").arg(value));
    onSolve();
}

void MainWindow::onRightTempChanged(int value) {
    right_temp_label_->setText(QString("Right Edge: %1 °C").arg(value));
    onSolve();
}

void MainWindow::onMaterialChanged(int) {
    onSolve();
}

void MainWindow::onSolverTypeChanged(int) {
    onSolve();
}

void MainWindow::onStartTransient() {
    createMesh();

    transient_solver_ = std::make_unique<fem::TransientSolver>(*mesh_);

    switch (material_combo_->currentIndex()) {
    case 0: transient_solver_->set_material(0, fem::Material::aluminum()); break;
    case 1: transient_solver_->set_material(0, fem::Material::steel()); break;
    case 2: transient_solver_->set_material(0, fem::Material::copper()); break;
    case 3: transient_solver_->set_material(0, fem::Material::aluminum_nonlinear()); break;
    }

    transient_solver_->set_time_step(dt_spinbox_->value());
    transient_solver_->set_end_time(end_time_spinbox_->value());
    transient_solver_->set_initial_temperature(initial_temp_spinbox_->value());
    transient_solver_->enable_openmp(parallel_checkbox_->isChecked());

    switch (theta_combo_->currentIndex()) {
    case 0: transient_solver_->set_theta(0.0); break;
    case 1: transient_solver_->set_theta(0.5); break;
    case 2: transient_solver_->set_theta(1.0); break;
    }

    transient_solver_->solve();
    temperature_history_ = transient_solver_->temperature_history();
    current_animation_step_ = 0;

    transient_running_ = true;
    start_transient_button_->setEnabled(false);
    stop_transient_button_->setEnabled(true);

    int interval = 200 - animation_speed_slider_->value() * 2;
    animation_timer_->start(std::max(10, interval));
}

void MainWindow::onStopTransient() {
    animation_timer_->stop();
    transient_running_ = false;
    start_transient_button_->setEnabled(true);
    stop_transient_button_->setEnabled(false);
}

void MainWindow::onResetTransient() {
    onStopTransient();
    current_animation_step_ = 0;
    progress_bar_->setValue(0);
    time_label_->setText("Time: 0.000 s");

    if (!temperature_history_.empty()) {
        mesh_widget_->setTemperatures(temperature_history_[0]);
    }
}

void MainWindow::onAnimationStep() {
    if (current_animation_step_ >= temperature_history_.size()) {
        onStopTransient();
        return;
    }

    mesh_widget_->setTemperatures(temperature_history_[current_animation_step_]);

    double time = current_animation_step_ * dt_spinbox_->value();
    time_label_->setText(QString("Time: %1 s").arg(time, 0, 'f', 4));

    int progress = static_cast<int>(100.0 * current_animation_step_ / temperature_history_.size());
    progress_bar_->setValue(progress);

    current_animation_step_++;
}

void MainWindow::onColorMapChanged(int index) {
    mesh_widget_->setColorMap(static_cast<fem::ColorMap>(index));
}

void MainWindow::onShowMeshToggled(bool checked) {
    mesh_widget_->setShowMesh(checked);
}

void MainWindow::onShowNodesToggled(bool checked) {
    mesh_widget_->setShowNodes(checked);
}

void MainWindow::onRunBenchmarks() {
    benchmark_log_->clear();
    appendLog("Running scaling benchmarks...\n");

    std::vector<size_t> resolutions = { 10, 20, 50, 100, 150, 200 };

    benchmark_results_ = fem::Benchmark::run_scaling_test(resolutions, true);

    appendLog(QString("%-15s %-10s %-12s %-12s %-12s %-10s")
        .arg("Name").arg("Nodes").arg("Elements").arg("Asm(ms)").arg("Solve(ms)").arg("Total(ms)"));
    appendLog(QString(80, '-'));

    for (const auto& r : benchmark_results_) {
        appendLog(QString("%-15s %-10d %-12d %-12.2f %-12.2f %-10.2f")
            .arg(QString::fromStdString(r.name))
            .arg(r.num_nodes)
            .arg(r.num_elements)
            .arg(r.assembly_time_ms)
            .arg(r.solve_time_ms)
            .arg(r.total_time_ms));
    }

    appendLog("\nSpeedups:");
    for (size_t i = 0; i < benchmark_results_.size(); i += 2) {
        if (i + 1 < benchmark_results_.size()) {
            double speedup = fem::Benchmark::compute_speedup(
                benchmark_results_[i], benchmark_results_[i + 1]);
            appendLog(QString("  %1 nodes: %.2fx")
                .arg(benchmark_results_[i].num_nodes)
                .arg(speedup));
        }
    }

    export_results_button_->setEnabled(true);
}

void MainWindow::onExportResults() {
    QString filename = QFileDialog::getSaveFileName(this, "Export Results", "", "CSV Files (*.csv)");
    if (!filename.isEmpty()) {
        fem::Benchmark::export_to_csv(benchmark_results_, filename.toStdString());
        QMessageBox::information(this, "Export", "Results exported successfully!");
    }
}

void MainWindow::onBCTypeChanged(int) {
}

void MainWindow::onNodeClicked(size_t nodeIndex, double x, double y) {
    QString msg = QString("Node %1 at (%2, %3)")
        .arg(nodeIndex)
        .arg(x, 0, 'f', 3)
        .arg(y, 0, 'f', 3);

    if (nodeIndex < static_cast<size_t>(solver_->temperatures().size())) {
        msg += QString(" - T = %1 °C").arg(solver_->temperatures()(nodeIndex), 0, 'f', 2);
    }

    statusBar()->showMessage(msg, 3000);
}

void MainWindow::updateStatusBar() {
    const auto& stats = solver_->stats();
    QString msg = QString("Nodes: %1 | Elements: %2 | Assembly: %3 ms | Solve: %4 ms | Total: %5 ms")
        .arg(mesh_->num_nodes())
        .arg(mesh_->num_elements())
        .arg(stats.assembly_time_ms, 0, 'f', 2)
        .arg(stats.solve_time_ms, 0, 'f', 2)
        .arg(stats.total_time_ms, 0, 'f', 2);

    status_label_->setText(msg);
}

void MainWindow::appendLog(const QString& message) {
    benchmark_log_->append(message);
}