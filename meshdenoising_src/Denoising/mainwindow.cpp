#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "glviewer.h"
#include "datamanager.h"
#include "parameterset.h"
#include "parametersetwidget.h"
#include "calculationthread.h"
#include "iothread.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    setWindowTitle("MeshDenoising");
    setWindowIcon(QIcon(":Icons/MeshDenoise.ico"));

    setGeometry(250, 50, 880, 660);

    init();
    CreateActions();
    CreateMenus();
    CreateToolBars();
    CreateStatusBar();

    QVBoxLayout *layout_left_ = new QVBoxLayout;
    layout_left_->addStretch(2);

    QHBoxLayout *layout_main = new QHBoxLayout;
    layout_main->addLayout(layout_left_);
    layout_main->addWidget(opengl_viewer_);
    layout_main->setStretch(1,1);
    this->centralWidget()->setLayout(layout_main);
}

MainWindow::~MainWindow()
{
    CloseWidget();

    if(calculation_thread_ != NULL)
        delete calculation_thread_;
    calculation_thread_ = NULL;

    if(opengl_viewer_ != NULL)
        delete opengl_viewer_;
    opengl_viewer_ = NULL;

    if(data_manager_ != NULL)
        delete data_manager_;
    data_manager_ = NULL;

    if(parameter_set_ != NULL)
        delete parameter_set_;
    parameter_set_ = NULL;
}

void MainWindow::init()
{
    opengl_viewer_ = new GLViewer(this);
    data_manager_ = new DataManager();
    parameter_set_ = new ParameterSet();
    parameter_set_widget_ = NULL;
    io_thread_ = new ioThread(data_manager_);
    io_thread_->io_type_ = ioThread::kNon;
    connect(io_thread_, SIGNAL(needToUpdateGL(bool)), this, SLOT(needToUpdateGL(bool)));
    connect(io_thread_, SIGNAL(setActionAndWidget(bool, bool)), this, SLOT(setActionAndWidget(bool, bool)));
    calculation_thread_ = new CalculationThread();
    calculation_thread_->algorithms_type_ = CalculationThread::kNon;
    connect(calculation_thread_, SIGNAL(needToUpdateGL(bool)), this, SLOT(needToUpdateGL(bool)));
    connect(calculation_thread_, SIGNAL(setActionAndWidget(bool,bool)), this, SLOT(setActionAndWidget(bool,bool)));
}

void MainWindow::CreateActions()
{
    action_import_mesh_ = new QAction(QIcon(":/Icons/open.ico"),tr("Import Mesh"),this);
    action_import_mesh_->setStatusTip("Import Mesh.");
    connect(action_import_mesh_, SIGNAL(triggered()), this, SLOT(ImportMesh()));

    action_export_mesh_ = new QAction(QIcon(":/Icons/save.ico"),tr("Export Mesh"), this);
    action_export_mesh_->setStatusTip("Export Mesh.");
    action_export_mesh_->setEnabled(false);
    connect(action_export_mesh_, SIGNAL(triggered()), this, SLOT(ExportMesh()));

    action_exit_ = new QAction(tr("Exit"), this);
    action_exit_->setStatusTip("Exit.");
    connect(action_exit_, SIGNAL(triggered()), this, SLOT(close()));

    action_render_points_ = new QAction(QIcon(":/Icons/points.ico"), tr("Render Points"), this);
    action_render_points_->setStatusTip("Render Points.");
    action_render_points_->setCheckable(true);
    connect(action_render_points_, SIGNAL(toggled(bool)), opengl_viewer_, SLOT(setDrawPointsStatus(bool)));

    action_render_edges_ = new QAction(QIcon(":/Icons/edges.ico"), tr("Render Edges"), this);
    action_render_edges_->setStatusTip("Render Edges.");
    action_render_edges_->setCheckable(true);
    connect(action_render_edges_, SIGNAL(toggled(bool)), opengl_viewer_, SLOT(setDrawEdgesStatus(bool)));

    action_render_faces_ = new QAction(QIcon(":/Icons/faces.ico"), tr("Render Faces"), this);
    action_render_faces_->setStatusTip("Render Faces.");
    action_render_faces_->setCheckable(true);
    action_render_faces_->setChecked(true);
    connect(action_render_faces_, SIGNAL(toggled(bool)), opengl_viewer_, SLOT(setDrawFacesStatus(bool)));

    action_set_background_color_ = new QAction(QIcon(":/Icons/background.ico"),tr("Change Background Color"), this);
    action_set_background_color_->setStatusTip("Change Background Color.");
    connect(action_set_background_color_, SIGNAL(triggered()), opengl_viewer_, SLOT(setBackgroundColor()));

    action_set_model_color_ = new QAction(QIcon(":/Icons/model_color.ico"), tr("Change Model Color"), this);
    action_set_model_color_->setStatusTip("Change Model Color.");
    action_set_model_color_->setCheckable(true);
    //action_set_model_color_->setChecked(true);
    connect(action_set_model_color_, SIGNAL(toggled(bool)), opengl_viewer_, SLOT(setDrawModelColor(bool)));

    action_to_original_mesh_ = new QAction(this);
    action_to_original_mesh_->setText("Original Mesh");
    action_to_original_mesh_->setStatusTip("Show Original Mesh.");
    connect(action_to_original_mesh_, SIGNAL(triggered()), this, SLOT(TransToOriginalMesh()));

    action_to_noisy_mesh_ = new QAction(this);
    action_to_noisy_mesh_->setText("Noisy Mesh");
    action_to_noisy_mesh_->setStatusTip("Show Noisy Mesh.");
    connect(action_to_noisy_mesh_, SIGNAL(triggered()), this, SLOT(TransToNoisyMesh()));

    action_to_denoised_mesh_ = new QAction(this);
    action_to_denoised_mesh_->setText("Denoised Mesh");
    action_to_denoised_mesh_->setStatusTip("Show Denoised Mesh.");
    connect(action_to_denoised_mesh_, SIGNAL(triggered()), this, SLOT(TransToDenoisedMesh()));

    action_clear_mesh_ = new QAction(this);
    action_clear_mesh_->setText("Clear");
    action_clear_mesh_->setStatusTip("Clear Mesh.");
    connect(action_clear_mesh_, SIGNAL(triggered()), this, SLOT(ClearMesh()));

    action_noise_ = new QAction(tr("Noise"), this);
    action_noise_->setStatusTip("Add Noise to Mesh using Algorithm -- Noise.");
    connect(action_noise_, SIGNAL(triggered()), this, SLOT(ShowNoiseWidget()));

    action_mesh_denoising_via_l0_minimization_ = new QAction(tr("Mesh Denoising via L0 Minimization"), this);
    action_mesh_denoising_via_l0_minimization_->setStatusTip("Denoise Mesh using Algorithm -- Mesh Denoising via L0 Minimization.");
    connect(action_mesh_denoising_via_l0_minimization_, SIGNAL(triggered()), this, SLOT(ShowMeshDenoisingViaL0MinimizationWidget()));

    action_mesh_denoising_via_l0_minimization_for_frames_ = new QAction(tr("Mesh Denoising via L0 Minimization For Frames"), this);
    action_mesh_denoising_via_l0_minimization_for_frames_->setStatusTip("Denoise Mesh using Algorithm -- Mesh Denoising via L0 Minimization For Frames.");
    connect(action_mesh_denoising_via_l0_minimization_for_frames_, SIGNAL(triggered()), this, SLOT(ShowMeshDenoisingViaL0MinimizationForFramesWidget()));

    action_about_ = new QAction(QIcon(":/Icons/about.ico"),tr("About"), this);
    action_about_->setStatusTip("Information about this UI.");
    connect(action_about_, SIGNAL(triggered()), this, SLOT(About()));
}

void MainWindow::CreateMenus()
{
    menu_file_ = menuBar()->addMenu(tr("File"));
    menu_file_->addAction(action_import_mesh_);
    menu_file_->addAction(action_export_mesh_);
    menu_file_->addSeparator();
    menu_file_->addAction(action_exit_);

    menu_algorithms_ = menuBar()->addMenu(tr("Algorithms"));
    menu_algorithms_->addAction(action_noise_);
    menu_algorithms_->addSeparator();
    menu_algorithms_->addAction(action_mesh_denoising_via_l0_minimization_);
    menu_algorithms_->addAction(action_mesh_denoising_via_l0_minimization_for_frames_);
    menu_algorithms_->setEnabled(false);
    menu_help_ = menuBar()->addMenu(tr("Help"));
    menu_help_->addAction(action_about_);
}

void MainWindow::CreateToolBars()
{
    toolbar_file_ = addToolBar(tr("File"));
    toolbar_file_->addAction(action_import_mesh_);
    toolbar_file_->addAction(action_export_mesh_);

    toolbar_opengl_info_ = addToolBar(tr("OpenGL"));
    toolbar_opengl_info_->addAction(action_render_points_);
    toolbar_opengl_info_->addAction(action_render_edges_);
    toolbar_opengl_info_->addAction(action_render_faces_);
    toolbar_opengl_info_->addSeparator();
    toolbar_opengl_info_->addAction(action_set_background_color_);
    toolbar_opengl_info_->addAction(action_set_model_color_);

    toolbar_mesh_status_ = addToolBar(tr("Status"));
    toolbar_mesh_status_->addAction(action_to_original_mesh_);
    toolbar_mesh_status_->addAction(action_to_noisy_mesh_);
    toolbar_mesh_status_->addAction(action_to_denoised_mesh_);
    toolbar_mesh_status_->addSeparator();
    toolbar_mesh_status_->addAction(action_clear_mesh_);
    toolbar_mesh_status_->setEnabled(false);
}

void MainWindow::CreateStatusBar()
{
    statusBar()->setStyleSheet(QString("QStatusBar::item{border: 0px}"));
    label_operation_info_ = new QLabel();
    label_operation_info_->setAlignment(Qt::AlignCenter);
    label_operation_info_->setMinimumSize(label_operation_info_->sizeHint());

    statusBar()->addWidget(label_operation_info_);
    connect(io_thread_, SIGNAL(statusShowMessage(QString)), label_operation_info_, SLOT(setText(QString)));
    connect(calculation_thread_, SIGNAL(statusShowMessage(QString)), label_operation_info_, SLOT(setText(QString)));
}

void MainWindow::CloseWidget()
{
    if(parameter_set_widget_ != NULL)
        delete parameter_set_widget_;
    parameter_set_widget_ = NULL;
}

void MainWindow::ShowWidget()
{
    calculation_thread_->initAlgorithm(data_manager_, parameter_set_);
    parameter_set_widget_ = new ParameterSetWidget(this, parameter_set_);

    connect(parameter_set_widget_, SIGNAL(ReadyToApply(QString)), this, SLOT(ApplyAlgorithms(QString)));
    addDockWidget(Qt::RightDockWidgetArea, parameter_set_widget_);
    parameter_set_widget_->showWidget();
}

void MainWindow::SetActionStatus(bool value)
{
    action_import_mesh_->setEnabled(value);
    action_export_mesh_->setEnabled(value);

    menu_algorithms_->setEnabled(value);

    toolbar_mesh_status_->setEnabled(value);
}

void MainWindow::ImportMesh()
{
    QString filename = QFileDialog::getOpenFileName(this, tr("Import Mesh"), "../Models", tr(" *.obj *.off *.ply"));

    if(filename.isNull())
    {
        return;
    }

    io_thread_->setFileName(filename);
    io_thread_->io_type_ = ioThread::kImport;
    io_thread_->start();
}

void MainWindow::ExportMesh()
{
    QString filename = QFileDialog::getSaveFileName(this, tr("Export Mesh"), ".", tr(" OBJ File (*.obj);;OFF File (*.off);;PLY File (*.ply)"));

    if(filename.isNull())
    {
        return;
    }

    io_thread_->setFileName(filename);
    io_thread_->io_type_ = ioThread::kExport;
    io_thread_->start();
}

void MainWindow::TransToNoisyMesh()
{
    data_manager_->MeshToNoisyMesh();
    opengl_viewer_->resetMesh(data_manager_->getMesh());
    opengl_viewer_->updateGL();
}

void MainWindow::TransToOriginalMesh()
{
    data_manager_->MeshToOriginalMesh();
    opengl_viewer_->resetMesh(data_manager_->getMesh());
    opengl_viewer_->updateGL();
}

void MainWindow::TransToDenoisedMesh()
{
    data_manager_->MeshToDenoisedMesh();
    opengl_viewer_->resetMesh(data_manager_->getMesh());
    opengl_viewer_->updateGL();
}

void MainWindow::ClearMesh()
{
    CloseWidget();
    data_manager_->ClearMesh();
    opengl_viewer_->resetMesh(data_manager_->getMesh());
    opengl_viewer_->updateGL();
    SetActionStatus(false);
    action_import_mesh_->setEnabled(true);
}

void MainWindow::ApplyAlgorithms(QString algorithm_name)
{
    calculation_thread_->setAlgorithmName(algorithm_name);
    calculation_thread_->start();
}

void MainWindow::ShowNoiseWidget()
{
    calculation_thread_->algorithms_type_ = CalculationThread::kNoise;
    CloseWidget();
    ShowWidget();
}

void MainWindow::ShowMeshDenoisingViaL0MinimizationWidget()
{
    calculation_thread_->algorithms_type_ = CalculationThread::kMeshDenoisingViaL0Minimization;
    CloseWidget();
    ShowWidget();
}

void MainWindow::ShowMeshDenoisingViaL0MinimizationForFramesWidget()
{
    calculation_thread_->algorithms_type_ = CalculationThread::kMeshDenoisingViaL0MinimizationForFrames;
    CloseWidget();
    ShowWidget();
}

void MainWindow::setActionAndWidget(bool value1, bool value2)
{
    SetActionStatus(value1);

    if(parameter_set_widget_ != NULL && value2){
        CloseWidget();
        ShowWidget();
    }
}

void MainWindow::needToUpdateGL(bool value)
{
    opengl_viewer_->resetMesh(data_manager_->getMesh(), value);
    opengl_viewer_->updateGL();
}

void MainWindow::About()
{
    QMessageBox::about(this, tr("About MeshDenoise UI"),
                tr("MeshDenoisingViaL0MinimizationForFrames<br>"
                   "The program is based on Qt, OpenMesh, ANN and Eigen Library. <br>"
                   "If you have any question, please contact <br>"
                   "<a href='mailto:shenbowei0@gmail.com'>shenbowei0@gmail.com</a>."));
}
