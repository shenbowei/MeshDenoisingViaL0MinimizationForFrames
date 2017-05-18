#include "calculationthread.h"

CalculationThread::CalculationThread()
{
    noise_ = NULL;
    mesh_denoise_base_ = NULL;
    algorithm_name_ = "";
}

CalculationThread::~CalculationThread()
{
    if(noise_ != NULL)
        delete noise_;
    noise_ = NULL;

    if(mesh_denoise_base_ != NULL)
        delete mesh_denoise_base_;
    mesh_denoise_base_ = NULL;
}

void CalculationThread::initAlgorithm(DataManager *_data_manager, ParameterSet *_parameter_set)
{
    switch (algorithms_type_) {
    case kNoise:
        noise_ = new Noise(_data_manager, _parameter_set);
        break;
    case kMeshDenoisingViaL0Minimization:
        mesh_denoise_base_ = new MeshDenoisingViaL0Minimization(_data_manager, _parameter_set);
        break;
    case kMeshDenoisingViaL0MinimizationForFrames:
        mesh_denoise_base_ = new MeshDenoisingViaL0MinimizationForFrames(_data_manager, _parameter_set);
        break;
    default:
        break;
    }

    connect(mesh_denoise_base_, SIGNAL(statusMessage(QString)), this, SLOT(sendMessage(QString)));
}

void CalculationThread::run()
{
    emit(statusShowMessage("Now applying algorithm --" + algorithm_name_ + "-- ..."));
    emit(setActionAndWidget(false, false));
    if(algorithms_type_ == kNoise)
        noise_->addNoise();
    else
        mesh_denoise_base_->denoise();
    emit(setActionAndWidget(true, false));
    //emit(statusShowMessage("Applying algorithm --" + algorithm_name_ + "-- done."));

    emit(needToUpdateGL(false));
}

void CalculationThread::sendMessage(QString message)
{
    emit(statusShowMessage(message));
}
