#ifndef CALCULATIONTHREAD_H
#define CALCULATIONTHREAD_H

#include <QThread>
#include "Algorithms/Noise.h"
#include "Algorithms/MeshDenoisingBase.h"
#include "Algorithms/MeshDenoisingViaL0MinimizationForFrames.h"
#include "Algorithms/MeshDenoisingViaL0Minimization.h"


class CalculationThread : public QThread
{
    Q_OBJECT
public:
    CalculationThread();
    ~CalculationThread();

    enum AlgorithmsType{kNon, kNoise,kMeshDenoisingViaL0Minimization, kMeshDenoisingViaL0MinimizationForFrames};

    void initAlgorithm(DataManager *_data_manager, ParameterSet *_parameter_set);
    void setAlgorithmName(const QString algorithm_name) {algorithm_name_ = algorithm_name;}

    void run();

signals:
    void needToUpdateGL(bool value);
    void setActionAndWidget(bool, bool);
    void statusShowMessage(QString message);

public slots:
    void sendMessage(QString message);

public:
    QString algorithm_name_;
    AlgorithmsType algorithms_type_;
    Noise *noise_;
    MeshDenoisingBase *mesh_denoise_base_;
};

#endif // CALCULATIONTHREAD_H
