#-------------------------------------------------
#
# Project created by QtCreator 2015-03-31T15:50:41
#
#-------------------------------------------------

QT       += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = MeshDenoising
TEMPLATE = app

DEFINES += _USE_MATH_DEFINES

SOURCES += main.cpp\
        mainwindow.cpp \
    glviewer.cpp \
    glexaminer.cpp \
    meshexaminer.cpp \
    datamanager.cpp \
    parameterset.cpp \
    parametersetwidget.cpp \
    calculationthread.cpp \
    Algorithms/Noise.cpp \
    Algorithms/MeshDenoisingBase.cpp \
    iothread.cpp \
    Algorithms/MeshDenoisingViaL0MinimizationForFrames.cpp \
    Algorithms/MeshDenoisingViaL0Minimization.cpp

HEADERS  += mainwindow.h \
    glviewer.h \
    glwrapper.h \
    glexaminer.h \
    meshexaminer.h \
    mesh.h \
    datamanager.h \
    parameterset.h \
    parametersetwidget.h \
    calculationthread.h \
    Algorithms/Noise.h \
    Algorithms/MeshDenoisingBase.h \
    iothread.h \
    Algorithms/MeshDenoisingViaL0MinimizationForFrames.h \
    Algorithms/MeshDenoisingViaL0Minimization.h

FORMS    += mainwindow.ui

include(../LibsInclude.pri)

RESOURCES += \
    mainwindow.qrc

INCLUDEPATH += $$PWD/../ANN/include/
#LIBS += E:/reconstruction3D/denosing/GuidedDenoising-master/src/ANN/lib/ANN.lib
#LIBS += -L../ANN/lib -lANN
#LIBS += -L$$PWD/../ANN/lib/ -lANN

#DESTDIR = $$PWD/bin
#OBJECTS_DIR = $$PWD/obj
#MOC_DIR = $$PWD/moc
