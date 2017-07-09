#include "mainwindow.h"
#include <QApplication>
#include <glut>

int main(int argc, char *argv[])
{
    QApplication::setStyle(QStyleFactory::create("cleanlooks"));
    QApplication a(argc, argv);
    MainWindow w;
    w.show();

    return a.exec();
}
