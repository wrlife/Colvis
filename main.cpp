#include <vtkAutoInit.h>
#include "mainwindow.h"
#include <QApplication>


#include <QVTKWidget.h>



VTK_MODULE_INIT(vtkRenderingOpenGL2);

int main(int argc, char *argv[])
{

    QApplication a(argc, argv);
    MainWindow w;

    w.resize(800,600);
    w.show();
    w.addlight();





    return a.exec();
}
