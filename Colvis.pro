#-------------------------------------------------
#
# Project created by QtCreator 2016-09-27T10:43:16
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = Colvis
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    genesyndata.cpp \
    filemanager.cpp \
    lightdialog.cpp

HEADERS  += mainwindow.h \
    genesyndata.h \
    filemanager.h \
    lightdialog.h

QMAKE_CXXFLAGS += -Wno-deprecated
INCLUDEPATH += /usr/local/include/vtk-7.0

LIBS += -L/usr/local/lib -ldl -lGL -lvtkCommonColor-7.0 -lvtkCommonCore-7.0 \
  -lvtkCommonComputationalGeometry-7.0 -lvtkCommonDataModel-7.0 -lvtkCommonExecutionModel-7.0 \
  -lvtkCommonMath-7.0 -lvtkCommonMisc-7.0 -lvtkCommonSystem-7.0 -lvtkCommonTransforms-7.0 -lvtkFiltersCore-7.0 \
  -lvtkFiltersExtraction-7.0 -lvtkFiltersGeneral-7.0 -lvtkFiltersGeometry-7.0 -lvtkFiltersHybrid-7.0 \
  -lvtkFiltersSources-7.0 -lvtkFiltersStatistics-7.0 -lvtkImagingCore-7.0 -lvtkImagingColor-7.0 \
  -lvtkImagingFourier-7.0 -lvtkImagingMath-7.0 -lvtkImagingMorphological-7.0 -lvtkImagingStencil-7.0 \
  -lvtkIOCore-7.0 -lvtkIOGeometry-7.0 -lvtkIOImage-7.0 -lvtkIOPLY-7.0 -lvtkInteractionStyle-7.0 \
  -lvtkInteractionWidgets-7.0 -lvtkGUISupportQt-7.0 -lvtkRenderingCore-7.0 -lvtkRenderingOpenGL2-7.0 \
  -lvtkRenderingQt-7.0 -lvtkRenderingVolume-7.0 -lvtkalglib-7.0 -lvtkjsoncpp-7.0 -lvtksys-7.0 -lvtkzlib-7.0 \
  -lvtkIOXML-7.0


FORMS    += mainwindow.ui \
    lightdialog.ui
