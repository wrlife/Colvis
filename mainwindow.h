#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QFileDialog>
#include "filemanager.h"
#include "genesyndata.h"
#include "lightdialog.h"
#include <QTimer>
#include <vtkLightCollection.h>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    void addlight();

    void movecamaround();


    ~MainWindow();

private slots:

    void on_action_New_file_triggered();

    void on_action_Load_camera_path_triggered();

    void updatecamera(int c_step);

    void on_action_Lighting_triggered();

    void on_actionModify_triggered();

    void on_actionParametricBoy_triggered();

private:
    Ui::MainWindow *ui;
    Filemanager * m_filemanager;
    //Start Generate synthetic data
    Genesyndata* m_syndata;
    QTimer *m_timer;

    lightDialog m_lightdialog;

    vtkLightCollection* m_originalLights;
};

#endif // MAINWINDOW_H
