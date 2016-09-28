#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QFileDialog>
#include "filemanager.h"
#include "genesyndata.h"
#include <QTimer>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:

    void on_action_New_file_triggered();

    void on_action_Load_camera_path_triggered();

    void updatecamera();

private:
    Ui::MainWindow *ui;
    Filemanager * m_filemanager;
    //Start Generate synthetic data
    Genesyndata* m_syndata;
    QTimer *m_timer;
};

#endif // MAINWINDOW_H
