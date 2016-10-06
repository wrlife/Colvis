#ifndef LIGHTDIALOG_H
#define LIGHTDIALOG_H

#include <QDialog>
#include "filemanager.h"


class Genesyndata;
class MainWindow;

namespace Ui {
class lightDialog;
}

class lightDialog : public QDialog
{
    Q_OBJECT

public:
    explicit lightDialog(QWidget *parent = 0);

    void setscene(Genesyndata* t_scene);

    void setwindow(MainWindow * m_window);

    ~lightDialog();

private slots:


    void on_radioButton_5_clicked();

    void on_radioButton_3_clicked();

    void on_radioButton_4_clicked();

    void on_horizontalSlider_valueChanged(int value);

    void on_horizontalSlider_2_valueChanged(int value);

    void on_horizontalSlider_3_valueChanged(int value);

private:
    Ui::lightDialog *ui;
    Genesyndata* m_scene;
    MainWindow * m_window;

};

#endif // LIGHTDIALOG_H
