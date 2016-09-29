#ifndef LIGHTDIALOG_H
#define LIGHTDIALOG_H

#include <QDialog>
#include "filemanager.h"


class Genesyndata;

namespace Ui {
class lightDialog;
}

class lightDialog : public QDialog
{
    Q_OBJECT

public:
    explicit lightDialog(QWidget *parent = 0);

    void setscene(Genesyndata* t_scene);

    ~lightDialog();

private slots:


    void on_radioButton_5_clicked();

    void on_radioButton_3_clicked();

    void on_radioButton_4_clicked();

private:
    Ui::lightDialog *ui;
    Genesyndata* m_scene;

};

#endif // LIGHTDIALOG_H
