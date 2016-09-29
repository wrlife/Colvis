#include "lightdialog.h"
#include "genesyndata.h"
#include "ui_lightdialog.h"

lightDialog::lightDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::lightDialog)
{
    ui->setupUi(this);
}



void lightDialog::setscene(Genesyndata *t_scene)
{
    m_scene=t_scene;
}


lightDialog::~lightDialog()
{
    delete ui;
}



void lightDialog::on_radioButton_5_clicked()
{
    m_scene->setconstantlight(this->findChild<QSlider*>("Attenvalue")->value());
}

void lightDialog::on_radioButton_3_clicked()
{
    m_scene->setlinearlight(this->findChild<QSlider*>("Attenvalue")->value());
}

void lightDialog::on_radioButton_4_clicked()
{
    m_scene->setqudraticlight(this->findChild<QSlider*>("Attenvalue")->value());
}
