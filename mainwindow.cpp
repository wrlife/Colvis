#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    m_filemanager= new Filemanager;
    m_syndata = new Genesyndata;

    QVTKWidget* widget = this->findChild<QVTKWidget*>("qvtk");

    widget->GetRenderWindow()->AddRenderer( m_syndata->getrenderer() );

}

MainWindow::~MainWindow()
{
    delete ui;
    delete m_filemanager;
    delete m_syndata;
}


//==================
//Load new model
//==================
void MainWindow::on_action_New_file_triggered()
{
    //Open file dialog
    QString filePath = QFileDialog::getOpenFileName(
                       this, tr("Open File"), "",
                       tr("3Dmodels (*.stl *.vtp *.off)"));


    if (filePath.isEmpty()){

        return;
    }

    //Read file using filemanager
    m_filemanager->loadnewfile(filePath);

    //Render file
    m_syndata->rendermodel(m_filemanager->getfile());
    QVTKWidget* widget = this->findChild<QVTKWidget*>("qvtk");
    widget->GetRenderWindow()->Render();

}


//==================
//Load camera path
//==================
void MainWindow::on_action_Load_camera_path_triggered()
{
    //Open file dialog
    QString filePath = QFileDialog::getOpenFileName(
                       this, tr("Open File"), "",
                       tr("3Dmodels (*.stl *.vtp *.off)"));


    if (filePath.isEmpty()){

        return;
    }

    //Read file using filemanager
    m_filemanager->loadnewcamera(filePath);

    //Render file
    m_syndata->loadcamerapath(m_filemanager->getfile());

    m_timer=new QTimer(this);

    connect(m_timer,SIGNAL(timeout()),this,SLOT(updatecamera()));

    m_timer->start(300);


}



//==================
//Camera animation
//==================
void MainWindow::updatecamera()
{
    m_syndata->updatecamera();
    QVTKWidget* widget = this->findChild<QVTKWidget*>("qvtk");
    widget->GetRenderWindow()->Render();

}

