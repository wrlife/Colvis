#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <vtkCallbackCommand.h>
#include <vtkCommand.h>




void KeypressCallbackFunction (
  vtkObject* caller,
  long unsigned int eventId,
  void* clientData,
  void* callData );


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    m_filemanager= new Filemanager;
    m_syndata = new Genesyndata;

    QVTKWidget* widget = this->findChild<QVTKWidget*>("qvtk");

    widget->GetRenderWindow()->AddRenderer( m_syndata->getrenderer() );
    m_lightdialog.setscene(m_syndata);
    m_lightdialog.setwindow(this);

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
                       tr("3Dmodels (*.stl *.vtp *.off *.ply)"));


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

    this->updatecamera(0);

//    m_timer=new QTimer(this);

//    connect(m_timer,SIGNAL(timeout()),this,SLOT(updatecamera()));

//    m_timer->start(300);


}



//==================
//Camera animation
//==================
void MainWindow::updatecamera(int c_step)
{
    m_syndata->updatecamera(c_step);
    QVTKWidget* widget = this->findChild<QVTKWidget*>("qvtk");
    widget->GetRenderWindow()->Render();

}


void MainWindow::addlight()
{
    QVTKWidget* widget = this->findChild<QVTKWidget*>("qvtk");
    widget->GetRenderWindow()->Render();

    m_syndata->getrenderer()->RemoveAllLights();


    m_syndata->addlight();
}


void MainWindow::on_action_Lighting_triggered()
{
    m_lightdialog.exec();
}



void MainWindow::on_actionModify_triggered()
{
    //m_syndata->modiflight();
    QVTKWidget* widget = this->findChild<QVTKWidget*>("qvtk");

//    int stepsize=20;

//    for (int i=0;i<m_syndata->get_num_cams()*stepsize-2;i++)
//    {
//        this->updatecamera(i);
//        m_syndata->get_z_values(widget->GetRenderWindow());
//        //this->movecamaround();
//    }

    m_syndata->get_orthognal_normal_view(m_filemanager->getfile(),widget->GetRenderWindow());
    //m_syndata->get_orthognal_normal_view(m_filemanager->getfile(),widget->GetRenderWindow());
}


void MainWindow::movecamaround()
{
    QVTKWidget* widget = this->findChild<QVTKWidget*>("qvtk");

    m_syndata->get_z_values(widget->GetRenderWindow());
    srand((unsigned)time(NULL));
    float x,y,z;
    float elevation,azimuth;
    for (int i=0;i<50;i++)
    {
         x= static_cast <float> (rand()) / static_cast <float> (RAND_MAX)*2-1;
         y= static_cast <float> (rand()) / static_cast <float> (RAND_MAX)*2-1;
         z= static_cast <float> (rand()) / static_cast <float> (RAND_MAX)*2-1;m_syndata->get_z_values(widget->GetRenderWindow());

         elevation= static_cast <float> (rand()) / static_cast <float> (RAND_MAX)*360;
         azimuth= static_cast <float> (rand()) / static_cast <float> (RAND_MAX)*360;

         m_syndata->get_z_values(widget->GetRenderWindow());
         m_syndata->randomcampos(x,y,z,elevation,azimuth);
         widget->GetRenderWindow()->Render();

    }
}

void MainWindow::on_actionParametricBoy_triggered()
{
    m_filemanager->renderparametricmodel();

    m_syndata->rendermodel(m_filemanager->getfile());



}


void KeypressCallbackFunction ( vtkObject* caller, long unsigned int vtkNotUsed(eventId), void* clientData, void* vtkNotUsed(callData) )
{

    MainWindow* t_window =
      static_cast<MainWindow*>(clientData);
    QVTKWidget* widget = t_window->findChild<QVTKWidget*>("qvtk");
    t_window->m_syndata->get_orthognal_normal_view(t_window->m_filemanager->getfile(),widget->GetRenderWindow());

    //std::cout<<"call back"<<std::endl;

}


void MainWindow::on_actionAnalyse_contour_triggered()
{
    vtkSmartPointer<vtkCallbackCommand> keypressCallback =
      vtkSmartPointer<vtkCallbackCommand>::New();
    keypressCallback->SetCallback ( KeypressCallbackFunction );
    keypressCallback->SetClientData(this);

    QVTKWidget* widget = this->findChild<QVTKWidget*>("qvtk");

    widget->GetRenderWindow()->GetInteractor()->AddObserver ( vtkCommand::KeyPressEvent, keypressCallback );

}
