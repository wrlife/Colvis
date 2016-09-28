#include "filemanager.h"

Filemanager::Filemanager()
{
    m_vtp = vtkSmartPointer<vtkXMLPolyDataReader>::New();

    m_stl=vtkSmartPointer<vtkSTLReader>::New();

    m_camera=vtkSmartPointer<vtkXMLPolyDataReader>::New();
}


void Filemanager::loadnewfile(QString filename)
{


    if (filename.endsWith("vtp")){
        m_vtp->SetFileName(filename.toStdString().c_str());
        m_vtp->Update();
        m_polydata=m_vtp->GetOutput();
    }
    else if(filename.endsWith("stl")){
        m_stl->SetFileName(filename.toStdString().c_str());
        m_stl->Update();
        m_polydata=m_stl->GetOutput();
    }
}


void Filemanager::loadnewcamera(QString filename)
{
    m_camera->SetFileName(filename.toStdString().c_str());
    m_camera->Update();
    m_polydata=m_camera->GetOutput();
}

Filemanager::~Filemanager()
{

}