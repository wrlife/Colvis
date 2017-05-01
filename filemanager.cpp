#include "filemanager.h"

Filemanager::Filemanager()
{
    m_vtp = vtkSmartPointer<vtkXMLPolyDataReader>::New();

    m_stl=vtkSmartPointer<vtkSTLReader>::New();

    m_camera=vtkSmartPointer<vtkXMLPolyDataReader>::New();

    m_ply=vtkSmartPointer<vtkPLYReader>::New();
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
    else if(filename.endsWith("ply")){
        m_ply->SetFileName(filename.toStdString().c_str());
        m_ply->Update();
        m_polydata=m_ply->GetOutput();
    }
}



void Filemanager::loadnewcamera(QString filename)
{
    m_camera->SetFileName(filename.toStdString().c_str());
    m_camera->Update();
    m_polydata=m_camera->GetOutput();
}


//==================
//Render 3D parametric model
//==================
void Filemanager::renderparametricmodel()
{

    parametricObjects = vtkSmartPointer<vtkParametricBoy>::New();
    parametricFunctionSources= vtkSmartPointer<vtkParametricFunctionSource>::New();

    parametricFunctionSources->SetParametricFunction(parametricObjects);
    parametricFunctionSources->SetUResolution(101);
    parametricFunctionSources->SetVResolution(101);
    parametricFunctionSources->SetWResolution(101);
    parametricFunctionSources->Update();

    m_polydata = parametricFunctionSources->GetOutput();

}

Filemanager::~Filemanager()
{

}
