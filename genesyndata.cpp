#include "genesyndata.h"
#include <iostream>
#include <vtkCenterOfMass.h>


Genesyndata::Genesyndata()
{

    //Mapper
    m_Mapper =
        vtkSmartPointer<vtkPolyDataMapper>::New();

    //Actor
    m_Actor =
        vtkSmartPointer<vtkActor>::New();
    m_Actor->SetMapper(m_Mapper);

    //light

    //Renderer
    m_renderer =
        vtkSmartPointer<vtkRenderer>::New();
    m_renderer->AddActor(m_Actor);

    counter=0;

}


//==================
//Render 3D model
//==================
void Genesyndata::rendermodel(vtkSmartPointer<vtkPolyData> t_model)
{

    m_Mapper->SetInputData(t_model);

    m_renderer->ResetCamera();

}


//==================
//Store camera path
//==================
void Genesyndata::loadcamerapath(vtkSmartPointer<vtkPolyData> t_model)
{

    //m_renderer->ResetCamera();


    m_camerapath=t_model;



}


//==================
//Set current camera position and look at
//==================
void Genesyndata::updatecamera()
{
    //Get inital camera position
    double p[3];
    m_camerapath->GetPoint(counter,p);


    vtkSmartPointer<vtkCamera> t_camera=m_renderer->GetActiveCamera();
    t_camera->SetPosition(p[0],p[1],p[2]);

    m_camerapath->GetPoint(counter+1,p);

    t_camera->SetFocalPoint(p[0],p[1],p[2]);

    counter++;
}


Genesyndata::~Genesyndata(){}

