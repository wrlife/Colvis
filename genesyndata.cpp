#include "genesyndata.h"
#include <iostream>
#include <vtkCenterOfMass.h>
#include <vtkProperty.h>
#include <vtkPolyDataNormals.h>
#include <vtkLightCollection.h>


Genesyndata::Genesyndata()
{

    m_sphereSource =
      vtkSmartPointer<vtkSphereSource>::New();
    //Mapper
    m_Mapper =
        vtkSmartPointer<vtkPolyDataMapper>::New();
      m_Mapper->SetInputConnection(m_sphereSource->GetOutputPort());

    //Actor
    m_Actor =
        vtkSmartPointer<vtkActor>::New();
    m_Actor->SetMapper(m_Mapper);

    m_Actor->GetProperty()->SetInterpolationToGouraud();

    //light
    double lightPosition[3] = {0, 0, 2};

    // Create a light
    double lightFocalPoint[3] = {1,0,0};
    double *p;

    m_light = vtkSmartPointer<vtkLight>::New();

    m_light->SetLightTypeToHeadlight();


    //m_light->SetPosition(lightPosition[0], lightPosition[1], lightPosition[2]);


   // m_light->SetFocalPoint(lightFocalPoint[0], lightFocalPoint[1], lightFocalPoint[2]);



    // Display where the light is
    //m_lightActor = vtkSmartPointer<vtkLightActor>::New();
    //m_lightActor->SetLight(m_light);


    //Renderer
    m_renderer =
        vtkSmartPointer<vtkRenderer>::New();
    m_light->SetAttenuationValues(0,0,0.3);

    m_light->SetPositional(true); // required for vtkLightActor below

    m_light->SetConeAngle(180);

    m_light->SetDiffuseColor(255,255,255);
    m_light->SetAmbientColor(255,255,255);
    m_light->SetSpecularColor(255,255,255);

    m_light->SetLightTypeToHeadlight();


    //p=m_renderer->GetAmbient();

    //p=m_light->GetAttenuationValues();

    //m_renderer->AddViewProp(m_lightActor);
    //originalLights->GetNextItem()->SetLightTypeToHeadlight();

    m_renderer->AddActor(m_Actor);


    counter=0;

}

//==================
//Render 3D model
//==================
void Genesyndata::rendermodel(vtkSmartPointer<vtkPolyData> t_model)
{

    // calculate normals
     vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
     normalGenerator->SetInputData(t_model);
     normalGenerator->ComputePointNormalsOn();
     normalGenerator->ComputeCellNormalsOn();
     normalGenerator->Update();

     t_model = normalGenerator->GetOutput();


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


void Genesyndata::addlight()
{
    m_renderer->AddLight(m_light);
}

void Genesyndata::setconstantlight(int value)
{
    m_light->SetAttenuationValues(float(value)/100,0,0);

}

void Genesyndata::setlinearlight(int value)
{
    m_light->SetAttenuationValues(0,float(value)/100,0);

}

void Genesyndata::setqudraticlight(int value)
{
    m_light->SetAttenuationValues(0,0,float(value)/100);

}


Genesyndata::~Genesyndata(){}

