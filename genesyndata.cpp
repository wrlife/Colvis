#include "genesyndata.h"
#include <iostream>
#include <vtkCenterOfMass.h>
#include <vtkProperty.h>
#include <vtkPolyDataNormals.h>
#include <vtkPointData.h>
#include <vtkLightCollection.h>

#include <vtkBMPWriter.h>
#include <vtkWindowToImageFilter.h>
#include <vtkImageShiftScale.h>
#include <vtkPNGWriter.h>
#include <fstream>
#include <vtkMatrix4x4.h>
#include <vtkFloatArray.h>
#include <string>
#include <sstream>
#include <vtkCoordinate.h>


#include <vtkImageCanvasSource2D.h>
#include <vtkJPEGWriter.h>







Genesyndata::Genesyndata()
{

//    m_sphereSource =
//      vtkSmartPointer<vtkSphereSource>::New();
    //Mapper
    m_Mapper =
        vtkSmartPointer<vtkPolyDataMapper>::New();
//      m_Mapper->SetInputConnection(m_sphereSource->GetOutputPort());

    //Actor
    m_Actor =
        vtkSmartPointer<vtkActor>::New();
    m_Actor->SetMapper(m_Mapper);
    m_Actor->GetProperty()->SetInterpolationToGouraud();

    //light
    m_light = vtkSmartPointer<vtkLight>::New();
    m_light->SetLightTypeToHeadlight();
    m_light->SetAttenuationValues(0,0.2,0);
    m_light->SetPositional(true); // required for vtkLightActor below
    m_light->SetConeAngle(180);
    m_light->SetDiffuseColor(50,50,50);
    m_light->SetAmbientColor(10,10,10);
    m_light->SetSpecularColor(50,50,50);
    m_light->SetLightTypeToHeadlight();

    //Renderer
    m_renderer =
        vtkSmartPointer<vtkRenderer>::New();

    m_renderer->AddActor(m_Actor);

    counter=0;
    totalcount=0;

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

     t_model->DeepCopy( normalGenerator->GetOutput());

     m_Mapper->SetInputData(t_model);

     m_renderer->ResetCamera();

}

//==================
//Render 3D parametric model
//==================
void Genesyndata::renderparametricmodel()
{

    parametricObjects = vtkSmartPointer<vtkParametricBoy>::New();
    parametricFunctionSources= vtkSmartPointer<vtkParametricFunctionSource>::New();

    parametricFunctionSources->SetParametricFunction(parametricObjects);
    parametricFunctionSources->SetUResolution(51);
    parametricFunctionSources->SetVResolution(51);
    parametricFunctionSources->SetWResolution(51);
    parametricFunctionSources->Update();

    m_Mapper->SetInputConnection(parametricFunctionSources->GetOutputPort());

    m_renderer->ResetCamera();

}



//==================
//Store camera path
//==================
void Genesyndata::loadcamerapath(vtkSmartPointer<vtkPolyData> t_model)
{

    //m_renderer->ResetCamera();
    m_camerapath=t_model;
    m_numcams=t_model->GetNumberOfPoints();

    my_camera->setcam(m_camerapath);

}




//==================
//Set current camera position and look at
//==================
void Genesyndata::updatecamera(int c_step)
{
    //Get inital camera position

    double* p;
    //m_camerapath->GetPoint(counter,p);

    my_camera->updatepos(c_step);

    p=my_camera->getcampos();

    vtkSmartPointer<vtkCamera> t_camera=m_renderer->GetActiveCamera();
    t_camera->SetPosition(p[0],p[1],p[2]);

    p=my_camera->getcamdirection();

    t_camera->SetFocalPoint(p[0],p[1],p[2]);
    totalcount+=1;

}


void Genesyndata::randomcampos(float x,float y,float z, float elevation,float azimuth)
{
    double p[3];
    m_camerapath->GetPoint(counter-1,p);
    vtkSmartPointer<vtkCamera> t_camera=m_renderer->GetActiveCamera();


    t_camera->SetPosition(p[0]+x,p[1]+y,p[2]+z);

    m_camerapath->GetPoint(counter,p);

    t_camera->SetFocalPoint(p[0]+x,p[1]+y,p[2]+z);

    t_camera->Elevation(elevation);
    t_camera->Azimuth(azimuth);

    totalcount+=1;
}


//===============
// Obtain the depth map observed at current camera view
//===============

void Genesyndata::get_z_values(vtkRenderWindow*t_renderwin)
{


    vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter =
      vtkSmartPointer<vtkWindowToImageFilter>::New();
    windowToImageFilter->SetInput(t_renderwin);
    windowToImageFilter->SetMagnification(1); //set the resolution of the output image (3 times the current resolution of vtk render windowtotalcount)
    windowToImageFilter->SetInputBufferTypeToRGBA(); //also record the alpha (transparency) channel
    windowToImageFilter->ReadFrontBufferOff(); // read from the back buffer
    windowToImageFilter->Update();

    vtkSmartPointer<vtkPNGWriter> writer =
      vtkSmartPointer<vtkPNGWriter>::New();

    std::ostringstream oss;

    int t_count=totalcount;
    oss<<"frame"<<setfill('0')<<setw(5)<<t_count<<".png";

    std::string filename=oss.str();

    writer->SetFileName( ("./images/"+filename).c_str());
    writer->SetInputConnection(windowToImageFilter->GetOutputPort());
    writer->Write();

    vtkSmartPointer<vtkCamera> t_camera=m_renderer->GetActiveCamera();


    double* cliprange=t_camera->GetClippingRange();

    double nearclip=cliprange[0];
    double farclip=cliprange[1];


    int * t_sz=t_renderwin->GetSize();
    int t_width=t_sz[0];
    int t_height=t_sz[1];


    float *zvalues = new float[t_width * t_height];

    t_renderwin->GetZbufferData(0,0,t_width-1,t_height-1,zvalues);

    //rescale z from [0,1] to [-1,1]
    for (int i=0;i<t_width*t_height;i++){
        zvalues[i]=zvalues[i]*2-1;
        zvalues[i]=2*nearclip*farclip/(farclip+nearclip-zvalues[i]*(farclip-nearclip));
    }

    std::string depthfilename="./depth/"+filename+".bin";
    ofstream myFile(depthfilename.c_str(),ios::out | ios::binary);
    myFile.write((char*)zvalues,sizeof(float)*t_width*t_height);

    delete zvalues;

}


//===============
// Obtain points where viewing ray is perpendicular to surface normal
//===============
void Genesyndata::get_orthognal_normal_view(vtkSmartPointer<vtkPolyData> t_model,vtkRenderWindow*t_renderwin)
{
    double* cp;
    vtkSmartPointer<vtkCamera> t_camera=m_renderer->GetActiveCamera();
    cp = t_camera->GetPosition();

    vtkFloatArray* normalDataFloat =
      vtkFloatArray::SafeDownCast(t_model->GetPointData()->GetArray("Normals"));


    //Compute dot product


    std::vector<mpoint> contour;
    if(normalDataFloat)
    {

        for(vtkIdType i = 0; i < normalDataFloat->GetNumberOfTuples(); i++)
        {
          double p[3];
          double NV[3];
          double NS[3];
          normalDataFloat->GetTuple(i, NS);

          t_model->GetPoint(i,p);

          //view normal
          NV[0]=cp[0]-p[0];
          NV[1]=cp[1]-p[1];
          NV[2]=cp[2]-p[2];


          //dot product between NV and NS
          int NV_dot_NS = NV[0]*NS[0]+NV[1]*NS[1]+NV[2]*NS[2];


          if(std::abs(NV_dot_NS)<0.0001)
          {
              mpoint temppoint;
              temppoint.x=p[0];
              temppoint.y=p[1];
              temppoint.z=p[2];
              contour.push_back(temppoint);
          }

        }

        vtkSmartPointer<vtkCoordinate> coordinate =
          vtkSmartPointer<vtkCoordinate>::New();

        coordinate->SetCoordinateSystemToWorld();

        //Draw those points on a image



        int * t_sz=t_renderwin->GetSize();
        int t_width=t_sz[0];
        int t_height=t_sz[1];

        int extent[6] = {0,t_width,0,t_height,0,0};

        vtkSmartPointer<vtkImageCanvasSource2D> imageSource =
          vtkSmartPointer<vtkImageCanvasSource2D>::New();
        imageSource->SetExtent( extent );
        imageSource->SetScalarTypeToUnsignedChar(); // vtkJPEGWriter only accepts unsigned char input
        imageSource->SetNumberOfScalarComponents( 3 ); // 3 color channels: Red, Green and Blue

        imageSource->SetDrawColor(0.0, 0.0, 0.0);
        imageSource->FillBox(extent[0],extent[1],extent[2],extent[3]);


        imageSource->SetDrawColor( 0, 127, 255 );


        for(int i=0;i<contour.size();i++)
        {

          coordinate->SetValue(contour[i].x,contour[i].y,contour[i].z);
          int* val;
          val = coordinate->GetComputedDisplayValue(m_renderer);

          imageSource->DrawPoint(val[0],val[1]);


        }

        std::string outputFilename = "output.jpg";
        vtkSmartPointer<vtkJPEGWriter> writer =
          vtkSmartPointer<vtkJPEGWriter>::New();
        writer->SetFileName( outputFilename.c_str() );
        writer->SetInputConnection( imageSource->GetOutputPort() );
        writer->Write();


    }


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

void Genesyndata::setambient(int value)
{
    m_light->SetAmbientColor(value,value,value);

}

void Genesyndata::setdiffuse(int value)
{
    m_light->SetDiffuseColor(value,value,value);

}

void Genesyndata::setspecular(int value)
{
    m_light->SetSpecularColor(value,value,value);

}



Genesyndata::~Genesyndata(){}

