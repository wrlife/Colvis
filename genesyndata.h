#ifndef GENESYNDATA_H
#define GENESYNDATA_H

#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkCamera.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkLight.h>
#include <vtkLightActor.h>
#include <vtkWindowToImageFilter.h>
#include "mycamera.h"
#include <vtkParametricBoy.h>
#include <vtkParametricFunctionSource.h>

#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>



using namespace cv;
using namespace std;

typedef struct {
  double x, y, z;
} mpoint;

class Genesyndata
{
public:
    Genesyndata();
    ~Genesyndata();
    vtkSmartPointer<vtkRenderer> getrenderer() {return m_renderer;}
    void addlight();


    void rendermodel(vtkSmartPointer<vtkPolyData> t_model);

    void loadcamerapath(vtkSmartPointer<vtkPolyData> t_model);

    void get_z_values(vtkRenderWindow*t_renderwin);

    int get_num_cams(){return m_numcams;}

    vtkSmartPointer<vtkPolyData> getparametricdata() {return parametricFunctionSources->GetOutput();}


    void updatecamera(int c_step);
    void setconstantlight(int value);
    void setlinearlight(int value);
    void setqudraticlight(int value);
    void setambient(int value);
    void setdiffuse(int value);
    void setspecular(int value);
    void randomcampos(float x,float y,float z, float elevation,float azimuth);

    void get_orthognal_normal_view(vtkSmartPointer<vtkPolyData> t_model, vtkRenderWindow *t_renderwin);
    void renderparametricmodel();

    void edgedetection(std::string filename);
    std::vector< double > computeCurvature(std::vector<cv::Point> vecContourPoints, int step, vector<int> *cornerindex, int minAngle);


private:
     vtkSmartPointer<vtkSphereSource> m_sphereSource;

    vtkSmartPointer<vtkLight> m_light;
    vtkSmartPointer<vtkPolyDataMapper> m_Mapper;
    vtkSmartPointer<vtkActor> m_Actor;
    vtkSmartPointer<vtkRenderer> m_renderer;
    vtkSmartPointer<vtkPolyData> m_camerapath;
    vtkSmartPointer<vtkLightActor> m_lightActor;

    vtkSmartPointer<vtkParametricBoy> parametricObjects;
    vtkSmartPointer<vtkParametricFunctionSource> parametricFunctionSources;

    Mycamera* my_camera=new Mycamera(20);

    int counter;
    int totalcount;
    int m_numcams;

    int t_width;
    int t_height;


};

#endif // GENESYNDATA_H
