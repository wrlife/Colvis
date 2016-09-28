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

class Genesyndata
{
public:
    Genesyndata();
    ~Genesyndata();
    vtkSmartPointer<vtkRenderer> getrenderer() {return m_renderer;}


    void rendermodel(vtkSmartPointer<vtkPolyData> t_model);

    void loadcamerapath(vtkSmartPointer<vtkPolyData> t_model);


    void updatecamera();


private:
    vtkSmartPointer<vtkPolyDataMapper> m_Mapper;
    vtkSmartPointer<vtkActor> m_Actor;
    vtkSmartPointer<vtkRenderer> m_renderer;
    vtkSmartPointer<vtkPolyData> m_camerapath;

    int counter;


};

#endif // GENESYNDATA_H
