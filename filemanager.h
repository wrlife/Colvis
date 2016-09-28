#ifndef FILELOADER_H
#define FILELOADER_H
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <QString>
#include <vtkSTLReader.h>

class Filemanager
{
public:
    Filemanager();
    ~Filemanager();

    void loadnewfile(QString filename);

    void loadnewcamera(QString filename);

    vtkSmartPointer<vtkPolyData> getfile() {return m_polydata;}


private:
     vtkSmartPointer<vtkXMLPolyDataReader> m_vtp;
     vtkSmartPointer<vtkSTLReader> m_stl;
     vtkSmartPointer<vtkXMLPolyDataReader> m_camera;
     vtkSmartPointer<vtkPolyData> m_polydata;

};

#endif // FILELOADER_H
