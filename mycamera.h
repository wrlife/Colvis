#ifndef MYCAMERA_H
#define MYCAMERA_H

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>


class Mycamera
{
public:
    Mycamera(int steps);

    void setcam(vtkSmartPointer<vtkPolyData> campath){m_camerapath=campath;}
    void updatepos(int c_step);

    double* getcampos(){return m_position;}
    double* getcamdirection(){return m_direction;}

private:
    vtkSmartPointer<vtkPolyData> m_camerapath;
    double m_vel[3];
    double m_vel_next[3];
    int m_steps;
    double m_position[3];
    double m_direction[3];
};

#endif // MYCAMERA_H
