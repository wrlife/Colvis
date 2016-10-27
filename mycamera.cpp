#include "mycamera.h"

Mycamera::Mycamera(int steps)
{
    m_steps=steps;
    for (int i=0;i<3;i++){
        m_vel[i]=0;
        m_vel_next[i]=0;
    }
}



void Mycamera::updatepos(int c_step)
{
    int c_point=(int)floor((c_step+1) / m_steps);

    if((c_step+1) % m_steps==0||c_step==0)
    {
        double p[3],p_n[3];
        m_camerapath->GetPoint(c_point,p);
        m_camerapath->GetPoint(c_point+1,p_n);

        m_camerapath->GetPoint(c_point,m_position);
        m_camerapath->GetPoint(c_point+1,m_direction);

        for(int i=0;i<3;i++)
        {
            m_vel[i]=(p_n[i]-p[i])/m_steps;

        }



        m_camerapath->GetPoint(c_point+2,p);


        for(int i=0;i<3;i++)
        {
            m_vel_next[i]=-(p_n[i]-p[i])/m_steps;

        }


    }

    else{
        for (int i=0;i<3;i++){
           m_position[i]+=m_vel[i];
           m_direction[i]+=m_vel_next[i];

        }

    }

}
