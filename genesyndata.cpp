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
#include <vtkCellLocator.h>
#include <vtkPolygonalSurfacePointPlacer.h>


#include <vtkFeatureEdges.h>
#include <vtkPolyDataMapper.h>


#include <vtkImageGaussianSmooth.h>
#include <vtkImageCast.h>
#include <vtkImageGradient.h>
#include <vtkImageMagnitude.h>
#include <vtkImageNonMaximumSuppression.h>
#include <vtkImageConstantPad.h>
#include <vtkImageToStructuredPoints.h>
#include <vtkLinkEdgels.h>
#include <vtkThreshold.h>
#include <vtkGeometryFilter.h>
#include <vtkSubPixelPositionEdgels.h>
#include <vtkStripper.h>
#include <vtkPolyDataMapper.h>
#include <vtkImageLuminance.h>
#include <vtkProperty.h>
#include <math.h>




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
    m_light->SetAmbientColor(1,1,1);
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

     //Project and draw onto 3D surface
     assignColorAttribute(t_model);

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

        int * t_sz=t_renderwin->GetSize();
        t_width=t_sz[0];
        t_height=t_sz[1];

        vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter =
          vtkSmartPointer<vtkWindowToImageFilter>::New();
        windowToImageFilter->SetInput(t_renderwin);
        windowToImageFilter->SetMagnification(1); //set the resolution of the output image (3 times the current resolution of vtk render windowtotalcount)
        windowToImageFilter->SetInputBufferTypeToRGB(); //also record the alpha (transparency) channel
        windowToImageFilter->ReadFrontBufferOff(); // read from the back buffer
        windowToImageFilter->Update();


        vtkSmartPointer<vtkPNGWriter> writer =
          vtkSmartPointer<vtkPNGWriter>::New();

        std::string filename = "contour.png";
        writer->SetFileName(filename.c_str());
        writer->SetInputConnection(windowToImageFilter->GetOutputPort());
        writer->Write();


        //Opencv edge detection
        int goodedge = edgedetection(filename);


        if(goodedge==1){
            //Color convex points
            double contourcolor[3],surfacecolor[3];
            contourcolor[0]=255;contourcolor[1]=255;contourcolor[2]=255;
            surfacecolor[0]=0;surfacecolor[1]=0;surfacecolor[2]=255;
            coloronsurface(t_model,contourcolor,surfacecolor);

            //Color concave points
            contourcolor[0]=0;contourcolor[1]=255;contourcolor[2]=0;
            surfacecolor[0]=0;surfacecolor[1]=255;surfacecolor[2]=0;
            coloronsurface(t_model,contourcolor,surfacecolor);

            //Color parabolic points

            contourcolor[0]=0;contourcolor[1]=0;contourcolor[2]=255;
            surfacecolor[0]=255;surfacecolor[1]=0;surfacecolor[2]=0;
            coloronsurface(t_model,contourcolor,surfacecolor);
        }

}



//===============
//Initialization process, give the surface an color attribute so we can color it
//===============
void Genesyndata::assignColorAttribute(vtkSmartPointer<vtkPolyData> outputPolyData){
    vtkSmartPointer<vtkUnsignedCharArray> colors =
      vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors->SetNumberOfComponents(3);
    colors->SetName("Colors");
    for(int i = 0; i < outputPolyData->GetNumberOfPoints(); i++)
      {
      double p[3];
      outputPolyData->GetPoint(i,p);

      unsigned char color[3];
      color[0]=255;
      color[1]=255;
      color[2]=255;

      colors->InsertNextTupleValue(color);
      }

      outputPolyData->GetPointData()->SetScalars(colors);
}



//===============
// Assign different color to the surface point according to different shape type
//===============
void Genesyndata::coloronsurface(vtkSmartPointer<vtkPolyData> outputPolyData,double contourcolor[3],double surfacecolor[3])
{


    double color[3];

    //===================
    //Find closest point in 3D
    //====================
    // Create the tree
    vtkSmartPointer<vtkCellLocator> cellLocator =
      vtkSmartPointer<vtkCellLocator>::New();
    cellLocator->SetDataSet(outputPolyData);
    cellLocator->BuildLocator();
    vtkSmartPointer<vtkPolygonalSurfacePointPlacer> mplacer = vtkSmartPointer<vtkPolygonalSurfacePointPlacer>::New();
    mplacer->AddProp(m_Actor);


    for(int i = 0; i < m_contour.size(); i++){


        Vec3b tmpcolor = imageSource.at<Vec3b>(Point(m_contour[i].x, m_contour[i].y));
        if(tmpcolor[0]==contourcolor[0]&&tmpcolor[1]==contourcolor[1]&&tmpcolor[2]==contourcolor[2]){
            color[0]=surfacecolor[0];
            color[1]=surfacecolor[1];
            color[2]=surfacecolor[2];


            double pixel[2];

            int * t_sz = m_renderer->GetSize();
            t_width=t_sz[0];
            t_height=t_sz[1];
            pixel[0]=m_contour[i].x;
            pixel[1]=t_height-m_contour[i].y;

            double world[3];
            double worldOrient[9];

            int gotpoint = mplacer->ComputeWorldPosition(m_renderer,pixel,world,worldOrient);

            if(gotpoint==0)
            {
                continue;
            }

            //std::cout << "World coordinate: " << world[0] << ", " << world[1] << ", " << world[2] << std::endl;


            double closestPoint[3];//the coordinates of the closest point will be returned here
            double closestPointDist2; //the squared distance to the closest point will be returned here
            vtkIdType cellId; //the cell id of the cell containing the closest point will be returned here
            int subId; //this is rarely used (in triangle strips only, I believe)
            cellLocator->FindClosestPoint(world, closestPoint, cellId, subId, closestPointDist2);
            vtkSmartPointer<vtkIdList> cellPointIds =
              vtkSmartPointer<vtkIdList>::New();
            outputPolyData->GetCellPoints(cellId, cellPointIds);


            //view normal

            double* cp;
            double NV[3],NV_dot_NS[3];
            vtkSmartPointer<vtkCamera> t_camera=m_renderer->GetActiveCamera();
            cp = t_camera->GetPosition();
            NV[0]=cp[0]-world[0];
            NV[1]=cp[1]-world[1];
            NV[2]=cp[2]-world[2];

            //Distance to each point
            double NSS[3][3];
            double smallestdot;
            int perpcellid;
            double dist[3];

            vtkFloatArray* normalDataFloat =
              vtkFloatArray::SafeDownCast(outputPolyData->GetPointData()->GetArray("Normals"));

            for(vtkIdType j = 0; j < cellPointIds->GetNumberOfIds(); j++)
            {
                double p[3];
                outputPolyData->GetPoint(cellPointIds->GetId(j),p);

                dist[j] = sqrt((p[0]-world[0])*(p[0]-world[0])+(p[1]-world[1])*(p[1]-world[1])+(p[2]-world[2])*(p[2]-world[2]));

                normalDataFloat->GetTuple(cellPointIds->GetId(j), NSS[j]);
                NV_dot_NS[j] = NV[0]*NSS[j][0]+NV[1]*NSS[j][1]+NV[2]*NSS[j][2];
                if(j==0)
                {
                    smallestdot=dist[j];//NV_dot_NS[j];
                    perpcellid=j;
                }
                else if(NV_dot_NS[j]<smallestdot)
                {
                    smallestdot = dist[j];//NV_dot_NS[j];
                    perpcellid=j;
                }

            }
             outputPolyData->GetPointData()->GetScalars()->SetTuple3(cellPointIds->GetId(perpcellid),color[0],color[1],color[2]);
        }


    }

    outputPolyData->Modified();
    outputPolyData->GetPointData()->Modified();
    outputPolyData->GetPointData()->GetScalars()->Modified();

}



//===============
// Detecting occlusion contour and compute its signed curvature
//===============
int Genesyndata::edgedetection(std::string filename)
{
    Mat src; Mat src_gray; Mat binary;
    RNG rng(12345);
    /// Load source image and convert it to gray
    src = cv::imread( filename.c_str(), 1 );

    /// Convert image to gray and blur it
    cvtColor( src, src_gray, CV_BGR2GRAY );

    threshold( src_gray, binary, 5, 255,THRESH_BINARY|THRESH_BINARY );


    int erosion_type;
    int erosion_size = 3;
    erosion_type = MORPH_ELLIPSE;
    Mat element = getStructuringElement( erosion_type,
                                         Size( 2*erosion_size + 1, 2*erosion_size+1 ),
                                         Point( erosion_size, erosion_size ) );

    /// Apply the erosion operation
    dilate( binary, binary, element );
    erode(binary,binary,element);

    erosion_size=1;
    element = getStructuringElement( erosion_type,
                                         Size( 2*erosion_size + 1, 2*erosion_size+1 ),
                                         Point( erosion_size, erosion_size ) );
    erode(binary,binary,element);

    imwrite( "binary.png", binary );

    blur( src_gray, src_gray, Size(3,3) );

    Mat canny_output;
    vector<vector<Point> > contours;
    vector<Vec4i> hierarchy;

    int thresh = 100;

    /// Detect edges using canny
    Canny( binary, canny_output, thresh, thresh*2, 3 );
    /// Find contours
    findContours( canny_output, contours, hierarchy, CV_RETR_CCOMP, CV_CHAIN_APPROX_NONE, Point(0, 0) );


    int largestcontour=0;
    int largestcontourid=0;
    for(int i=0;i<contours.size();i++)
    {
        if (contours[i].size()>largestcontour)
        {
            largestcontour=contours[i].size();
            largestcontourid = i;
        }

    }

    if(hierarchy[largestcontourid][2] == -1)
    {
        return -1;
    }


    //Find out corner points in the contour
    vector< double > vecCurvature;
    vector<int> cornerindex;
    vecCurvature=computeCurvature(contours[largestcontourid],30,&cornerindex,30);


    //Draw those points on a image


    imageSource = Mat(t_height,t_width, CV_8UC3, cvScalar(0,0,0));

    for(int i=0;i<contours[largestcontourid].size();i++)
    {
        if (vecCurvature[i]<0)
        {
        //cout<<"x: "<<contours[largestcontourid][i].x<<"y: "<<contours[largestcontourid][i].y<<endl;

            imageSource.at<Vec3b>(Point(contours[largestcontourid][i].x, contours[largestcontourid][i].y))[0]=255;
            imageSource.at<Vec3b>(Point(contours[largestcontourid][i].x, contours[largestcontourid][i].y))[1]=255;
            imageSource.at<Vec3b>(Point(contours[largestcontourid][i].x, contours[largestcontourid][i].y))[2]=255;

        }
        if (vecCurvature[i]>0&&vecCurvature[i]!=std::numeric_limits<double>::infinity())
        {

            imageSource.at<Vec3b>(Point(contours[largestcontourid][i].x, contours[largestcontourid][i].y))[1]=255;

        }

    }

    for(int i = 0;i<cornerindex.size();i++)
    {

        imageSource.at<Vec3b>(Point(contours[largestcontourid][cornerindex[i]].x,contours[largestcontourid][cornerindex[i]].y))[0]=255;

    }


    postprocessCurvature(&imageSource,contours[largestcontourid]);

    m_contour = contours[largestcontourid];

    // Show in a window
    namedWindow( "Contours", CV_WINDOW_AUTOSIZE );
    imshow( "Contours", imageSource );

    return 1;


}



//===============
//Post processing of the signed curvature.
//===============
void Genesyndata::postprocessCurvature(Mat *imageSource,vector<Point> contours)
{
    //Filter out noisy curvature
    int kernalsize = 10;

    Vec3b color;

    for(int i=0; i<contours.size();i++)
    {
        int concavecount=0;
        int convexcount=0;
        color = imageSource->at<Vec3b>(Point(contours[i].x, contours[i].y));
        if(color[0] ==255&&color[1]==0&&color[2]==0)
        {
            continue;
        }
        for(int j=1;j<kernalsize;j++)
        {

            int iminus = i-j;
            iminus = iminus < 0 ? iminus + contours.size() : iminus;
            int iplus = i+j;
            iplus = iplus >= contours.size() ? iplus - contours.size() : iplus;

            color = imageSource->at<Vec3b>(Point(contours[iminus].x, contours[iminus].y));
            if(color[0]==255&&color[1]==255&&color[2]==255)
            {
                convexcount+=1;
            }
            else if (color[0]==0&&color[1]==255&&color[2]==0)
            {
                concavecount+=1;
            }
            color = imageSource->at<Vec3b>(Point(contours[iplus].x, contours[iplus].y));
            if(color[0]==255&&color[1]==255&&color[2]==255)
            {
                convexcount+=1;
            }
            else if (color[0]==0&&color[1]==255&&color[2]==0)
            {
                concavecount+=1;
            }
        }
        if(convexcount>1.5*concavecount)
        {
            imageSource->at<Vec3b>(Point(contours[i].x, contours[i].y))[0]=255;
            imageSource->at<Vec3b>(Point(contours[i].x, contours[i].y))[1]=255;
            imageSource->at<Vec3b>(Point(contours[i].x, contours[i].y))[2]=255;
        }
        else if(concavecount>1.5*convexcount)
        {
            imageSource->at<Vec3b>(Point(contours[i].x, contours[i].y))[0]=0;
            imageSource->at<Vec3b>(Point(contours[i].x, contours[i].y))[1]=255;
            imageSource->at<Vec3b>(Point(contours[i].x, contours[i].y))[2]=0;
        }

        //Find parabolic point
        if(abs(concavecount-convexcount)<3)
        {
            imageSource->at<Vec3b>(Point(contours[i].x, contours[i].y))[0]=0;
            imageSource->at<Vec3b>(Point(contours[i].x, contours[i].y))[1]=0;
            imageSource->at<Vec3b>(Point(contours[i].x, contours[i].y))[2]=255;
        }


    }
}




//===============
//Compute signed curvature of a given contour
//===============
vector< double > Genesyndata::computeCurvature(vector<Point> vecContourPoints,int step, vector<int> *cornerindex,int minAngle)
{
      std::vector< double > vecCurvature( vecContourPoints.size() );

      if (vecContourPoints.size() < step)
        return vecCurvature;

      auto frontToBack = vecContourPoints.front() - vecContourPoints.back();

      bool isClosed = ((int)std::max(std::abs(frontToBack.x), std::abs(frontToBack.y))) <= 1;


      if (isClosed)
      {
          cout<<"contourclose"<<endl;
      }



      //Find out corner points

      for (int i = 0; i < vecContourPoints.size(); i++ )
      {
          int iminus = i-3;
          int iplus = i+3;
          const cv::Point2f& pos = vecContourPoints[i];

          cv::Point2f pminus_corner = vecContourPoints[iminus < 0 ? iminus + vecContourPoints.size() : iminus];
          cv::Point2f pplus_corner = vecContourPoints[iplus > vecContourPoints.size() ? iplus - vecContourPoints.size() : iplus];

          //Calculate angle between current and previous point
          double adotb = (pos.y - pminus_corner.y)*(pplus_corner.y-pos.y)+(pos.x - pminus_corner.x)*(pplus_corner.x-pos.x);
          double anormdotbnorm = sqrt(pow(pplus_corner.x-pos.x,2)+pow(pplus_corner.y-pos.y,2))*sqrt(pow(pos.y - pminus_corner.y,2)+pow(pos.x - pminus_corner.x,2));

          double costheta = adotb/anormdotbnorm;
          if(costheta>1) costheta=1;

          double diffangle = abs(acos(costheta)*180.0/M_PI);

          if (diffangle>minAngle)
          {
            cornerindex->push_back(i);
          }
      }


      cv::Point2f f1stDerivative, f2ndDerivative;
      for (int i = 0; i < vecContourPoints.size(); i++ )
      {
          cv::Point2f pplus, pminus;
          cv::Point2f pos = vecContourPoints[i];

          int maxStep = step;

          //Central diff
          int iminus = i-maxStep;
          int iplus = i+maxStep;
          pminus = vecContourPoints[iminus < 0 ? iminus + vecContourPoints.size() : iminus];
          pplus = vecContourPoints[iplus > vecContourPoints.size() ? iplus - vecContourPoints.size() : iplus];

        int hascorner=0;
        //Check if corner points are within neighbor
        for(int u=0;u<cornerindex->size();u++)
        {
            if(iminus<(*cornerindex)[u]&&iplus>(*cornerindex)[u])
            {
              if((*cornerindex)[u]>i)
              {
                  hascorner=1; //Use backward difference
                  pplus=vecContourPoints[(*cornerindex)[u]];
                  pos = vecContourPoints[round((*cornerindex)[u]+iminus)/2];
                  iplus=(*cornerindex)[u];



              }
              else if ((*cornerindex)[u]<i)
              {
                  hascorner=2; //Use forwad difference
                  pminus=vecContourPoints[(*cornerindex)[u]];
                  pos = vecContourPoints[round((*cornerindex)[u]+iplus)/2];
                  iminus=(*cornerindex)[u];
              }
              else
              {
                  hascorner = 3;
                  pminus=pos;
                  pplus=pos;
              }
              break;
            }

        }


        f1stDerivative.x =   (pplus.x -        pminus.x) / (iplus-iminus);
        f1stDerivative.y =   (pplus.y -        pminus.y) / (iplus-iminus);
        f2ndDerivative.x = (pplus.x - 2*pos.x + pminus.x) / ((iplus-iminus)/2*(iplus-iminus)/2);
        f2ndDerivative.y = (pplus.y - 2*pos.y + pminus.y) / ((iplus-iminus)/2*(iplus-iminus)/2);

        double curvature2D;
        double divisor = f1stDerivative.x*f1stDerivative.x + f1stDerivative.y*f1stDerivative.y;
        if ( std::abs(divisor) > 10e-8 )
        {
          curvature2D =  (f2ndDerivative.y*f1stDerivative.x - f2ndDerivative.x*f1stDerivative.y) /
                pow(divisor, 3.0/2.0 )  ;
        }
        else
        {
          curvature2D = std::numeric_limits<double>::infinity();
        }

        vecCurvature[i] = curvature2D;

      }

      return vecCurvature;

}



void Genesyndata::addlight()
{
    //m_renderer->AddLight(m_light);
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

