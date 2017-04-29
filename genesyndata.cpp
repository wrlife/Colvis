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



//    double* cp;
//    vtkSmartPointer<vtkCamera> t_camera=m_renderer->GetActiveCamera();
//    cp = t_camera->GetPosition();

//    vtkFloatArray* normalDataFloat =
//      vtkFloatArray::SafeDownCast(t_model->GetPointData()->GetArray("Normals"));


//    //Compute dot product

//    int * t_sz=t_renderwin->GetSize();
//    int t_width=t_sz[0];
//    int t_height=t_sz[1];

//    // Create the tree
//    vtkSmartPointer<vtkCellLocator> cellLocator =
//      vtkSmartPointer<vtkCellLocator>::New();
//    cellLocator->SetDataSet(t_model);
//    cellLocator->BuildLocator();


//    vtkSmartPointer<vtkPolygonalSurfacePointPlacer> mplacer = vtkSmartPointer<vtkPolygonalSurfacePointPlacer>::New();
//    mplacer->AddProp(m_Actor);

//    std::vector<mpoint> contour;

//    for(int i = 0;i<t_width;i++)
//    {
//        for(int j=0; j<t_height;j++)
//        {

//            double pixel[2];

//            pixel[0]=i;
//            pixel[1]=j;

//            double world[3];
//            double worldOrient[9];

//            mplacer->ComputeWorldPosition(m_renderer,pixel,world,worldOrient);

//           // std::cout << "World coordinate: " << world[0] << ", " << world[1] << ", " << world[2] << std::endl;


//            double closestPoint[3];//the coordinates of the closest point will be returned here
//            double closestPointDist2; //the squared distance to the closest point will be returned here
//            vtkIdType cellId; //the cell id of the cell containing the closest point will be returned here
//            int subId; //this is rarely used (in triangle strips only, I believe)
//            cellLocator->FindClosestPoint(world, closestPoint, cellId, subId, closestPointDist2);

////            std::cout << "Coordinates of closest point: " << closestPoint[0] << " " << closestPoint[1] << " " << closestPoint[2] << std::endl;
////            std::cout << "Squared distance to closest point: " << closestPointDist2 << std::endl;
////            std::cout << "CellId: " << cellId << std::endl;

//            vtkSmartPointer<vtkIdList> cellPointIds =
//              vtkSmartPointer<vtkIdList>::New();
//            t_model->GetCellPoints(cellId, cellPointIds);


//            //Distance to each point
//            double dist[3];
//            double NSS[3][3];

//            for(vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++)
//            {
//                double p[3];
//                t_model->GetPoint(cellPointIds->GetId(i),p);

//                dist[i] = sqrt((p[0]-world[0])*(p[0]-world[0])+(p[1]-world[1])*(p[1]-world[1])+(p[2]-world[2])*(p[2]-world[2]));

////                std::cout << "Cell point "<<i<<  ": "<< p[0] << " " << p[1] << " " << p[2] << std::endl;
////                std::cout << "Dist: " << dist[i] << std::endl;

//                normalDataFloat->GetTuple(cellPointIds->GetId(i), NSS[i]);

//            }

//            double sumdist = dist[0]+dist[1]+dist[2];

//            double NS[3];
//            double NV[3];

//            NS[0] = dist[0]/sumdist*NSS[0][0]+dist[1]/sumdist*NSS[1][0]+dist[2]/sumdist*NSS[2][0];
//            NS[1] = dist[0]/sumdist*NSS[0][1]+dist[1]/sumdist*NSS[1][1]+dist[2]/sumdist*NSS[2][1];
//            NS[2] = dist[0]/sumdist*NSS[0][2]+dist[1]/sumdist*NSS[1][2]+dist[2]/sumdist*NSS[2][2];



//            //view normal
//            NV[0]=cp[0]-world[0];
//            NV[1]=cp[1]-world[1];
//            NV[2]=cp[2]-world[2];

//            //dot product between NV and NS
//            int NV_dot_NS = NV[0]*NS[0]+NV[1]*NS[1]+NV[2]*NS[2];
//            if(std::abs(NV_dot_NS)<0.0001)
//            {
//                mpoint temppoint;
//                temppoint.x=world[0];
//                temppoint.y=world[1];
//                temppoint.z=world[2];
//                contour.push_back(temppoint);
//            }

//        }
//    }


//        vtkSmartPointer<vtkCoordinate> coordinate =
//          vtkSmartPointer<vtkCoordinate>::New();

//        coordinate->SetCoordinateSystemToWorld();

////        //Draw those points on a image
//        int extent[6] = {0,t_width,0,t_height,0,0};

//        vtkSmartPointer<vtkImageCanvasSource2D> imageSource =
//          vtkSmartPointer<vtkImageCanvasSource2D>::New();
//        imageSource->SetExtent( extent );
//        imageSource->SetScalarTypeToUnsignedChar(); // vtkJPEGWriter only accepts unsigned char input
//        imageSource->SetNumberOfScalarComponents( 3 ); // 3 color channels: Red, Green and Blue

//        imageSource->SetDrawColor(0.0, 0.0, 0.0);
//        imageSource->FillBox(extent[0],extent[1],extent[2],extent[3]);


//        imageSource->SetDrawColor( 0, 127, 255 );


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
        edgedetection(filename);

        //Find inflection of the contour point








}


void Genesyndata::edgedetection(std::string filename)
{
    Mat src; Mat src_gray; Mat binary;
    RNG rng(12345);
    /// Load source image and convert it to gray
    src = cv::imread( filename.c_str(), 1 );

    /// Convert image to gray and blur it
    cvtColor( src, src_gray, CV_BGR2GRAY );

    threshold( src_gray, binary, 10, 255,THRESH_BINARY|THRESH_BINARY );

    imshow( "Contours", binary );

    blur( src_gray, src_gray, Size(3,3) );

    Mat canny_output;
    vector<vector<Point> > contours;
    vector<Vec4i> hierarchy;

    int thresh = 100;

    /// Detect edges using canny
    Canny( binary, canny_output, thresh, thresh*2, 3 );
    /// Find contours
    findContours( canny_output, contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_NONE, Point(0, 0) );

    /// Draw contours
    Mat drawing = Mat::zeros( canny_output.size(), CV_8UC3 );


    Scalar color = Scalar( rng.uniform(0, 255), rng.uniform(0,255), rng.uniform(0,255) );

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

    drawContours( drawing, contours, largestcontourid, color, 1, 8, hierarchy, 0, Point() );


    /// Show in a window
    namedWindow( "Contours", CV_WINDOW_AUTOSIZE );
    imshow( "Contours", drawing );


    //Find out corner points in the contour



    vector< double > vecCurvature;
    vector<int> cornerindex;
    vecCurvature=computeCurvature(contours[largestcontourid],20,&cornerindex,40);


    //Draw those points on a image

    int extent[6] = {0,t_width,0,t_height,0,0};

    vtkSmartPointer<vtkImageCanvasSource2D> imageSource =
      vtkSmartPointer<vtkImageCanvasSource2D>::New();
    imageSource->SetExtent( extent );
    imageSource->SetScalarTypeToUnsignedChar(); // vtkJPEGWriter only accepts unsigned char input
    imageSource->SetNumberOfScalarComponents( 3 ); // 3 color channels: Red, Green and Blue

    imageSource->SetDrawColor(0.0, 0.0, 0.0);
    imageSource->FillBox(extent[0],extent[1],extent[2],extent[3]);



    for(int i=0;i<contours[largestcontourid].size();i++)
    {
        if (vecCurvature[i]<0)
        {
            imageSource->SetDrawColor( 0, 127, 255 );
            imageSource->DrawPoint(contours[largestcontourid][i].x,contours[largestcontourid][i].y);
        }
        if (vecCurvature[i]>0)
        {
            imageSource->SetDrawColor( 255, 255, 255 );
            imageSource->DrawPoint(contours[largestcontourid][i].x,contours[largestcontourid][i].y);

        }
        if (vecCurvature[i]==0)
        {
            imageSource->SetDrawColor( 255, 0, 0 );
            imageSource->DrawPoint(contours[largestcontourid][i].x,contours[largestcontourid][i].y);
        }
        if (vecCurvature[i]==std::numeric_limits<double>::infinity())
        {
            imageSource->SetDrawColor( 0, 255, 0 );
            imageSource->DrawPoint(contours[largestcontourid][i].x,contours[largestcontourid][i].y);
        }

    }

    for(int i = 0;i<cornerindex.size();i++)
    {
        //imageSource->SetDrawColor( 255, 0, 0 );
        //imageSource->DrawPoint(contours[largestcontourid][cornerindex[i]].x,contours[largestcontourid][cornerindex[i]].y);
    }

    vtkSmartPointer<vtkPNGWriter> writer =
      vtkSmartPointer<vtkPNGWriter>::New();

    std::string mfilename = "test.png";
    writer->SetFileName(mfilename.c_str());
    writer->SetInputConnection(imageSource->GetOutputPort());
    writer->Write();



}

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
          int iminus = i-5;
          int iplus = i+5;
          const cv::Point2f& pos = vecContourPoints[i];

          cv::Point2f pminus_corner = vecContourPoints[iminus < 0 ? iminus + vecContourPoints.size() : iminus];
          cv::Point2f pplus_corner = vecContourPoints[iplus > vecContourPoints.size() ? iplus - vecContourPoints.size() : iplus];

          //Calculate angle between current and previous point
          double previousAngle = abs(atan2( pos.y - pminus_corner.y, pos.x - pminus_corner.x ) * 180 / M_PI);
          double currentAngle = abs(atan2( pplus_corner.y-pos.y, pplus_corner.x-pos.x ) * 180 / M_PI);

          double diffangle = abs(previousAngle-currentAngle);


          if (diffangle>minAngle)
          {
            cornerindex->push_back(i);
          }
      }


      cv::Point2f f1stDerivative, f2ndDerivative;
      for (int i = 0; i < vecContourPoints.size(); i++ )
      {
          vector<cv::Point2f> pplus, pminus;
          const cv::Point2f& pos = vecContourPoints[i];

          int maxStep = step;
          if (!isClosed)
            {
              maxStep = std::min(std::min(step, i), (int)vecContourPoints.size()-1-i);
              if (maxStep == 0)
                {
                  vecCurvature[i] = std::numeric_limits<double>::infinity();
                  continue;
                }
            }

//          int iminus = i-maxStep;
//          int iplus = i+maxStep;
          int iminus;
          int iplus;
          for(int j=1;j<4;j++)
          {
              iminus = i-maxStep*j;
              iplus = i+maxStep*j;
              pminus.push_back( vecContourPoints[iminus < 0 ? iminus + vecContourPoints.size() : iminus]);
              pplus.push_back( vecContourPoints[iplus > vecContourPoints.size() ? iplus - vecContourPoints.size() : iplus]);



          }

          int hascorner=0;
          //Check if corner points are within neighbor
          for(int u=0;u<cornerindex->size();u++)
          {
              if(iminus<(*cornerindex)[u]&&iplus>(*cornerindex)[u])
              {
                if((*cornerindex)[u]>i)
                {
                    hascorner=1; //Use backward difference
                    for(int v=0;v<3;v++){
                        pplus[v]=pos;
                    }

                }
                else if ((*cornerindex)[u]<i)
                {
                    hascorner=2; //Use forwad difference

                    for(int v=0;v<3;v++){
                        pminus[v]=pos;
                    }
                }
                else
                {
                    hascorner = 3;
                    for(int v=0;v<3;v++){
                        pminus[v]=pos;
                        pplus[v]=pos;
                    }

                }
                break;
              }
          }



           //cout<<pplus[2].x<<" "<<pplus[1].x<< " "<<pplus[0].x<<endl;

//          pminus = vecContourPoints[iminus < 0 ? iminus + vecContourPoints.size() : iminus];
//          pplus = vecContourPoints[iplus > vecContourPoints.size() ? iplus - vecContourPoints.size() : iplus];


          f1stDerivative.x =   (pplus[2].x-9*pplus[1].x+45*pplus[0].x -  45*pminus[0].x+9*pminus[1].x-pminus[2].x) / (60*maxStep);//(iplus-iminus);
          f1stDerivative.y =   (pplus[2].y-9*pplus[1].y+45*pplus[0].y -  45*pminus[0].y+9*pminus[1].y-pminus[2].y) / (60*maxStep);//(pplus.y -        pminus.y) / (iplus-iminus);
          f2ndDerivative.x = (2*pplus[2].x-27*pplus[1].x+270*pplus[0].x-490*pos.x+270*pminus[0].x-27*pminus[1].x+2*pminus[2].x)/(180*maxStep*maxStep);//(pplus.x - 2*pos.x + pminus.x) / ((iplus-iminus)/2*(iplus-iminus)/2);
          f2ndDerivative.y = (2*pplus[2].y-27*pplus[1].y+270*pplus[0].y-490*pos.y+270*pminus[0].y-27*pminus[1].y+2*pminus[2].y)/(180*maxStep*maxStep);//(pplus.y - 2*pos.y + pminus.y) / ((iplus-iminus)/2*(iplus-iminus)/2);

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

