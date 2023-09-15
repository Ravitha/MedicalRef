/*
  Copyright (c) 2019 AMIGO RESEARCH GROUP, <lalvarez@ulpgc.es>

  This program is free software: you can redistribute it and/or modify it under
  the terms of the GNU Affero General Public License as published by the Free
  Software Foundation, either version 3 of the License, or (at your option) any
  later version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
  details.

  You should have received a copy of the GNU Affero General Public License along
  with this program. If not, see <http://www.gnu.org/licenses/>.
*/

/**
 * \file main3DCurveSmoothing.cpp
 * \brief Main program to run the method for 3D curve smoothing
 * \author Luis Alvarez and Daniel Santana-Cedr�s
*/

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cstring>
#include <string.h>

#include "../auxiliary/image3D.h"
#include "curve3D.h"
#include "../auxiliary/obj3D.h"

using namespace std;

/** \brief Function to manage IPOL error demo file
 *
 * \param[in] mes: char vector with the message to write in the file
 * \return The error demo file with the message error or show an error if the
 *         file cannot be opened
 *
 */
void fprintf_demo_failure(char mes[300])
{
  FILE * f;
  f = fopen ("demo_failure.txt", "w");
  if(f==NULL) {
    printf("We cannot open the file demo_failure.txt");
    return;
  }
  fprintf(f,"%s\n",mes);
  fclose(f);
  return;
}

//==============================================================================

/** \brief Main program to compute a smoothed 3D curve
 *
 * \param[in] argc: number of provided parameters
 * \param[in] argv: vector with the following arguments:
 *    argv[1]: Name of ASCII file with the 3D curve point coordinates.
 *    argv[2]: Weight parameter w to balance the energy.
 *    argv[3]: Tolerance to stop iterations.
 *    argv[4]: Max number of iterations.
 *    argv[5]: Filename for the output smoothed 3D curve.
 *    argv[6]: output OBJ file to compare the original and smoothed curves.
 *    argv[7]: (OPTIONAL) 3D image (unsigned char in Analyze format) with a
 *             segmentation of a 3D object for which the given 3D curve is the
 *             centerline. If the image value is equal to zero the point is
 *             outside the 3D object.
 *    argv[8]: (Mandatory if you) 3D image (unsigned char in Analyze format) with a
 *             segmentation of a 3D object for which the given 3D curve is the
 *             centerline. If the image value is equal to zero the point is
 *             outside the 3D object.
 * \return Returns 0 in case of succeeding, 1 if the parameters checking fails,
 *         2 in case the curve file cannot be read, and 3 if the segmentation
 *         reading fails
 *
 */
int main(int argc, char *argv[])
{
  // WE CHECK THE NUMBER OF PARAMETERS
  if(argc<7) {
    char mes[300];
    sprintf(mes,"ERROR: argc=%d and it should be >=7\n",argc);
    printf("%s\n",mes);

    //Call description and parameters
    cout << "exe_file input_curve.txt w tol max_iter output_smoothed_curve.txt"
         << " 3Dcomparison.obj segmentation.hdr segmentation.img" << endl << endl;
    cout << "PARAMETERS:" << endl;
    cout << "\t* exe_file: exe file (called ./3DCurveSmoothing)." << endl;
    cout << "\t* input_curve.txt: name of ASCII file with the 3D curve point"
         << " coordinates." << endl;
    cout << "\t* w: weight parameter to balance the energy." << endl;
    cout << "\t* tol: tolerance to stop iterations." << endl;
    cout << "\t* max_iter: maximum number of iterations." << endl;
    cout << "\t* output_smoothed_curve.txt: filename for the output smoothed"
         << " 3D curve." << endl;
    cout << "\t* 3Dcomparison.obj: obj file to compare the original and"
         << " smoothed curves." << endl;
    cout << "\t* segmentation.hdr: optional 3D image (unsigned char in Analyze"
         << " format) with a segmentation of a 3D object for which the given 3D"
         << " curve is the centerline. If the image value is equal to zero the"
         << " point is outside the 3D object." << endl;
    cout << "\t* segmentation.img: Required if you use segmentation.hdr." << endl << endl;

    //Call examples (without/with the last parameter)
    cout << "EXAMPLES:" << endl;
    cout << "\tExample1: ./3DCurveSmoothing example/input_curve.txt"
         << " 1 0.00001 1000 example/output_curve_smoothed.txt"
         << " example/3DCurvesComparison_w1.obj" << endl;
    cout << "\tExample2: ./3DCurveSmoothing example/input_curve.txt"
         << " 1 0.00001 1000 example/output_curve_smoothed.txt"
         << " example/3DCurvesComparison_w1.obj"
         << " example/Segmentation.hdr"
         << " example/Segmentation.img" << endl;
    //Print error file
    fprintf_demo_failure(mes);
    return 1;
  }

  // WE READ THE FUNCTION PARAMETERS
  float w = atof(argv[2]);
  if(argc<8) {
    // Empiric adjustment of w parameter (in case no segmentation is provided)
    w *= 20;
  }
  float TOL   = atof(argv[3]);
  int MaxIter = atoi(argv[4]);

  // WE READ THE 3D CURVE
  vector <point3d<double> > C=read_curve3D(argv[1]);
  if(C.size()==0) {
    char mes[300];
    sprintf(mes,"Problem reading file %s\n",argv[1]);
    printf("%s\n",mes);
    fprintf_demo_failure(mes);
    return 2;
  }

  // WE READ OR BUILD THE 3D VOLUME AND COMPUTE THE DISTANCE FUNCTION
  image3D<unsigned char> I; // 3D volume
  image3D<float> u; // 3D image for the result of the distance function
  point3d<double> t(0.,0.,0.); // Translation vector to apply for coordinates
  // normalization (in case no volume is provided)
  double zoom = 1.; // Zoom factor to apply for coordinates normalization
                    // (in case no volume is provided)
  point3d<double> pmin = C[0];
  point3d<double> pmax = C[0];

  // In case no segmentation volume is provided
  if(argc<8) {
    // COMPUTATION OF MIN AND MAX OF POINT COORDINATES VALUES
    for(unsigned int k=1; k<C.size(); k++) {
      if(C[k].x>pmax.x)
        pmax.x = C[k].x;
      if(C[k].y>pmax.y)
        pmax.y = C[k].y;
      if(C[k].z>pmax.z)
        pmax.z = C[k].z;
      if(C[k].x<pmin.x)
        pmin.x = C[k].x;
      if(C[k].y<pmin.y)
        pmin.y = C[k].y;
      if(C[k].z<pmin.z)
        pmin.z = C[k].z;
    }
    printf("pmin=(%lf,%lf,%lf) pmax=(%lf,%lf,%lf)\n",
           pmin.x,pmin.y,pmin.z,pmax.x,pmax.y,pmax.z);
    double Size    = (0.001+pmax.x-pmin.x)*(0.001+pmax.y-pmin.y)*
                     (0.0001+pmax.z-pmin.z);
    double SizeMin = 1000000.;
    double SizeMax = 10000000.;

    // WE NORMALIZE THE COORDINATES OF THE 3D CURVE
    if(Size<SizeMin) {
      zoom = pow(SizeMin/Size,1./3.);
    } else if(Size>SizeMax) {
      zoom = pow(SizeMax/Size,1./3.);
    }
    printf("Size=%lf zoom=%lf\n",Size,zoom);
    if(zoom!=1) {
      pmin = pmin*zoom;
      pmax = pmax*zoom;
    }

    t = point3d<double>(10.,10.,10.)-pmin;
    for(unsigned int k=0; k<C.size(); k++) {
      C[k] = C[k]*zoom+t;
    }
  } else {
    
    if(argc<9){
      char mes[300];
      sprintf(mes,"Expected IMG gile %s\n",argv[7]);
      printf("%s\n",mes);
      fprintf_demo_failure(mes);
      return 3;
    }
    
    // WE READ THE SEGMENTATION VOLUME
    I = image3D<unsigned char>(argv[7], argv[8]);
    if(I.size()==0) {
      char mes[300];
      sprintf(mes,"Problem reading files %s or %s\n",argv[7], argv[8]);
      printf("%s\n",mes);
      fprintf_demo_failure(mes);
      return 3;
    }
    // WE NORMALIZE SEGMENTATION VALUES
    for(int k=0; k<I.size(); k++) {
      // We check if the point is inside (>0) the segmentation
      // (output points == 0)
      if(I[k]>0)
        I[k] = 1;
    }
  }

  // 3D CURVE REPARAMETERIZATION
  C = curve3D_reparameterization(C,1.);

  // WE COMPUTE THE 3D VOLUME IN CASE NO SEGMENTATION IS PROVIDED
  if(argc<8) {
    I = image3D<unsigned char>(20+pmax.x-pmin.x,20+pmax.y-pmin.y,
                               20+pmax.z-pmin.z,1.,1.,1.,0);
    for(unsigned int k=0; k<C.size(); k++) {
      I((int) round(C[k].x),(int) round(C[k].y),(int) round(C[k].z))=1;
    }
  }

  // WE COMPUTE THE DISTANCE FUNCTION
    u = distance_function3D(I,30.,30.,1);

  // WE SMOOTH THE curve3D
  vector <point3d<double> >  C2 = C;
  curve3D_smoothing(C2, u, w,  MaxIter, TOL);

  // WE UNDO THE POINT COORDINATE NORMALIZATION
  if(argc<8) {
    for(unsigned int k=0; k<C.size(); k++) {
      C[k] = (C[k]-t)/zoom;
    }
    for(unsigned int k=0; k<C2.size(); k++) {
      C2[k] = (C2[k]-t)/zoom;
    }
  }

  // WE WRITE THE OUTPUT SMOOTHED CURVE
  write_curve3D(C2,argv[5]);

  // WE BUILD AND SAVE AN OBJ FILE WITH THE RESULTS
  obj3D curve3Ds;

  vector<color> vColor(2);
  vColor[0] = color(255,0,0);
  vColor[1] = color(0,0,255);

  if(argc<8) {
    curve3Ds.addLineMM(C,vColor[0],0.2);
    curve3Ds.addLineMM(C2,vColor[1],0.2);
  } else {
    curve3Ds = *(new obj3D(I));

    vector<point3d<double> > aux(C.size());
    for(unsigned int j=0; j<aux.size(); j++)
      aux[j] = C[j]/I.voxel_size();
    curve3Ds.addLine(aux,vColor[0],0.2);

    aux.resize(C2.size());
    for(unsigned int j=0; j<aux.size(); j++)
      aux[j] = C2[j]/I.voxel_size();
    curve3Ds.addLine(aux,vColor[1],0.2);
    curve3Ds.addSegmentation(I,color(128,128,128),0.2);
  }


  curve3Ds.saveOBJ(argv[6]);

  return 0;

}
