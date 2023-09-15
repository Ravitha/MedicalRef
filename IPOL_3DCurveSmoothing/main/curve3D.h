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
 * \file curve3D.h
 * \brief Methods to manage 3D curves
 * \author Luis Alvarez and Daniel Santana-Cedrés
 */


#ifndef CURVE3D_H
#define CURVE3D_H
#include <iostream>
#include <cstdlib>
#include "../auxiliary/image3D.h"

// curve3D
vector <point3d<double> > curve3D_reparameterization(vector <point3d<double> >
    &Cl,const double length_segment,const double max_length=1e100);
//------------------------------------------------------------------------------
void curve3D_smoothing(vector <point3d<double> > &Cl,const image3D<float> &u,
                       float w=1.,int MaxIter=100000,float TOL=1e-4,
                       vector <point3d<double> > *GT=NULL);
//------------------------------------------------------------------------------
double length(const vector <point3d<double> > &Cl, int k0=0, int k1=9999999);
//------------------------------------------------------------------------------
vector <point3d<double> > read_curve3D(char name[400]);
//------------------------------------------------------------------------------
void write_curve3D(vector <point3d<double> > &C, char name[400]);
//------------------------------------------------------------------------------
// AVERAGE DISTANCE BETWEEN 2 CURVES
double distance_curve2curve(vector< point3d<double> > &c1,
                            vector< point3d<double> > &c2);
//------------------------------------------------------------------------------
// COMPUTATION OF 3D SIGNED DISTANCE FUNCTION FROM AN IMAGE AS REGION INDICATOR
template <class  T>
image3D<float> distance_function3D(const image3D<T> &Is,
                                   const float &MaxDistancePositive,
                                   const float &MaxDistanceNegative,
                                   const int &neighborhood_type);

// POINT NEIGHBORHOOD STRUCTURES IN 3D
vector<int> Neighborhood(int width,int size2d,int type);
vector<double> Distace2Neighbors(const double &vx,const double &vy,
                                 const double &vz,const int &type);
void change_order_neighborhood(vector<int> &N,vector<double> &Nd);

//==============================================================================

/** \brief Computation of 3D signed distance function from an image as region
 *         indicator
 *
 * \param[in,out] Is: Region indicator Is[m]>0 means that the voxel 'm' is
 *                    inside the region, Is[m]<=0 outside the region.
 * \param[in] MaxDistancePositive: Maximum distance computed inside the region
 * \param[in] MaxDistanceNegative: Maximum distance computed outside the region
 * \param[in] neighborhood_type: Neighborhood type (==0 : 6 voxels,
 *                               ==1 : 18 voxels, ==2 : 26 voxels)
 * \return It returns the distance function to the region boundary (positive
 *         inside the region and negative outside it)
 *
 */
template <class  T>
image3D<float> distance_function3D(const image3D<T>&Is,
                                   const float &MaxDistancePositive,
                                   const float &MaxDistanceNegative,
                                   const int &neighborhood_type)
{
  printf("init distance_function3D()...");

  // AUXILIARY VARIABLES
  int width   = Is.width();
  int height  = Is.height();
  int nSlices = Is.nSlices();
  int size2d  = Is.size2d();

  printf("nSlices=%d\n",nSlices);

  // DEFINITION OF VOXEL NEIGBORHOOD AND DISTANCE
  vector<int> N     = Neighborhood(Is.width(),Is.size2d(),neighborhood_type);
  vector<double> Nd = Distace2Neighbors(Is.voxel_size_x(),Is.voxel_size_y(),
                                      Is.voxel_size_z(),neighborhood_type);
  if(Is.voxel_size_z()<Is.voxel_size_x())
    change_order_neighborhood(N,Nd);

  // INT IMAGE TO CONTROL VISITED VOXELS
  image3D<int> visited(Is.width(),Is.height(),Is.nSlices(),Is.voxel_size_x(),
                       Is.voxel_size_y(),Is.voxel_size_z(),1e8);

  // OUTPUT DISTANCE FUNCTION
  image3D<float> u(width,height,nSlices,Is.voxel_size_x(),Is.voxel_size_y(),
                   Is.voxel_size_z());
  for(int k=0; k<u.size(); k++)
    u[k]=Is[k]>0?MaxDistancePositive:-MaxDistanceNegative;

  // ESTIMATION OF THE INITIAL CONTOUR, DISTANCE FUNCTION AND visited ARRAY
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for(int z=1; z<nSlices-1; z++) {
    for(int x=1; x<width-1; x++) {
      for(int y=1; y<height-1; y++) {
        int m=x+y*width+z*size2d;
        if(Is[m]>0) {
          for(unsigned int k=0; k<N.size(); k++) {
            if(Is[m+N[k]]<=0) {
              u[m]       = 0.5*Nd[k];
              visited[m] = 1;
              break;
            }
          }
        } else {
          for(unsigned int k=0; k<N.size(); k++) {
            if(Is[m+N[k]]>0) {
              u[m]       = -0.5*Nd[k];
              visited[m] = -1;
              break;
            }
          }
        }
      }
    }
  }
  // COMPUTATION OF EVOLUTION FOR POSITIVE PART
  if(MaxDistancePositive>0.) {
    // COMPUTATION OF THE DISTANCE FUNCTION IN AN ITERATIVE WAY
    float maxDistance = 0.;
    int iter          = 0;
    while(maxDistance<MaxDistancePositive) {
      iter++;
      float maxDistance2 = maxDistance;
      printf("%1.2lf ",maxDistance);
      #ifdef _OPENMP
      #pragma omp parallel for
      #endif
      for(int z=1; z<nSlices-1; z++) {
        for(int x=1; x<width-1; x++) {
          for(int y=1; y<height-1; y++) {
            int m = x+y*width+z*size2d;
            if(Is[m]<=0 || visited[m]!=iter)
              continue;
            for(unsigned int k=0; k<N.size(); k++) {
              if(Is[m+N[k]]>0) {
                if(visited[m+N[k]]==1e8) {
                  visited[m+N[k]] = iter+1;
                  u[m+N[k]]       = u[m]+Nd[k];
                } else {
                  float temp = u[m]+Nd[k];
                  if(u[m+N[k]]>temp)
                    u[m+N[k]] = temp;
                }
                if(u[m+N[k]]>maxDistance) {
                  maxDistance = u[m+N[k]];
                }
              }
            }
          }
        }
      }
      if(maxDistance2==maxDistance)
        break;
    }
  }

  // COMPUTATION OF EVOLUTION FOR NEGATIVE PART
  if(MaxDistanceNegative>0.) {
    // COMPUTATION OF THE DISTANCE FUNCTION IN AN ITERATIVE WAY
    float maxDistance = 0.;
    int iter          = 0;
    while(maxDistance<MaxDistanceNegative) {
      iter--;
      float maxDistance2 = maxDistance;
      printf("%1.2lf ",maxDistance);
      #ifdef _OPENMP
      #pragma omp parallel for
      #endif
      for(int z=1; z<nSlices-1; z++) {
        for(int x=1; x<width-1; x++) {
          for(int y=1; y<height-1; y++) {
            int m=x+y*width+z*size2d;
            if(Is[m]>0 || visited[m]!=iter)
              continue;
            for(unsigned int k=0; k<N.size(); k++) {
              if(Is[m+N[k]]<=0) {
                if(visited[m+N[k]]==1e8) {
                  visited[m+N[k]] = iter-1;
                  u[m+N[k]]       = u[m]-Nd[k];
                } else {
                  float temp = u[m]-Nd[k];
                  if(u[m+N[k]]<temp)
                    u[m+N[k]] = temp;
                }
                if(-u[m+N[k]]>maxDistance) {
                  maxDistance = -u[m+N[k]];
                }
              }
            }
          }
        }
      }
      if(maxDistance2==maxDistance)
        break;
    }
  }


  // BOUNDARY
  for(int z=1; z<nSlices-1; z++) {
    int m = z*size2d;
    for(int x=0; x<width; x++) {
      u[m+x]                  = u[m+x+width];
      u[m+x+(height-1)*width] = u[m+x+(height-2)*width];
    }
    for(int y=1; y<height-1; y++) {
      u[m+y*width]         = u[m+y*width+1];
      u[m+y*width+width-1] = u[m+y*width+width-2];
    }
  }
  for(int x=0; x<width; x++) {
    for(int y=1; y<height-1; y++) {
      u[x+y*width]                    = u[x+y*width+size2d];
      u[x+y*width+(nSlices-1)*size2d] = u[x+y*width+(nSlices-2)*size2d];
    }
  }

  printf(" end distance_function3D()\n");
  return u;
}


#endif
