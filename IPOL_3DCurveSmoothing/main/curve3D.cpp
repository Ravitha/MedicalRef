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
 * \file curve3D.cpp
 * \brief Methods to manage 3D curves
 * \author Luis Alvarez and Daniel Santana-Cedrés
 */


#include "curve3D.h"
#include <math.h>

//#define DEBUG

/** \brief Function to define a neighborhood index vector for an image3D
 *
 * \param[in] width: image width
 * \param[in] size2d: size of an axial slice
 * \param[in] type: neighborhood type
 * \return Returns a vector with the neighborhood indexes
 *
 */
vector<int> Neighborhood(int width,int size2d,int type)
{
  vector<int> N(6);
  int k=0;
  N[k++] = 1;                   //0
  N[k++] = -1;                  //1
  N[k++] = width;               //2
  N[k++] = -width;              //3
  N[k++] = size2d;              //4
  N[k++] = -size2d;             //5
  if(type==0)
    return N;
  N.push_back(1+width);         //6
  N.push_back(1-width);         //7
  N.push_back(-1+width);        //8
  N.push_back(-1-width);        //9
  N.push_back(1+size2d);        //10
  N.push_back(1-size2d);        //11
  N.push_back(-1+size2d);       //12
  N.push_back(-1-size2d);       //13
  N.push_back(width+size2d);    //14
  N.push_back(width-size2d);    //15
  N.push_back(-width+size2d);   //16
  N.push_back(-width-size2d);   //17
  if(type==1)
    return N;
  N.push_back(1+width+size2d);  //18
  N.push_back(1+width-size2d);  //19
  N.push_back(1-width+size2d);  //20
  N.push_back(1-width-size2d);  //21
  N.push_back(-1+width+size2d); //22
  N.push_back(-1+width-size2d); //23
  N.push_back(-1-width+size2d); //24
  N.push_back(-1-width-size2d); //25

  return N;
}

//==============================================================================

/** \brief Function to compute the distance in mm to the neighbors according to
 *         the voxel size
 *
 * \param[in] vx: voxel size in the x dimension
 * \param[in] vy: voxel size in the y dimension
 * \param[in] vz: voxel size in the z dimension
 * \param[in] type: neighborhood type
 * \return Returns a vector with the distances computed
 *
 */
vector<double> Distace2Neighbors(const double &vx,const double &vy,
                                 const double &vz,const int &type)
{
  vector<double> N(6);
  int k=0;
  N[k++] = vx;       //0
  N[k++] = vx;       //1
  N[k++] = vy;       //2
  N[k++] = vy;       //3
  N[k++] = vz;       //4
  N[k++] = vz;       //5
  if(type==0)
    return N;
  double vxy = sqrt(vx*vx+vy*vy);
  double vxz = sqrt(vx*vx+vz*vz);
  double vyz = sqrt(vz*vz+vy*vy);
  N.push_back(vxy);  //6
  N.push_back(vxy);  //7
  N.push_back(vxy);  //8
  N.push_back(vxy);  //9
  N.push_back(vxz);  //10
  N.push_back(vxz);  //11
  N.push_back(vxz);  //12
  N.push_back(vxz);  //13
  N.push_back(vyz);  //14
  N.push_back(vyz);  //15
  N.push_back(vyz);  //16
  N.push_back(vyz);  //17
  if(type==1)
    return N;
  double vxyz = sqrt(vx*vx+vy*vy+vz*vz);
  N.push_back(vxyz); //18
  N.push_back(vxyz); //19
  N.push_back(vxyz); //20
  N.push_back(vxyz); //21
  N.push_back(vxyz); //22
  N.push_back(vxyz); //23
  N.push_back(vxyz); //24
  N.push_back(vxyz); //25
  #ifdef DEBUG
    for(unsigned int k=0; k<N.size(); k++)
      printf("N[%d]=%lf\n",k,N[k]);
  #endif // DEBUG
  return N;
}

//==============================================================================

/** \brief Procedure to change neighborhood order to go through the neighborhood
 *         in a different order
 *
 * \param[in] N:
 * \param[in] Nd:
 * \return Updates the neighborhood order
 *
 */
void change_order_neighborhood(vector<int> &N,vector<double> &Nd)
{
  if(N.size()<6 || N.size()!=Nd.size())
    return;
  for(int k=0; k<2; k++) {
    int ti    = N[k];
    double td = Nd[k];
    N[k]      = N[k+4];
    Nd[k]     = Nd[k+4];
    N[k+4]    = ti;
    Nd[k+4]   = td;
  }
  if(N.size()<18)
    return;
  for(int k=6; k<10; k++) {
    int ti    = N[k];
    double td = Nd[k];
    N[k]      = N[k+4];
    Nd[k]     = Nd[k+4];
    N[k+4]    = ti;
    Nd[k+4]   = td;
  }
}

//==============================================================================

/** \brief Computation of a point in a segment with a prefixed distance to a
 *         given point. This function is called when
 *         (q-p0).norm()<d && (q-p1).norm()>d
 *
 * \param[in] q: Prefixed point to evaluate the distance
 * \param[in] p0,p1: Segment points
 * \param[in] d: Distance from 'q' to the point on the segment
 * \return Returns the point at the prefixed distance
 *
 */
point3d<double> fix_distance_point2segment(point3d<double> q,
                                           point3d<double> p0,
                                           point3d<double> p1, double d)
{
  point3d<double> v = p1-p0;
  point3d<double> p = p0-q;
  double a          = v.by(v);
  double b          = p.by(v);
  double c          = -d*d+p.by(p);
  double det        = b*b-a*c;
  if(det<0 || a==0)
    return point3d<double>();
  det       = sqrt(det);
  double r0 = (-b+det)/(a);
  return p0+v*r0;
}

//==============================================================================

/** \brief curve3D reparametrization with fixed length segments. The length used
*          is the min of the voxel size.
 *
 * \param[in] Cl: input curve3D
 * \param[in] length_segment: length of segment
 * \param[in] max_length: max length of the new curve
 * \return It returns the new parametrization. We do not change the initial and
 *         final points of the original curve3D
 *
 */
vector <point3d<double> > curve3D_reparameterization(
                                                  vector <point3d<double> > &Cl,
                                                  const double length_segment,
                                                  const double max_length)
{
  #ifdef DEBUG
    printf("Cl.size()=%d\n",(int)Cl.size());
  #endif // DEBUG

  // New parametrization
  vector <point3d<double> > Cn;

  // Adding points to the new parametrization
  Cn.push_back(Cl[0]);
  unsigned int k = 1;
  while(true) {
    if(Cn.size()*length_segment>=max_length)
      break;
    point3d<double> v = Cl[k]-Cn[Cn.size()-1];
    double norm       = v.norm();
    if(norm>=length_segment) {
      // Adding points in the segment Cn[Cn.size()-1],Cl[k]
      Cn.push_back(Cn[Cn.size()-1]+v*(length_segment/norm));
      continue;
    }
    k++;
    // Looking for points Cl[k] such that |Cn[Cn.size()-1]-Cl[k]|>length_segment
    while(k<Cl.size()&& norm<length_segment) {
      point3d<double> v = Cl[k]-Cn[Cn.size()-1];
      norm              = v.norm();
      if(norm<length_segment)
        k++;
    }
    if(k==Cl.size())
      break;
    // Adding a point in the segment Cl[k-1],Cl[k]
    Cn.push_back(fix_distance_point2segment(Cn[Cn.size()-1],Cl[k-1],Cl[k],
                                            length_segment));

  }
  if( (Cl[Cl.size()-1]-Cn[Cn.size()-1]).norm()>0.)
    Cn.push_back(Cl[Cl.size()-1]);
  return Cn;
}

//==============================================================================

/** \brief curve3D smoothing. The curve3D points tend to move to points
 *         maximizing the distance to the segmentation contours but at, the same
 *         time, it keeps the curve regular. We use the curvature scheme
 *
 * \param[in,out] Cl: Input/output curve3D
 * \param[in] u: Distance function to the segmentation contour
 * \param[in] w: Weight to balance the smoothing and maximum distance terms
 * \param[in] MaxIter: Maximum number of iterations
 * \param[in] TOL: Tolerance to stop iterations
 * \param[in] GT: Optional parameter to compute the error between the curve we
 *            smooth and its ground-truth (both curves must have the same scale)
 * \return Updates Cl with the smoothed version of the curve
 *
 */
void curve3D_smoothing(vector <point3d<double> > &Cl, const image3D<float> &u,
                       float w, int MaxIter, float TOL,
                       vector <point3d<double> > *GT)
{
  point3d<double> v = u.voxel_size();

  #ifdef DEBUG
    FILE *f = NULL;
    if(GT!=NULL) {
      char buf[100];
      sprintf(buf,"curve_evolution_w_%.2lf.csv",w);
      f = fopen(buf,"w");
      fprintf(f,"w;iteration;length;energy;error\n");
    }
  #endif // DEBUG

  // Computation of the min voxel size
  double v_min = v.x<v.z?(v.x<v.y?v.x:v.y):(v.z<v.y?v.z:v.y);

  // Reparameterization of the curve3D and update of variables
  Cl = curve3D_reparameterization(Cl,1.,1e20);

  double length0 = length(Cl);

  // Computation of the initial energy
  double u_sum = 0;
  for(unsigned int k=1; k<Cl.size(); k++) {
    u_sum -= u(Cl[k]/v);
  }
  double energy0 = w*length0+u_sum;

  #ifdef DEBUG
    if(GT!=NULL) {
      fprintf(f,"%lf;%d;%lf;%lf;%lf\n",w,0,length0,energy0,
              sqrt(distance_curve2curve(*GT,Cl)));
    }
  #endif // DEBUG

  // Computation of the maximum step size
  double max_step = v_min*v_min/(M_PI*w);
  if(max_step>0.5*v_min)
    max_step = 0.5*v_min;
  printf("max_step=%lf\n",max_step);
  float step = max_step;
  step      *= 0.5;

  int iter = 0;
  while(iter++<MaxIter) {
    // Computation of the image gradient
    vector< point3d<double> > Gu(Cl.size()); // gradient
    Gu[Gu.size()-1] = Gu[0] = point3d<double>(0.,0.,0.);
    for(unsigned int k=1; k<Gu.size()-1; k++) {
      point3d<double> q = Cl[k]/v;
      double Eu = u(q.x,q.y,q.z);
      Gu[k].x   = (u(q.x+1,q.y,q.z)-Eu)/v.x;
      Gu[k].y   = (u(q.x,q.y+1,q.z)-Eu)/v.y;
      Gu[k].z   = (u(q.x,q.y,q.z+1)-Eu)/v.z;
    }
    // Computation of the normal to the curve and the curvature
    vector< point3d<double> > normal(Cl.size()); // Normal unit vector
    vector<double> curv(Cl.size(),0.); // Curvature
    for(unsigned int k=1; k<Cl.size()-1; k++) {
      normal[k]   = (Cl[k-1]+Cl[k+1])*0.5-Cl[k];
      double norm = normal[k].norm();
      if(norm>0.) {
        normal[k].x /= norm;
        normal[k].y /= norm;
        normal[k].z /= norm;
      }
      double temp = norm/v_min;
      if(temp<=1.) {
        curv[k] = (M_PI-2.*acos(temp))/v_min;
      }
    }

    double length1 = length0;
    // Computation of the new curve
    vector< point3d<double> > Cln(Cl.size());
    Cln[0]           = Cl[0];
    Cln[Cl.size()-1] = Cl[Cl.size()-1];
    for(unsigned int k=1; k<Cl.size()-1; k++) {
      Cln[k] = Cl[k]+Gu[k]*step+normal[k]*curv[k]*w*step;
    }
    // Computation of the new energy
    Cln   = curve3D_reparameterization(Cln,1.,1e20);
    u_sum = 0;
    for(unsigned int k=1; k<Cln.size(); k++) {
      u_sum -= u(Cln[k]/v);
    }
    length1  = length(Cln);
    double energy = w*length1+u_sum;
    #ifdef DEBUG
      printf("iter=%d step=%e w*length(Cln)=%lf u_sum=%lf energie=%lf\n",
             iter,step,w*length1,u_sum,energy);
      if(GT!=NULL) {
        fprintf(f,"%lf;%d;%lf;%lf;%lf\n",w,iter,length1,energy,
                sqrt(distance_curve2curve(*GT,Cln)));
      }
    #endif // DEBUG

    Cl = Cln;
    printf("((energy0-energy)/(fabs(energy0)+1e-20))=%e\n",
           ((energy0-energy)/(fabs(energy0)+1e-20)));
    if(GT!=NULL)
      printf("Error(GT,Cln) = %lf\n",sqrt(distance_curve2curve(*GT,Cln)));
    if( (fabs(energy0-energy)/(fabs(energy0)+1e-20))<TOL)
      iter  = MaxIter+1;

    energy0 = energy;
  }
  #ifdef DEBUG
    if(GT!=NULL)
      fclose(f);
  #endif // DEBUG
}

//==============================================================================

/** \brief Function to compute of the curve3D length between 2 points
 *
 * \param[in] Cl: curve3D
 * \param[in] k0,k1: indexes of the points inside the curve
 * \return Returns the distance between the two points
 *
 */
double length(const vector <point3d<double> > &Cl, int k0, int k1)
{
  if(k1==9999999 || k1>=(int)Cl.size())
    k1 = (int)Cl.size()-1;
  if(k0==9999999 || k0>=(int)Cl.size())
    k0 = (int)Cl.size()-1;
  double sum=0.;
  if(k1>k0) {
    for(int k=k0+1; k<=k1; k++)
      sum += (Cl[k]-Cl[k-1]).norm();
  } else {
    for(int k=k1+1; k<=k0; k++)
      sum += (Cl[k]-Cl[k-1]).norm();
  }
  return sum;
}

//==============================================================================

/** \brief Function to read a curve3D from disk
 *
 * \param[in] filename: name of the file with the curve points
 * \return Returns the vector containing the curve points
 *
 */
vector <point3d<double> > read_curve3D(char filename[400])
{

  // Open the file
  FILE *f;
  f = fopen (filename, "r");
  if (f == NULL) {
    printf("Error opening file %s\n",filename);
    return vector <point3d<double> >();
  }

  // Reading 3D vector
  int N;

  if (feof(f))
    return vector <point3d<double> >();
  int val = fscanf(f,"%d\n",&N);
  if(val==EOF)
    printf("Error reading 3D vector\n");
  #ifdef DEBUG
    printf("N=%d\n",N);
  #endif // DEBUG

  // Creating a 3D curve
  vector <point3d<double> > C(N);

  // Reading the 3D points coordinates
  double x,y,z;
  for(int k=0; k<N; k++) {
    if (feof(f))
      return vector <point3d<double> >();
    int val = fscanf(f,"%lf %lf %lf\n",&x,&y,&z);
    if(val==EOF)
      printf("Error reading 3D point coordinate\n");
    #ifdef DEBUG
      printf("x=%lf y=%lf z=%lf\n",x,y,z);
    #endif // DEBUG
    C[k]=point3d<double>(x,y,z);
  }
  fclose(f);

  return C;
}

//==============================================================================

/** \brief Procedure to write a 3D curve in a file
 *
 * \param[in] C: vector with the curve points
 * \param[in] filename: name of the file in which the points will be stored
 * \return Writes in disk the curve points
 *
 */
void write_curve3D(vector <point3d<double> > &C, char filename[400])
{
  // Open the file
  FILE *f;
  f = fopen (filename, "w");
  if (f == NULL) {
    printf("Error opening file %s\n",filename);
    return;
  }
  fprintf(f,"%d\n",(int) C.size());
  for(int k=0; k<(int) C.size(); k++) {
    fprintf(f,"%lf %lf %lf\n",C[k].x,C[k].y,C[k].z);
  }
  fclose(f);
}

//==============================================================================

/** \brief Function to compute the distance between a 3D point and a segment
 *
 * \param[in] q: 3D point from which we compute the distance
 * \param[in] p1,p2: 3D points which describe the segment
 * \return Returns the distance between q and the segment described by p1 and p2
 *
 */
double distace_point3d2segment(point3d<double> q,point3d<double> p1,
                               point3d<double> p2)
{
  point3d<double> v = p2-p1;
  double lambda=-(p1.x*v.x+p1.y*v.y-q.x*v.x+p1.z*v.z-q.y*v.y-q.z*v.z)/v.norm2();
  if(lambda<0) return (p1-q).norm2();
  if(lambda>1) return (p2-q).norm2();
  return (p1+v*lambda-q).norm2();
}

//==============================================================================

/** \brief  Function to compute the distance between a 3D point to a curve
            taking into account the distance to the segments of the curve
 *
 * \param[in] p: 3D point from which we compute the distance
 * \param[in] c: reference to a vector of 3D points which describe the curve
 * \param[in,out] kmin: index of the curve to compute the distance
 * \return Returns the distance between the point and the curve
 *
 */
double distance_point3d2curve(point3d<double> &p,vector< point3d<double> > &c,
                              int &kmin)
{
  if (c.size()==0) return 1e20;
  kmin = 0;
  double d1 = 1e20, d2 = 1e20, dmin = (p-c[0]).norm2();
  for(int k=1; k<(int)c.size(); k++) {
    double aux = (p-c[k]).norm2();
    if(aux<dmin) {
      dmin = aux;
      kmin = k;
    }
  }

  if(kmin>0) {
   d1 = distace_point3d2segment(p,c[kmin],c[kmin-1]);
  }
  if(kmin<(int)c.size()-1) {
   d2 = distace_point3d2segment(p,c[kmin],c[kmin+1]);
  }

  return d1<d2?d1:d2;
}

//==============================================================================

/** \brief Function to compute the average distance between two curves
 *
 * \param[in] c1: vector of 3D points which represents the first curve
 * \param[in] c2: vector of 3D points which represents the second curve
 * \return Returns the average distance between the provided curves
 *
 */
double distance_curve2curve(vector< point3d<double> > &c1,
                            vector< point3d<double> > &c2)
{
  if(c1.size()<2 || c2.size()<2) return 1e20;
  int k2, k1 = c1.size()/2;
  double distance1 = distance_point3d2curve(c1[k1],c2,k2);
  int Np1 = 1;
  for(int k=k1+1; k<(int)c1.size() && k2<(int)c2.size()-1; k++) {
    distance1 += distance_point3d2curve(c1[k],c2,k2);
    Np1++;
  }
  for(int k=k1-1; k>=0 && k2>0; k--) {
    distance1 += distance_point3d2curve(c1[k],c2,k2);
    Np1++;
  }

  k2 = c2.size()/2;
  double distance2 = distance_point3d2curve(c2[k2],c1,k1);
  int Np2 = 1;
  for(int k=k2+1; k<(int)c2.size() && k1<(int)c1.size()-1; k++) {
    distance2 += distance_point3d2curve(c2[k],c1,k1);
    Np2++;
  }
  for(int k=k2-1; k>=0 && k1>0; k--) {
    distance2 += distance_point3d2curve(c2[k],c1,k1);
    Np2++;
  }

  return (distance1/Np1+distance1/Np2)*0.5;
}
