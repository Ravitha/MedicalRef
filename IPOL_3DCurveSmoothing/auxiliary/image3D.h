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
 * \file image3D.h
 * \brief Class  to store 3D image
 * \author Luis Alvarez
*/
#ifndef image3D_H
#define image3D_H

#include <fstream>
#include <cstring>
#include <stdlib.h>

#include "dbh.h"
#include "point3d.h"

//#define IMAGE3D_DEBUG
///\brief Macro to control unsigned short type
#define USHORT 0
///\brief Macro to control unsigned char type
#define UCHAR 1
///\brief Macro to control short type
#define SHORT 2
///\brief Macro ton control float type
#define FLOAT 3

using namespace std;

//==============================================================================

/**
 * \class image3D
 * \brief Class  to store 3D images in Analyze format and basic methods
 * \param[in] T: the data type of the image
 * \author Luis Alvarez
 */
template <class  T>
class image3D
{
  /// \brief std vector to allocate image3D "slice 1 + slice 2 + ..."
  std::vector <T> image3D_;
  /// \brief image width
  int width_;
  /// \brief image height
  int height_;
  /// \brief number of slices of the image
  int nSlices_;
  /// \brief voxel size in x
  double voxel_size_x_;
  /// \brief voxel size in y
  double voxel_size_y_;
  /// \brief voxel size in z
  double voxel_size_z_;
  /// \brief image 2D size == width*height
  int size2d_;

public:

  // CONSTRUCTORS - DESTRUCTOR METHODS
  image3D();
  image3D(const int width,const int height, const int nSlices,
          const float voxel_size_x,const float voxel_size_y,
          const float voxel_size_z,const T &a);
  image3D(const int width,const int height, const int nSlices,
          const float voxel_size_x,const float voxel_size_y,
          const float voxel_size_z);
  image3D(const char *fileName, const char *fileName_2);
  ~image3D();
  /** \brief Procedure to clear an image3D object
   *
   * \return Empties the image vector and sets its size to zero
   *
   */
  void clear()
  {
    image3D_.clear();
    image3D_.shrink_to_fit();
    width_ = height_ = nSlices_ = 0;
  };
  /** \brief Procedure to resize the number of slices
   *
   * \param[in] NslicesNew: the new number of slices
   * \return If the new number of slices is lower than the current one, it
   *         returns the resized image. Otherwise, the call does not affect the
   *         image
   *
   */

  void resize(const int &NslicesNew)
  {
    if(NslicesNew<nSlices_) {
      nSlices_ = NslicesNew;
      image3D_.resize(width_*height_*nSlices_);
    }
  }
//------------------------------------------------------------------------------

  // class ELEMENT ACCESS METHODS
  ///\brief Method to get image width
  inline int width() const
  {
    return width_;
  }
  ///\brief Method to get image height
  inline int height() const
  {
    return height_;
  }
  ///\brief Method to get the number of slices of an image
  inline int nSlices() const
  {
    return nSlices_;
  }
  ///\brief Method to get the size of the voxel in X
  inline float voxel_size_x() const
  {
    return voxel_size_x_;
  }
  ///\brief Method to get the size of the voxel in Y
  inline float voxel_size_y() const
  {
    return voxel_size_y_;
  }
  ///\brief Method to get the size of the voxel in Z
  inline float voxel_size_z() const
  {
    return voxel_size_z_;
  }
  ///\brief Method to get the complete size of an image (number of elements)
  inline int size() const
  {
    return image3D_.size();
  }
  ///\brief Method to get the size of a slice (number of elements)
  inline int size2d() const
  {
    return size2d_;
  }
  inline T& operator[](const int &i);
  inline T& operator()(const int &x,const int &y,const int &z);
  inline T& operator()(const int &pos,const int &dx,const int &dy,
                       const int &dz);
  inline double operator()(const double &x,const double &y,const double &z,
                           int boundary_type=0) const;
  /** \brief Method to get the image value in a point
   *
   * \param[in] p: point which represents the position in the image
   * \return Returns the image value in p
   *
   */
  inline double operator()(const point3d<double> &p) const
  {
    return (*this)(p.x,p.y,p.z);
  }
  inline const T& operator[](const int &i) const;
  /** \brief Method to get voxel size
   *
   * \return Returns the voxel size in X, Y, Z in a 3D point format
   *
   */
  inline point3d<double> voxel_size() const
  {
    return point3d<double>(voxel_size_x_,voxel_size_y_,voxel_size_z_);
  }
  /** \brief Method to get image dimensions
   *
   * \return Returns the with, height, and number of slices in a 3D point format
   *
   */
  inline point3d<double> image_size() const
  {
    return point3d<double>((double)width_,(double)height_,(double) nSlices_);
  }
  /** \brief Method to set the voxel size
   *
   * \param[in] vx: voxel size in X
   * \param[in] vy: voxel size in Y
   * \param[in] vz: voxel size in Z
   * \return Sets the voxel size to the 3D image object
   *
   */
  void set_voxel_size(const double &vx,const double &vy,const double &vz)
  {
    voxel_size_x_ = vx;
    voxel_size_y_ = vy;
    voxel_size_z_ = vz;
  }
  /** \brief Method to print information regarding the image
   *
   * \return Shows in the screen image information: dimensions, voxel size,
   *         complete size, and number of voxels
   *
   */
  void print()
  {
    printf("dim: %d,%d,%d voxel size: %lf %lf %lf size=%d n.voxels=%d\n",width_,
           height_,nSlices_,voxel_size_x_,voxel_size_y_,voxel_size_z_,
           width_*height_*nSlices_,image3D_.size());
  }
//------------------------------------------------------------------------------

  // BASIC OPERATION METHODS
  // write in HDR format
  int write_analyze(const char *name,const int type=SHORT) const;
  // write in ASCII format
  int write_ascii(const char *name) const;
  // operator =
  ///\brief Assignment operator
  template <class  U> image3D<T> & operator=(const image3D<U> &image3D2);
  // copy
  void copy(const image3D<T> &image3D2);

};

//==============================================================================

/**
 * \brief operator [] to access image value (const)
 * \param[in] i: position in the image vector from which we want to get the
 *            value
 * \return Returns the image value in such a position
 * \author Luis Alvarez
*/
template <class  T>
const T& image3D<T>::operator[](const int &i) const
{
  if(i<0 || i>=(int) (image3D_.size())) {
    printf("image<T>: bounds error vector access\n");
    printf("image size()=%d index to access the vector=%d\n",
           (int)image3D_.size(),i);
    int j;
    int val = scanf("%d",&j);
    if(val==EOF)
      printf("Error reading value from stdin\n");
    exit(EXIT_FAILURE);
  }
  return image3D_[i];
}

//==============================================================================

/**
 * \brief operator [] to access image value
 * \param[in] i: position in the image vector from which we want to get the
 *            value
 * \return Returns the image value in such a position
 * \author Luis Alvarez
*/
template <class  T>
T& image3D<T>::operator[](const int &i)
{
  if(i<0 || i>=(int) (image3D_.size())) {
    printf("image<T>: bounds error vector access\n");
    printf("image size()=%d index to access the vector=%d\n",
           (int)image3D_.size(),i);
    int j;
    int val = scanf("%d",&j);
    if(val==EOF)
      printf("Error reading value from stdin\n");
    exit(EXIT_FAILURE);
  }
  return  (image3D_.at(i));
}

//==============================================================================

/**
 * \brief operator () to access image value
 * \param[in] x: x coordinate of the image point
 * \param[in] y: y coordinate of the image point
 * \param[in] z: z coordinate of the image point
 * \return Returns the image value associated to the given 3D position
 * \author Luis Alvarez
*/
template <class  T>
T& image3D<T>::operator()(const int &x,const int &y,const int &z)
{
  unsigned int value=x+y*width_+z*size2d_;
  if(value>=(image3D_.size()) || value<0) {
    printf("image<T>: bounds error vector access\n");
    printf("image size()=%d index to access the vector=%d\n",
           (int)image3D_.size(),value);
    int j;
    int val = scanf("%d",&j);
    if(val==EOF)
      printf("Error reading value from stdin\n");
    exit(EXIT_FAILURE);
  }
  return image3D_[value];
}

//==============================================================================

/**
 * \brief operator () to access image value using a neighborhood
 * \param[in] pos: position in the image vector
 * \param[in] dx: displacement in x
 * \param[in] dy: displacement in y
 * \param[in] dz: displacement in z
 * \return Returns the value in the position of the neighborhood indicated by
 *         (dx,dy,dz)
 * \author Luis Alvarez
*/
template <class  T>
T& image3D<T>::operator()(const int &pos,const int &dx,const int &dy,
                          const int &dz)
{
  unsigned int value=pos+dx+dy*width_+dz*size2d_;
  if(value>=(image3D_.size()) || value<0) {
    printf("image<T>: bounds error vector access\n");
    printf("image size()=%d index to access the vector =%d\n",
           image3D_.size(),value);
    int j;
    scanf("%d",&j);
    exit(EXIT_FAILURE);
  }
  return image3D_[value];
}

//==============================================================================

/**
 * \brief operator () to access image value using bilinear interpolation
 * \param[in] x: x coordinate of the image point
 * \param[in] y: y coordinate of the image point
 * \param[in] z: z coordinate of the image point
 * \param[in] boundary_type: neighborhood type
 * \return Returns the interpolated value
 * \author Luis Alvarez
*/
template <class  T>
double image3D<T>::operator()(const double &x,const double &y,const double &z,
                              const int boundary_type) const
{
  if(nSlices_==1) {
    int ix=(int) x;
    int iy=(int) y;
    if(boundary_type==1 && (ix>=width_ || iy>=height_ || ix<-1|| iy<-1))
      return (-3000.);
    else {
      if(ix>=width_)
        ix=width_-1;
      else if(ix<-1)
        ix=0;
      if(iy>=height_)
        iy=height_-1;
      else if(iy<-1)
        iy=0;
    }
    double dx=x-ix;
    double dy=y-iy;
    int ix1=(ix<(width_-1))?ix+1:ix;
    int iy1=(iy<(height_-1))?iy+1:iy;
    if(ix==-1)
      ix=0;
    if(iy==-1)
      iy=0;
    return(
            (1-dx)*(1-dy)*(*this)[ix+iy*width_]+
            (dx)*(1-dy)*(*this)[ix1+iy*width_]+
            (dx)*(dy)*(*this)[ix1+iy1*width_]+
            (1-dx)*(dy)*(*this)[ix+iy1*width_]
          );
  }
  int ix=(int) x;
  int iy=(int) y;
  int iz=(int) z;
  if(boundary_type==1 && (ix>=width_ || iy>=height_ || ix<-1|| iy<-1 ||
                          iz>=nSlices_ || iz<-1))
    return (-3000.);
  else {
    if(ix>=width_)
      ix=width_-1;
    else if(ix<-1)
      ix=0;
    if(iy>=height_)
      iy=height_-1;
    else if(iy<-1)
      iy=0;
    if(iz>=nSlices_)
      iz=nSlices_-1;
    else if(iz<-1)
      iz=0;
  }
  double dx=x-ix;
  double dy=y-iy;
  double dz=z-iz;
  int ix1=(ix<(width_-1))?ix+1:ix;
  int iy1=(iy<(height_-1))?iy+1:iy;
  int iz1=(iz<(nSlices_-1))?iz+1:iz;
  if(ix==-1)
    ix=0;
  if(iy==-1)
    iy=0;
  if(iz==-1)
    iz=0;

  int p000=ix+iy*width_+iz*size2d_;
  int p100=p000,p010=p000,p001=p000,p110=p000,p101=p000,p011=p000,p111=p000;
  if(ix1>ix) {
    p100++;
    p110++;
    p101++;
    p111++;
  }
  if(iy1>iy) {
    p010+=width_;
    p110+=width_;
    p011+=width_;
    p111+=width_;
  }
  if(iz1>iz) {
    p001+=size2d_;
    p101+=size2d_;
    p011+=size2d_;
    p111+=size2d_;
  }

  return(
          (1-dx)*(1-dy)*(1-dz)*image3D_[p000]+
          dx*(1-dy)*(1-dz)*image3D_[p100]+
          dx*dy*(1-dz)*image3D_[p110]+
          dx*dy*dz*image3D_[p111]+
          dx*(1-dy)*dz*image3D_[p101]+
          (1-dx)*dy*(1-dz)*image3D_[p010]+
          (1-dx)*dy*dz*image3D_[p011]+
          (1-dx)*(1-dy)*dz*image3D_[p001]
        );

}

//==============================================================================

/**
 * \brief Constructor from an image file
 * \param [in] fileName Input filename.
 * \author Luis Alvarez
*/
template <class  T>
image3D<T>::image3D(const char *fileName, const char *fileName_2)
{
  // WE CHECK IF THE FILE HAS THE .hdr EXTENSION
  int N1=strlen(fileName);
  printf("N1=%d\n",N1);
  if( (fileName[N1-1]!='r' && fileName[N1-1]!='R') ||
      (fileName[N1-2]!='d' && fileName[N1-2]!='D') ||
      (fileName[N1-3]!='h' && fileName[N1-3]!='H')
    ) { // case of ASCII FILE
    // OPEN THE FILE
    FILE *f;
    f = fopen (fileName, "r");
    if (f == NULL) {
      printf("Error opening file %s\n",fileName);
      return;
    }
    int val = fscanf(f,"%d %d %d\n",&width_,&height_,&nSlices_);
    if(val==EOF)
      printf("Error reading width, height and number of slices\n");
    val = fscanf(f,"%lf %lf %lf\n",&voxel_size_x_,&voxel_size_y_,&voxel_size_z_);
    if(val==EOF)
      printf("Error reading voxel size\n");

    printf("%d %d %d\n",width_,height_,nSlices_);
    printf("%lf %lf %lf\n",voxel_size_x_,voxel_size_y_,voxel_size_z_);

    double Size=(double) width_*height_*nSlices_;
    if(Size<=0 || Size> 268435456.) {
      *this=image3D();
      return;
    }
    image3D_.resize(width_*height_*nSlices_);
    size2d_ = width_*height_;
    for(unsigned int k=0; k<image3D_.size(); k++) {
      if (feof(f)) {
        *this=image3D();
        return;
      }
      int temp;
      val = fscanf(f,"%d ",&temp);
      if(val==EOF)
        printf("Error reading temp\n");
      image3D_[k]=(T) temp;
    }
    fclose(f);
    return;
  }


  ifstream filein;
  filein.open(fileName, ios::binary);

  if(filein.fail()) {
    cout << "We can not open file : " << fileName << "\n";
    int j;
    int val = scanf("%d",&j);
    if(val==EOF)
      printf("Error reading value from stdin\n");
    exit(EXIT_FAILURE);
  }


  struct dsr hdr;
  filein.read((char *)(&hdr),sizeof(dsr)); //Read Analyze header
  filein.close();


  int size;
  if (hdr.hk.sizeof_hdr == 348) {
    width_            = hdr.dime.dim[1];    // Number of rows
    height_           = hdr.dime.dim[2];    // Number of columns
    nSlices_          = hdr.dime.dim[3];    // Number of slices
    voxel_size_x_     = hdr.dime.pixdim[1]; // Voxel size in X
    voxel_size_y_     = hdr.dime.pixdim[2]; // Voxel size in Y
    voxel_size_z_     = hdr.dime.pixdim[3]; // Voxel size in Z
    printf("dimension = (%d,%d,%d), voxel_size = (%lf,%lf,%lf)\n",
           width_,height_,nSlices_,voxel_size_x_,voxel_size_y_,voxel_size_z_);
    size2d_           = width_*height_;
    size              = size2d_*nSlices_;

    if (size == 0) {
      printf("Error reading file %s. Size = %d", fileName, size);
      exit(EXIT_FAILURE);
    }
    image3D_.resize(size);

  } else {
    cout << "Invalid Analyze format" << endl;
    exit(EXIT_FAILURE);
  }

  // READING RAW DATA
  string rut = fileName_2;
  //int n      = rut.size();
  //rut.erase(n-3,3);
  //rut.append("img");

  if(hdr.dime.datatype == DT_UNSIGNED_CHAR) {
    fflush(stdout);
    ifstream fichero;
    fichero.open(rut.c_str(), ios::binary);
    unsigned char temp[size2d_];
    for(int k=0; k<size; k+=size2d_) {
      fichero.read((char *) temp, size2d_*sizeof(unsigned char));
      for (int n=0; n<size2d_; n++)
        image3D_[k+n]=(T) temp[n];
    }
    fichero.close();

  } else if(hdr.dime.datatype == DT_SIGNED_SHORT) {
    ifstream fichero;
    fichero.open(rut.c_str(), ios::binary);

    for(int k=0; k<size; k++) {
      short temp;
      fichero.read((char *)(&temp), sizeof(short));
      image3D_[k]=(T) temp;
    }
    fichero.close();
  } else if(hdr.dime.datatype == DT_FLOAT) {
    fflush(stdout);
    ifstream fichero;
    fichero.open(rut.c_str(), ios::binary);

    for(int k=0; k<size; k++) {
      float temp;
      fichero.read((char *)(&temp), sizeof(float));
      image3D_[k]=(T) temp;
    }
    fichero.close();

  } else { // BY DEFAULT WE ASSUME THAT THE IMAGE IS short
    ifstream fichero;
    fichero.open(rut.c_str(), ios::binary);
    short *temp = new short[size2d_];
    for(int k=0; k<size; k+=size2d_) {
      fichero.read((char *) temp, size2d_*sizeof(short));
      for (int n=0; n<size2d_; n++)
        image3D_[k+n]=(T) temp[n];
    }
    delete[] temp;
    fichero.close();
  }
}

//==============================================================================

/**
 * \brief Function to write an image to disk in ASCII format
 * \param[in] fileName: name of the image.
 * \return Returns 0 if the writing operation ends successfully, -1 if the image
 *         size is 0, or -2 if the opening operation fails.
 * \author Luis Alvarez
*/
template <class  T>
int image3D<T>::write_ascii(const char *fileName /** image file name */) const
{
  if(image3D_.size()==0)
    return -1;
  // OPEN THE FILE
  FILE *f;
  f = fopen (fileName, "w");
  if (f == NULL) {
    printf("Error opening file %s\n",fileName);
    return -2;
  }


  fprintf(f,"%d %d %d\n",width_,height_,nSlices_);
  fprintf(f,"%lf %lf %lf\n",voxel_size_x_,voxel_size_y_,voxel_size_z_);
  for(unsigned int k=0; k<image3D_.size(); k++) {
    fprintf(f,"%d ",(int) image3D_[k]);
  }

  fclose(f);
  return 0;

}

//==============================================================================

/**
 * \brief Function to write an Analyze image to disk.
 * \param[in] fileName: name of the image to write.
 * \param[in] type: the data type of the image.
 * \return Returns 0 if the writing operation ends successfully, -1 if the image
 *         size is 0, or -2 if the opening operation fails.
 * \author Luis Alvarez
*/
template <class  T>
int image3D<T>::write_analyze(const char *fileName, const int type) const
{
  if(image3D_.size()==0)
    return -1;

  int imaSize = size2d_;
  int size    = imaSize*nSlices_;

  char fitxer_proj[400], fitxer_img[400], fitxer_hdr[400];
  memset(fitxer_proj, 0, 400);
  memset(fitxer_img,  0, 400);
  memset(fitxer_hdr,  0, 400);

  strcpy(fitxer_proj, fileName);
  strncpy(fitxer_hdr, fitxer_proj, strlen(fitxer_proj)-3);
  strcat(fitxer_hdr,"hdr");

  strncpy(fitxer_img, fitxer_hdr, strlen(fitxer_hdr)-3);
  strcat(fitxer_img,"img");

  int nitems;
  struct dsr hdr;
  FILE *file;
  memset(&hdr, 0, sizeof(hdr));



  if((file=fopen(fitxer_img,"wb"))==NULL) {
    cout << "Error opening output file " << fitxer_img << endl;
    int j;
    int val = scanf("%d",&j);
    if(val==EOF)
      printf("Error writing value\n");
    exit(EXIT_FAILURE);
  }

  hdr.hk.sizeof_hdr = sizeof(hdr);
  for (int i=0; i<10; i++)
    hdr.hk.data_type[i] = 0;
  for (int i=0; i<18; i++)
    hdr.hk.db_name[i]  = 0;
  hdr.hk.extents       = 0;
  hdr.hk.session_error = 0;
  hdr.hk.regular       = 'r';
  hdr.hk.hkey_un0      = 0;

  for (int i=0; i<8; i++)
    hdr.dime.dim[i] = 0;
  for (int i=0; i<4; i++)
    hdr.dime.vox_units[i] = 0;
  for (int i=0; i<8; i++)
    hdr.dime.cal_units[i] = 0;
  hdr.dime.unused1  = 0;
  hdr.dime.datatype = 0;
  hdr.dime.bitpix   = 0;
  hdr.dime.dim_un0  = 0;
  for (int i=0; i<8; i++)
    hdr.dime.pixdim[i] = 0.0;
  hdr.dime.vox_offset = 0.0;
  hdr.dime.funused1   = 0.0;
  hdr.dime.funused2   = 0.0;
  hdr.dime.funused3   = 0.0;
  hdr.dime.cal_max    = 0.0;
  hdr.dime.cal_min    = 0.0;
  hdr.dime.compressed = 0;
  hdr.dime.verified   = 0;
  hdr.dime.glmax      = 0;
  hdr.dime.glmin      = 0;

  hdr.dime.dim[0]     = 4;
  hdr.dime.dim[1]     = width_;
  hdr.dime.dim[2]     = height_;
  hdr.dime.dim[3]     = nSlices_;
  hdr.dime.dim[4]     = 1;

  hdr.dime.pixdim[1]  = voxel_size_x_;
  hdr.dime.pixdim[2]  = voxel_size_y_;
  hdr.dime.pixdim[3]  = voxel_size_z_;

  hdr.dime.funused1   = 1.0; // scale factor

  for (int i=0; i<80; i++)
    hdr.hist.descrip[i] = 0;
  for (int i=0; i<24; i++)
    hdr.hist.aux_file[i]= 0;
  hdr.hist.orient = 0;
  for (int i=0; i<10; i++)
    hdr.hist.originator[i] = 0;
  for (int i=0; i<10; i++)
    hdr.hist.generated[i] = 0;
  for (int i=0; i<10; i++)
    hdr.hist.scannum[i] = 0;
  for (int i=0; i<10; i++)
    hdr.hist.patient_id[i] = 0;
  for (int i=0; i<10; i++)
    hdr.hist.exp_date[i] = 0;
  for (int i=0; i<10; i++)
    hdr.hist.exp_time[i] = 0;
  for (int i=0; i<3; i++)
    hdr.hist.hist_un0[i] = 0;
  hdr.hist.views       = 0;
  hdr.hist.vols_added  = 0;
  hdr.hist.start_field = 0;
  hdr.hist.field_skip  = 0;
  hdr.hist.omax        = 0;
  hdr.hist.omin        = 0;
  hdr.hist.smax        = 0;
  hdr.hist.smin        = 0;



  if(type==UCHAR) {
    hdr.dime.bitpix   = 8;
    hdr.dime.datatype = DT_UNSIGNED_CHAR;

    unsigned char *array1D= new unsigned char[size];
    for (int k=0; k<size; k++) {
      array1D[k] = (unsigned char) (image3D_[k]<0 ? 0:(image3D_[k]>255 ? 255:
                                    image3D_[k]));
    }
    nitems = size;
    fwrite(array1D, sizeof(unsigned char), nitems, file);
    delete [] array1D;
  }

  else if(type==USHORT) {
    hdr.dime.bitpix   = 16;
    hdr.dime.datatype = DT_UNSIGNED_SHORT;

    short *array1D = new short[size];
    for (int k=0; k<size; k++) {
      array1D[k] = (short) (image3D_[k]<0 ? 0:(image3D_[k]>65024 ?
                                               65024:image3D_[k]));
    }
    nitems = size;
    fwrite(array1D, sizeof(short), nitems, file);
    delete [] array1D;
  }

  else if(type==SHORT) {

    hdr.dime.bitpix   = 16;
    hdr.dime.datatype = DT_SIGNED_SHORT;

    short *array1D = new short[size];
    for (int k=0; k<size; k++) {
      array1D[k] = (short) (image3D_[k]<-32512 ? -32512:(image3D_[k]>32512 ?
                                                         32512:image3D_[k]));
    }
    nitems = size;
    fwrite(array1D, sizeof(short), nitems, file);
    delete [] array1D;
  }

  else if(type==FLOAT) {

    hdr.dime.bitpix   = 32;
    hdr.dime.datatype = DT_FLOAT;

    float *array1D = new float[size];
    for (int k=0; k<size; k++) {
      array1D[k] = (float) (image3D_[k]);
    }
    nitems = size;
    fwrite(array1D, sizeof(float), nitems, file);
    delete [] array1D;
  } else {
    cout << "Invalid data type" << endl;
    fclose(file);
    exit(EXIT_FAILURE);
  }

  fclose(file);

  FILE *fp;
  if ((fp = fopen(fitxer_hdr,"wb")) ==NULL)
    return -2;
  fwrite (&hdr,1,sizeof(struct dsr),fp);
  fclose(fp);

  return 0;
}

//==============================================================================

/**
 *
 * \fn image3D<T>::~image3D()
 * \brief Basic destructor
 * \author Luis Alvarez
*/
template <class  T>
image3D<T>::~image3D()
{
  width_ = height_ = nSlices_ = size2d_ = 0;
  image3D_.clear();
}

//==============================================================================

/**
 * \brief Basic constructor
 * \param[in] width: width of the 3D image
 * \param[in] height: height of the 3D image
 * \param[in] nSlices: number of slices in the 3D image
 * \param[in] voxel_size_x: voxel size in X
 * \param[in] voxel_size_y: voxel size in Y
 * \param[in] voxel_size_z: voxel size in Z
 * \param[in] a: value to initialize the whole image
 * \author Luis Alvarez
*/
template <class  T>
image3D<T>::image3D(const int width, const int height, const int nSlices,
                    const float voxel_size_x, const float voxel_size_y,
                    const float voxel_size_z, const T &a)
{
  width_   = width;
  size2d_  = width*height;
  width_   = width;
  height_  = height;
  nSlices_ = nSlices;

  voxel_size_x_ = voxel_size_x;
  voxel_size_y_ = voxel_size_y;
  voxel_size_z_ = voxel_size_z;

  if(nSlices_==0)
    return;

  int size = size2d_*nSlices_;

  image3D_.resize(size);
  for(int k=0; k<size; k++) {
    image3D_[k] = a;
  }

}

//==============================================================================

/**
 * \brief Basic constructor
 * \param[in] width: width of the 3D image
 * \param[in] height: height of the 3D image
 * \param[in] nSlices: number of slices in the 3D image
 * \param[in] voxel_size_x: voxel size in X
 * \param[in] voxel_size_y: voxel size in Y
 * \param[in] voxel_size_z: voxel size in Z
 * \author Luis Alvarez
*/
template <class  T>
image3D<T>::image3D(const int width, const int height, const int nSlices,
                    const float voxel_size_x, const float voxel_size_y,
                    const float voxel_size_z)
{
  width_   = width;
  size2d_  = width*height;
  width_   = width;
  height_  = height;
  nSlices_ = nSlices;

  voxel_size_x_ = voxel_size_x;
  voxel_size_y_ = voxel_size_y;
  voxel_size_z_ = voxel_size_z;

  if(nSlices_==0)
    return;

  int size=size2d_*nSlices_;

  image3D_.resize(size);

}

//==============================================================================

/**
 * \fn image3D<T>::image3D()
 * \brief Default constructor
 * \author Luis Alvarez
*/
template <class  T>
image3D<T>::image3D() {}

//==============================================================================

/**
 * \brief Assignment operator
 * \param[in] image2: the 3D image to assign
 * \author Luis Alvarez
*/
template <class  T> template <class  U>
image3D<T> & image3D<T>::operator=(const image3D<U> &image2)
{

  width_    = image2.width();
  size2d_   = image2.width()*image2.height();
  height_   = image2.height();
  nSlices_  = image2.nSlices();

  voxel_size_x_ = image2.get_voxel_size_x();
  voxel_size_y_ = image2.get_voxel_size_y();
  voxel_size_z_ = image2.get_voxel_size_z();

  if(nSlices_==0)
    return *this;

  int size = size2d_*nSlices_;

  if(image3D_.size()!=size)
    image3D_.resize(size);

  #ifdef _OPENMP
  #pragma omp parallel for shared(image2)
  #endif
  for(int k=0; k<size; k++) {
    #ifdef IMAGE3D_DEBUG
      printf("k=%d\n",k);
    #endif // IMAGE3D_DEBUG
    image3D_[k] = (T) image2[k];
  }

  return *this;
}

//==============================================================================

/** \brief Method to make a copy of a 3D image
 *
 * \param[in] I2: original 3D image to copy
 * \return Copies in the current object the 3D image provided
 *
 */

template <class  T>
void image3D<T>::copy(const image3D<T> &I2)
{
  width_   = I2.width();
  height_  = I2.height();
  nSlices_ = I2.nSlices();
  size2d_  = width_*height_;

  voxel_size_x_ = I2.voxel_size_x();
  voxel_size_y_ = I2.voxel_size_y();
  voxel_size_z_ = I2.voxel_size_z();

  if(image3D_.size()!=I2.size())
    image3D_.resize(I2.size());

  #ifdef _OPENMP
  #pragma omp parallel for shared(I2)
  #endif
  for(int k=0; k<I2.size(); k++) {
    #ifdef IMAGE3D_DEBUG
      printf("k=%d\n",(int) k); system("pause");
    #endif // IMAGE3D_DEBUG
    image3D_[k] = I2[k];
  }
}


#endif
