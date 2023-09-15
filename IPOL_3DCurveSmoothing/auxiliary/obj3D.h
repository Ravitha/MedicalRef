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
 * \file obj3D.h
 * \brief Class to generate Wavefront .obj files
 * \author Daniel Santana-Cedrés
*/

#ifndef OBJ3D_H_INCLUDED
#define OBJ3D_H_INCLUDED

#include "image3D.h"
#include "color.h"

/** \brief Class to create and manage Wavefront .obj files by means of
 *         operations to draw volumes (given by segmentations in image3D
 *         objects) or lines.
 */
class obj3D
{
  // Dimensions of the image3D used for building the object
  /// \brief obj width
  int width,
  /// \brief obj height
  height,
  /// \brief obj number of slices
  nSlices,
  /// \brief obj size 2D
  size2d;
  // Voxel size attributes
  /// \brief voxel size in x
  double vsx,
  /// \brief voxel size in y
  vsy,
  /// \brief voxel size in z
  vsz;
  /// \brief point 3D containing the voxel size
  point3d<double> vs;
  /// \brief OBJ content
  string content;
  /// \brief Material id
  int id;
//==============================================================================
private:
  // Private method to add materials
  void addMaterial(color c, float alpha);
  // Private method to update material
  void updateContentMaterial(string s, color c=color(0,0,0), float alpha=-1,
                             int id=-1);
  // Private method to update geometry
  void updateContentGeometry(string s, int id=-1,
                             point3d<double> p=point3d<double>(-1.,-1.,-1.),
                             float dx=-1, float dy=-1, float dz=-1);
  // Private method to update geometry without using the voxel size
  void updateContentGeometryMM(string s, int id=-1,
                               point3d<double> p=point3d<double>(-1.,-1.,-1.),
                               float dx=-1, float dy=-1, float dz=-1);
//==============================================================================
public:
  // Constructor without parameters
  obj3D();
  // Constructor with parameters
  obj3D(image3D<unsigned char> &i);
  // Constructor using another obj3D object
  obj3D(obj3D &obj);
//------------------------------------------------------------------------------
  // Set/Get methods
  /// \brief Set obj width
  void setWidth(int w)
  {
    width=w;
  }
  /// \brief Get obj width
  int  getWidth()
  {
    return width;
  }
  /// \brief Set obj height
  void setHeight(int h)
  {
    height=h;
  }
  /// \brief Get obj height
  int  getHeight()
  {
    return height;
  }
  /// \brief Set obj number of slices
  void setnSlices(int nS)
  {
    nSlices=nS;
  }
  /// \brief Get obj number of slices
  int  getnSlices()
  {
    return nSlices;
  }
  /// \brief Set obj size 2D
  void setSize2D(int s2d)
  {
    size2d=s2d;
  }
  /// \brief Get obj size 2D
  int  getSize2D()
  {
    return size2d;
  }
//------------------------------------------------------------------------------
  // Voxel size
  /// \brief Set voxel size in x
  void  setVSX(float vs)
  {
    vsx=vs;
  }
  /// \brief Get voxel size in x
  float getVSX()
  {
    return vsx;
  }
  /// \brief Set voxel size in y
  void  setVSY(float vs)
  {
    vsy=vs;
  }
  /// \brief Get voxel size in y
  float getVSY()
  {
    return vsy;
  }
  /// \brief Set voxel size in z
  void  setVSZ(float vs)
  {
    vsz=vs;
  }
  /// \brief Get voxel size in z
  float getVSZ()
  {
    return vsz;
  }
  /// \brief Set voxel size
  void  setVS(point3d<double> p)
  {
    vs=p;
  }
  /// \brief Get voxel size
  point3d<double> getVS()
  {
    return vs;
  }
//------------------------------------------------------------------------------
  //CONTENT
  /// \brief Set obj content
  void setContent(string c)
  {
    content=c;
  }
  /// \brief Get obj content
  string getContent()
  {
    return content;
  }
//------------------------------------------------------------------------------
  /* Method to add a segmentation considering the values of the neighbors in
     the three main directions of a "true" voxel */
  void addSegmentation(image3D<unsigned char> &S, color c, float alpha);
//------------------------------------------------------------------------------
  // Method to add a line given by a vector of contiguous points
  void addLine(vector<point3d<double> > &cl, color c, float alpha);
//------------------------------------------------------------------------------
  // Method to add a line given by a vector of contiguous points in mm
  void addLineMM(vector<point3d<double> > &cl, color c, float alpha);
//------------------------------------------------------------------------------
  // Assignment operator
  obj3D operator=(obj3D &src);
//------------------------------------------------------------------------------
  // Method to save the obj3D content to a file
  void saveOBJ(string filename);
//------------------------------------------------------------------------------
  //Method to reset an obj3D object
  void resetOBJ();

};


#endif // OBJ3D_H_INCLUDED
