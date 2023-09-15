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

#include "obj3D.h"

///Face vertices
static string face("f -4 -3 -2 -1\n");

/** \brief Private method used to update a new material in the content string
 *
 * \param[in] s: a formated string
 * \param[in] c: material color
 * \param[in] alpha: color transparency
 * \param[in] id: material id
 * \return Updates the content string
 *
 */
void obj3D::updateContentMaterial(string s, color c, float alpha, int id)
{
  char buf[1000];

  //We check if there are format specifiers
  if(s.find("%")!=string::npos) {
    if(id==-1 && alpha==-1) //Color
      sprintf(buf,s.c_str(), (float)c[0]/255., (float)c[1]/255.,
              (float)c[2]/255.);
    else { //Check if it is material id or transparency
      if(id!=-1)
        sprintf(buf, s.c_str(), id);
      if(alpha!=-1)
        sprintf(buf, s.c_str(), alpha);
    }
    content += buf;
  } else //Otherwise, we're adding a string without format specifiers
    content += s;
}

//==============================================================================

/** \brief Private method used to update the geometry in the content string
 *
 * \param[in] s: a formated string
 * \param[in] id: material id
 * \param[in] p: geometry coordinates
 * \param[in] dx: displacement in x
 * \param[in] dy: displacement in y
 * \param[in] dz: displacement in z
 * \return Updates the content string
 *
 */
void obj3D::updateContentGeometry(string s, int id, point3d<double> p,
                                  float dx, float dy, float dz)
{
  char buf[1000];
  int maxz = nSlices-1;
  //We check if there are format specifiers
  if(s.find("%")!=string::npos) {
    if(id!=-1) //material
      sprintf(buf, s.c_str(), id);
    else { //coordinates
      sprintf(buf, s.c_str(), (p.x+dx)*vsx, (p.y+dy)*vsy, (maxz-(p.z+dz))*vsz);
    }
    content += buf;
  } else //Otherwise, we're adding a string without format specifiers
    content += s;
}

//==============================================================================

/** \brief Private method used to update the geometry in the content string in
 *         mm
 *
 * \param[in] s: a formated string
 * \param[in] id: material id
 * \param[in] p: geometry coordinates in mm
 * \param[in] dx: displacement in x
 * \param[in] dy: displacement in y
 * \param[in] dz: displacement in z
 * \return Updates the content string
 *
 */
void obj3D::updateContentGeometryMM(string s, int id, point3d<double> p,
                                    float dx, float dy, float dz)
{
  char buf[1000];
  int maxz = nSlices-1;
  //We check if there are format specifiers
  if(s.find("%")!=string::npos) {
    if(id!=-1) //material
      sprintf(buf, s.c_str(), id);
    else { //coordinates
      sprintf(buf, s.c_str(), (p.x+dx), (p.y+dy), (maxz-(p.z+dz)));
    }
    content += buf;
  } else //Otherwise, we're adding a string without format specifiers
    content += s;
}

//==============================================================================

/** \brief Private method to add a new material to the obj file (called by other
 *         add methods)
 *
 * \param[in] c: the color of the new material
 * \param[in] alpha: the transparency of the new material
 * \return Writes in the content string the new material.
 *
 */
void obj3D::addMaterial(color c, float alpha)
{
  id++;
  /// Add material
  updateContentMaterial(
               "# To visualize this file, drag to url http://3dviewer.net\n\n");
  updateContentMaterial("# Material:\n");
  updateContentMaterial("newmtl material%d\n",color(0,0,0),-1,id);
  updateContentMaterial("Ka 0.2 0.2 0.2\n");
  /// material color and transparency
  updateContentMaterial("Kd %f %f %f\n",c,-1.,-1);
  updateContentMaterial("Ks 0.2 0.2 0.2\n");
  updateContentMaterial("Ns 20  # shininess\n");
  updateContentMaterial("d %.2f # transparency\n", color(0,0,0), alpha, -1);
}

//==============================================================================

/** \brief Constructor without parameters
 * \return Returns an empty obj3D object
 *
 */
obj3D::obj3D()
{
  content = "";
  id      = width = height = nSlices = size2d = 0;
  vsx     = vsy   = vsz    = 0.;
  vs      = point3d<double>(0.,0.,0.);
}

//==============================================================================

/** \brief Constructor with parameters
 *
 * \param[in] i: an image3D used to get the image parameters
 * \return Returns an obj3D object
 *
 */
obj3D::obj3D(image3D<unsigned char> &i)
{
  content = "";
  width   = i.width();
  height  = i.height();
  nSlices = i.nSlices();
  size2d  = i.size2d();
  vsx     = i.voxel_size_x();
  vsy     = i.voxel_size_y();
  vsz     = i.voxel_size_z();
  vs      = i.voxel_size();
  id      = 0;
}

//==============================================================================

/** \brief Constructor using an object of class obj3D
 *
 * \param[in] obj: an obj3D object
 * \return Returns a new obj3D object created from the parameter
 *
 */
/// Does not copy the content, just create the object. Use '=' instead.
obj3D::obj3D(obj3D &obj)
{
  content = "";
  width   = obj.width;
  height  = obj.height;
  nSlices = obj.nSlices;
  vsx     = obj.vsx;
  vsy     = obj.vsy;
  vsz     = obj.vsz;
  vs      = obj.vs;
  size2d  = obj.size2d;
  id      = 0;

}

//==============================================================================

/** \brief Method to add a segmentation considering the values of the neighbors
 *         in the three main directions of a "true" voxel
 *
 * \param[in] S: reference to an image3D containing the segmentation
 * \param[in] c: color for adding the segmentation
 * \param[in] alpha: volume transparency
 * \return Adds to the content string the material and geometry for the
 *         segmentation
 *
 */
void obj3D::addSegmentation(image3D<unsigned char> &S, color c, float alpha)
{
  addMaterial(c,alpha);

  updateContentGeometry("\n# Geometry:\n");
  updateContentGeometry("g segmentation%d\n",id);
  updateContentGeometry("usemtl material%d\n",id);
  for(int z=1; z<nSlices-1; z++)
    for(int y=1; y<height-1; y++)
      for (int x=1; x<width-1; x++)
        if (S(x,y,z)>0) {
          point3d<double> p(x,y,z);
          if (S(x,y,z-1) == 0) {
            updateContentGeometry("v %f %f %f\n", -1, p, 0.5, 0.5, -0.5);
            updateContentGeometry("v %f %f %f\n", -1, p, 0.5, -0.5, -0.5);
            updateContentGeometry("v %f %f %f\n", -1, p, -0.5, 0.5, -0.5);
            updateContentGeometry("v %f %f %f\n", -1, p, -0.5, -0.5, -0.5);
            updateContentGeometry(face.c_str());
          }
          if (S(x,y,z+1) == 0) {
            updateContentGeometry("v %f %f %f\n", -1, p, -0.5, -0.5, 0.5);
            updateContentGeometry("v %f %f %f\n", -1, p, -0.5, 0.5, 0.5);
            updateContentGeometry("v %f %f %f\n", -1, p, 0.5, 0.5, 0.5);
            updateContentGeometry("v %f %f %f\n", -1, p, 0.5, -0.5, 0.5);
            updateContentGeometry(face.c_str());
          }
          if (S(x-1,y,z) == 0) {
            updateContentGeometry("v %f %f %f\n", -1, p, -0.5, -0.5, 0.5);
            updateContentGeometry("v %f %f %f\n", -1, p, -0.5, -0.5, -0.5);
            updateContentGeometry("v %f %f %f\n", -1, p, -0.5, 0.5, -0.5);
            updateContentGeometry("v %f %f %f\n", -1, p, -0.5, 0.5, 0.5);
            updateContentGeometry(face.c_str());
          }
          if (S(x+1,y,z) == 0) {
            updateContentGeometry("v %f %f %f\n", -1, p, 0.5, 0.5, 0.5);
            updateContentGeometry("v %f %f %f\n", -1, p, 0.5, 0.5, -0.5);
            updateContentGeometry("v %f %f %f\n", -1, p, 0.5, -0.5, -0.5);
            updateContentGeometry("v %f %f %f\n", -1, p, 0.5, -0.5, 0.5);
            updateContentGeometry(face.c_str());
          }
          if (S(x,y-1,z) == 0) {
            updateContentGeometry("v %f %f %f\n", -1, p, 0.5, -0.5, 0.5);
            updateContentGeometry("v %f %f %f\n", -1, p, 0.5, -0.5, -0.5);
            updateContentGeometry("v %f %f %f\n", -1, p, -0.5, -0.5, -0.5);
            updateContentGeometry("v %f %f %f\n", -1, p, -0.5, -0.5, 0.5);
            updateContentGeometry(face.c_str());
          }
          if (S(x,y+1,z) == 0) {
            updateContentGeometry("v %f %f %f\n", -1, p, -0.5, 0.5, 0.5);
            updateContentGeometry("v %f %f %f\n", -1, p, -0.5, 0.5, -0.5);
            updateContentGeometry("v %f %f %f\n", -1, p, 0.5, 0.5, -0.5);
            updateContentGeometry("v %f %f %f\n", -1, p, 0.5, 0.5, 0.5);
            updateContentGeometry(face.c_str());
          }
        }

}

//==============================================================================

/** \brief Method to add a line given by a vector of contiguous points
 *
 * \param[in] cl: reference to a vector containing the 3D points of the line
 * \param[in] c: the color for the curve3D
 * \param[in] alpha: curve3D transparency
 * \return Adds to the content string the material and geometry for the line
 *
 */
void obj3D::addLine(vector<point3d<double> > &cl, color c, float alpha)
{
  addMaterial(c,alpha);

  updateContentGeometry("\n# Geometry:\n");
  updateContentGeometry("g segmentation%d\n",id);
  updateContentGeometry("usemtl material%d\n",id);
  for(unsigned int i=0; i<cl.size(); i++) {
    point3d<double> p(cl[i].x,cl[i].y,cl[i].z);
    //S(x,y,z-1)
    updateContentGeometry("v %f %f %f\n", -1, p, 0.5, 0.5, -0.5);
    updateContentGeometry("v %f %f %f\n", -1, p, 0.5, -0.5, -0.5);
    updateContentGeometry("v %f %f %f\n", -1, p, -0.5, 0.5, -0.5);
    updateContentGeometry("v %f %f %f\n", -1, p, -0.5, -0.5, -0.5);
    updateContentGeometry(face.c_str());
    //S(x,y,z+1)
    updateContentGeometry("v %f %f %f\n", -1, p, -0.5, -0.5, 0.5);
    updateContentGeometry("v %f %f %f\n", -1, p, -0.5, 0.5, 0.5);
    updateContentGeometry("v %f %f %f\n", -1, p, 0.5, 0.5, 0.5);
    updateContentGeometry("v %f %f %f\n", -1, p, 0.5, -0.5, 0.5);
    updateContentGeometry(face.c_str());
    //S(x-1,y,z)
    updateContentGeometry("v %f %f %f\n", -1, p, -0.5, -0.5, 0.5);
    updateContentGeometry("v %f %f %f\n", -1, p, -0.5, -0.5, -0.5);
    updateContentGeometry("v %f %f %f\n", -1, p, -0.5, 0.5, -0.5);
    updateContentGeometry("v %f %f %f\n", -1, p, -0.5, 0.5, 0.5);
    updateContentGeometry(face.c_str());
    //S(x+1,y,z)
    updateContentGeometry("v %f %f %f\n", -1, p, 0.5, 0.5, 0.5);
    updateContentGeometry("v %f %f %f\n", -1, p, 0.5, 0.5, -0.5);
    updateContentGeometry("v %f %f %f\n", -1, p, 0.5, -0.5, -0.5);
    updateContentGeometry("v %f %f %f\n", -1, p, 0.5, -0.5, 0.5);
    updateContentGeometry(face.c_str());
    //S(x,y-1,z)
    updateContentGeometry("v %f %f %f\n", -1, p, 0.5, -0.5, 0.5);
    updateContentGeometry("v %f %f %f\n", -1, p, 0.5, -0.5, -0.5);
    updateContentGeometry("v %f %f %f\n", -1, p, -0.5, -0.5, -0.5);
    updateContentGeometry("v %f %f %f\n", -1, p, -0.5, -0.5, 0.5);
    updateContentGeometry(face.c_str());
    //S(x,y+1,z)
    updateContentGeometry("v %f %f %f\n", -1, p, -0.5, 0.5, 0.5);
    updateContentGeometry("v %f %f %f\n", -1, p, -0.5, 0.5, -0.5);
    updateContentGeometry("v %f %f %f\n", -1, p, 0.5, 0.5, -0.5);
    updateContentGeometry("v %f %f %f\n", -1, p, 0.5, 0.5, 0.5);
    updateContentGeometry(face.c_str());
  }
}

//==============================================================================

/** \brief Method to add a line given by a vector of contiguous points in mm
 *
 * \param[in] cl: reference to a vector containing the 3D points of the line in
 *                mm
 * \param[in] c: the color for the curve3D
 * \param[in] alpha: curve3D transparency
 * \return Adds to the content string the material and geometry for the line
 *         without using the voxel size
 *
 */
void obj3D::addLineMM(vector<point3d<double> > &cl, color c, float alpha)
{
  addMaterial(c,alpha);

  updateContentGeometryMM("\n# Geometry:\n");
  updateContentGeometryMM("g segmentation%d\n",id);
  updateContentGeometryMM("usemtl material%d\n",id);
  double scale = 1.;
  if(cl.size()>1.) {
    scale = (cl[1]-cl[0]).norm();
  }
  for(unsigned int i=0; i<cl.size(); i++) {
    point3d<double> p(cl[i].x,cl[i].y,cl[i].z);
    //S(x,y,z-1)
    updateContentGeometryMM("v %f %f %f\n", -1, p, 0.5*scale, 0.5*scale,
                            -0.5*scale);
    updateContentGeometryMM("v %f %f %f\n", -1, p, 0.5*scale, -0.5*scale,
                            -0.5*scale);
    updateContentGeometryMM("v %f %f %f\n", -1, p, -0.5*scale, 0.5*scale,
                            -0.5*scale);
    updateContentGeometryMM("v %f %f %f\n", -1, p, -0.5*scale, -0.5*scale,
                            -0.5*scale);
    updateContentGeometryMM(face.c_str());
    //S(x,y,z+1)
    updateContentGeometryMM("v %f %f %f\n", -1, p, -0.5*scale, -0.5*scale,
                            0.5*scale);
    updateContentGeometryMM("v %f %f %f\n", -1, p, -0.5*scale, 0.5*scale,
                            0.5*scale);
    updateContentGeometryMM("v %f %f %f\n", -1, p, 0.5*scale, 0.5*scale,
                            0.5*scale);
    updateContentGeometryMM("v %f %f %f\n", -1, p, 0.5*scale, -0.5*scale,
                            0.5*scale);
    updateContentGeometryMM(face.c_str());
    //S(x-1,y,z)
    updateContentGeometryMM("v %f %f %f\n", -1, p, -0.5*scale, -0.5*scale,
                            0.5*scale);
    updateContentGeometryMM("v %f %f %f\n", -1, p, -0.5*scale, -0.5*scale,
                            -0.5*scale);
    updateContentGeometryMM("v %f %f %f\n", -1, p, -0.5*scale, 0.5*scale,
                            -0.5*scale);
    updateContentGeometryMM("v %f %f %f\n", -1, p, -0.5*scale, 0.5*scale,
                            0.5*scale);
    updateContentGeometryMM(face.c_str());
    //S(x+1,y,z)
    updateContentGeometryMM("v %f %f %f\n", -1, p, 0.5*scale, 0.5*scale,
                            0.5*scale);
    updateContentGeometryMM("v %f %f %f\n", -1, p, 0.5*scale, 0.5*scale,
                            -0.5*scale);
    updateContentGeometryMM("v %f %f %f\n", -1, p, 0.5*scale, -0.5*scale,
                            -0.5*scale);
    updateContentGeometryMM("v %f %f %f\n", -1, p, 0.5*scale, -0.5*scale,
                            0.5*scale);
    updateContentGeometryMM(face.c_str());
    //S(x,y-1,z)
    updateContentGeometryMM("v %f %f %f\n", -1, p, 0.5*scale, -0.5*scale,
                            0.5*scale);
    updateContentGeometryMM("v %f %f %f\n", -1, p, 0.5*scale, -0.5*scale,
                            -0.5*scale);
    updateContentGeometryMM("v %f %f %f\n", -1, p, -0.5*scale, -0.5*scale,
                            -0.5*scale);
    updateContentGeometryMM("v %f %f %f\n", -1, p, -0.5*scale, -0.5*scale,
                            0.5*scale);
    updateContentGeometryMM(face.c_str());
    //S(x,y+1,z)
    updateContentGeometryMM("v %f %f %f\n", -1, p, -0.5*scale, 0.5*scale,
                            0.5*scale);
    updateContentGeometryMM("v %f %f %f\n", -1, p, -0.5*scale, 0.5*scale,
                            -0.5*scale);
    updateContentGeometryMM("v %f %f %f\n", -1, p, 0.5*scale, 0.5*scale,
                            -0.5*scale);
    updateContentGeometryMM("v %f %f %f\n", -1, p, 0.5*scale, 0.5*scale,
                            0.5*scale);
    updateContentGeometryMM(face.c_str());
  }
}

//==============================================================================

/** \brief Assignment operator
 *
 * \param[in] src: A reference to an obj3D object
 * \return A copy of the obj3D
 *
 */
obj3D obj3D::operator=(obj3D &src)
{
  content = src.content;
  width   = src.width;
  height  = src.height;
  nSlices = src.nSlices;
  vsx     = src.vsx;
  vsy     = src.vsy;
  vsz     = src.vsz;
  vs      = src.vs;
  size2d  = src.size2d;
  id      = src.id;
  return *this;
}

//==============================================================================

/** \brief Method to save the obj content in a file
 *
 * \param[in] filename: the path and name for the obj file
 * \return Writes in disc the obj3D object content
 *
 */
void obj3D::saveOBJ(string filename)
{
  std::ofstream f(filename.c_str(),std::ios::binary);
  f << content;
  f.close();
}

//==============================================================================

/** \brief Method to reset an obj3D object
 *
 * \return A clear obj3D object
 *
 */
void obj3D::resetOBJ()
{
  content = "";
  id      = width = height = nSlices = size2d = 0;
  vsx     = vsy   = vsz    = 0.;
  vs      = point3d<double>(0.,0.,0.);
}
