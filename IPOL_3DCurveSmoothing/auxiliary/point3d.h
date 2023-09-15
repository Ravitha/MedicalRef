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
 * \file point3d.h
 * \brief point3d class
 * \author Luis Alvarez
*/
#ifndef POINT3D_H
#define POINT3D_H
#include <iostream>
#include <math.h>
#include <vector>


/**
 * \class  point3d
 * \brief Class  to store 3D points and basic methods
 * \author Luis Alvarez
 */
template < class  T >
class point3d
{
public :
  /// \brief point x coordinate
  T x;
  /// \brief point y coordinate
  T y;
  /// \brief point z coordinate
  T z;
//------------------------------------------------------------------------------
  ///Default constructor
  point3d():x( (T) 0), y( (T) 0), z( (T) 0) {};
//------------------------------------------------------------------------------
  ///Default destructor
  ~point3d() {};
//------------------------------------------------------------------------------
  ///Constructor with parameters
  point3d(T xx, T yy, T zz)
  {
    x = xx;
    y = yy;
    z = zz;
  }
//------------------------------------------------------------------------------
  ///Constructor with two points as parameters
  point3d(const point3d<T> &p,const point3d<T> &q)
  {
    x = p.y*q.z-p.z*q.y;
    y = p.z*q.x-p.x*q.z;
    z = p.x*q.y-p.y*q.x;
  }
//------------------------------------------------------------------------------
  ///Addition assignment operator
  void operator+=(const point3d &p)
  {
    x += p.x;
    y += p.y;
    z += p.z;
  }
//------------------------------------------------------------------------------
  ///Assignment operator
  point3d<T> operator=(const T &scalar) const
  {
    return point3d(scalar,scalar,scalar);
  }
//------------------------------------------------------------------------------
  ///Sum operator
  point3d<T> operator+(const point3d &p)const
  {
    return point3d(x+p.x,y+p.y,z+p.z);
  }
//------------------------------------------------------------------------------
  ///Subtraction operator
  point3d<T> operator-(const point3d &p)const
  {
    return point3d(x-p.x,y-p.y,z-p.z);
  }
//------------------------------------------------------------------------------
  ///Product operator (between 3D point and scalar value)
  point3d<T> operator*(const T &a) const
  {
    return point3d(a*x,a*y,a*z);
  }
//------------------------------------------------------------------------------
  ///Division operator (between 3D point and scalar value)
  point3d<T> operator/(const T &a) const
  {
    return point3d(x/a,y/a,z/a);
  }
//------------------------------------------------------------------------------
  ///Product operator (between 3D points)
  point3d<T>  operator*(const point3d<T> &p)const
  {
    return point3d(x*p.x,y*p.y,z*p.z);
  }
//------------------------------------------------------------------------------
  ///Division operator (between 3D points)
  point3d<T>  operator/(const point3d<T> &p)const
  {
    return point3d(x/p.x,y/p.y,z/p.z);
  }
//------------------------------------------------------------------------------
  ///Equality operator
  bool operator==(const point3d &p)
  {
    return (x == p.x && y == p.y && z == p.z);
  }
//------------------------------------------------------------------------------
  ///Not equal operator
  bool operator!=(const point3d &p)
  {
    return (x != p.x && y != p.y && z != p.z);
  }
//------------------------------------------------------------------------------
  ///3D point norm
  inline float norm() const
  {
    float paso = x*x+y*y+z*z;
    return(paso>0.?sqrtf((float) paso):0.);
  }
//------------------------------------------------------------------------------
  ///3D point squared
  inline float norm2() const
  {
    return(x*x+y*y+z*z);
  }
//------------------------------------------------------------------------------
  ///Squared product between two 3D points
  inline float norm2(const float dx,const float dy,const float dz) const
  {
    return(dx*dx*x*x+dy*dy*y*y+dz*dx*z*z);
  }
//------------------------------------------------------------------------------
  ///3D point norm (providing point coordinates)
  inline float norm(const float dx,const float dy,const float dz) const
  {
    float paso = norm2(dx,dy,dz);
    return(paso>0.?sqrtf(paso):0.);
  }
//------------------------------------------------------------------------------
  ///Scalar product
  inline double by(const point3d<T> &p3) const
  {
    return(x*p3.x+y*p3.y+z*p3.z);
  }
//------------------------------------------------------------------------------
  ///Print 3D point coordinates
  void print()
  {
    std::cout << "point3d : (" << x << "," << y << "," << z <<")" << std::endl;
  }
//------------------------------------------------------------------------------
  ///Minimum computation
  T min()
  {
    if(x<=y && x<=z)
      return x;
    if(y<=x && y<=z)
      return y;
    return z;
  }
//------------------------------------------------------------------------------
  ///3D point normalization
  point3d<T> normalize() const
  {
    float mod = x*x+y*y+z*z;
    if(mod>0) {
      mod=sqrtf(mod);
      return point3d<T>(x / mod, y / mod, z / mod);
    } else
      return point3d<T>();
  }

};


#endif
