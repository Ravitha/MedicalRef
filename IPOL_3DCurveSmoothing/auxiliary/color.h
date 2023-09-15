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
 * \file color.h
 * \brief Class to manage color values
 * \author Luis Alvarez
*/


#ifndef COLOR_H_
#define COLOR_H_
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <math.h>

/**
 * \class color
 * \brief Class to store and manage RGB color values
 * \author Luis Alvarez
 */
class color
{
public :
  /** \brief red color channel */
  unsigned char R_;
  /** \brief green color channel */
  unsigned char G_;
  /** \brief blue color channel */
  unsigned char B_;
//------------------------------------------------------------------------------
  /** \brief color constructor with parameters
   *
   * \param[in] R: red color channel
   * \param[in] G: green color channel
   * \param[in] B: blue color channel
   * \return Builds and fills a color object with the parameters provided
   *
   */
  color(unsigned char R,unsigned char G,unsigned char B)
  {
    R_ = R;
    G_ = G;
    B_ = B;
  }
//------------------------------------------------------------------------------
  /** \brief Default constructor. */
  color()
  {
    R_ = (unsigned char) (255.*rand()/RAND_MAX);
    G_ = (unsigned char) (255.*rand()/RAND_MAX);
    B_ = (unsigned char) (255.*rand()/RAND_MAX);
  };
//------------------------------------------------------------------------------
  /** \brief Accessing color channels using operator []
   *
   * \param[in] k: index of the channel to access
   * \return Returns the channel value (unsigned char)
   *
   */
  inline unsigned char &operator[](int k)
  {
    switch (k) {
    case 0 :
      return(R_);
    case 1 :
      return(G_);
    default :
      return(B_);
    }
  }
};


#endif
