/*===========================================================================*\
 *                                                                           *
 *                            ACG Localizer                                  *
 *      Copyright (C) 2011 by Computer Graphics Group, RWTH Aachen           *
 *                           www.rwth-graphics.de                            *
 *                                                                           *
 *---------------------------------------------------------------------------* 
 *  This file is part of ACG Localizer                                       *
 *                                                                           *
 *  ACG Localizer is free software: you can redistribute it and/or modify    *
 *  it under the terms of the GNU General Public License as published by     *
 *  the Free Software Foundation, either version 3 of the License, or        *
 *  (at your option) any later version.                                      *
 *                                                                           *
 *  ACG Localizer is distributed in the hope that it will be useful,         *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *  GNU General Public License for more details.                             *
 *                                                                           *
 *  You should have received a copy of the GNU General Public License        *
 *  along with ACG Localizer.  If not, see <http://www.gnu.org/licenses/>.   *
 *                                                                           *
\*===========================================================================*/ 

#ifndef SIFT_KEYPOINT_HH
#define SIFT_KEYPOINT_HH

/**
 *    Simple class defining a SIFT-keypoint, excluding its descriptor.
 *
 *  author : Torsten Sattler (tsattler@cs.rwth-aachen.de)
 *  date : 11-02-2011
**/ 

class SIFT_keypoint
{
  public:
	float x;
	float y;
	float scale;
	float orientation;
	
	SIFT_keypoint( )
	{
	  x = y = scale = orientation = 0.0;
	}
	
	SIFT_keypoint( float x_, float y_, float scale_, float orientation_ )
	{
	  x = x_;
	  y = y_;
	  scale = scale_;
	  orientation = orientation_;
	}
};

enum SIFT_FORMAT{
  LOWE = 0,
  UNDEFINED = 3
};

#endif