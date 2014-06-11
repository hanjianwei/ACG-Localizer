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

#ifndef LOCALIZE_ME_EXIF_READER_HH
#define LOCALIZE_ME_EXIF_READER_HH

/**
 *    Class to read information from the EXIF-tag of an JPEG image.
 *    This class is basically a wrapper for functionality implemented in
 *    the jhead tool by Matthias Wandel (http://www.sentex.ca/~mwandel/jhead/).
 *  
 *  author : Torsten Sattler (tsattler@cs.rwth-aachen.de)
 *  date : 07-26-2010
**/

#include "jhead-2.90/jhead.hh"

class exif_reader
{
  
  public:
	
	//! read the exif tag from a jpeg file
	static void open_exif( const char * filename );

	//! get the focal length from the exif tag previously open by open_exif
	static float get_focal_length( );

	//! get the width of the CCD from the exif tag, returns -1.0 if not available
	static float get_CCD_width();
	
	//! get the width of the image from the exif tag
	static int get_image_width();
	
	//! get the height of the image from the exif tag
	static int get_image_height();
	
        //! check if the exif data cotnains gps tags
	static bool is_gps_data_present();
	
        //! get the latitude from the gps tag
	static double get_gps_latitude();
        //! get the longitude from the gps tag
	static double get_gps_longitude();
        //! get the height from the gps tag(in meters)
	static double get_gps_height();
	
	/**
	 * get the orientation of the image: 
	 * 1 - normal orientation
	 * 6 - rotated by 90°
	 * 3 - rotated by 180°
	 * 8 - rotated by 270°
	 * rest is mirroring
	**/
	static int get_image_orientation();

	//! free the memory previously used for reading the exif tag
	static void close_exif( );
  
};

#endif
