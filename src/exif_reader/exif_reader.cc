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

#include "exif_reader.hh"

#include <string>
#include <sstream>
#include <stdio.h>
#include <iostream>

void exif_reader::open_exif( const char * filename )
{
  // Code directly take from jhead.c
  // Many thanks to Matthias Wandel for writing this 
  // excellent piece of software
  ResetJpgfile();
  
  // Start with an empty image information structure.
  memset(&ImageInfo, 0, sizeof(ImageInfo));
  ImageInfo.FlashUsed = -1;
  ImageInfo.MeteringMode = -1;
  ImageInfo.Whitebalance = -1;
  
  ImageInfo.CCDWidth = -1.0;

  
  strncpy(ImageInfo.FileName, filename, PATH_MAX);

  ReadJpegFile(filename, READ_METADATA);
}
  
float exif_reader::get_focal_length( )
{
  return ImageInfo.FocalLength;
}
  
float exif_reader::get_CCD_width()
{
  return ImageInfo.CCDWidth;
}

int exif_reader::get_image_width()
{
  return ImageInfo.Width;
}

int exif_reader::get_image_height()
{
  return ImageInfo.Height;
}

double parse_double(const std::string & s)
{
  std::istringstream i(s);
   double x;
   if (!(i >> x))
     return 0;
   return x;

}

double convert_date_to_degree(const char data[31])
{
  std::string str(data, 31);
  
  size_t d_pos,m_pos,s_pos;
  
  double multiplier = 0;
  
  if(str[0] == 'S' || str[0] == 'W')
    multiplier = -1;
  else if(str[0] == 'E' || str[0] == 'N')
    multiplier = 1;
  else
    return 0;
  
  d_pos = str.find('d');
  m_pos = str.find('m');
  s_pos = str.find('s');
  
  if(d_pos == std::string::npos || m_pos == std::string::npos || s_pos == std::string::npos)
    return 0;
  
  std::string d_string, m_string, s_string;
  
  d_string = str.substr(2, d_pos-2);
  m_string = str.substr(d_pos+2, m_pos-d_pos-2);
  s_string = str.substr(m_pos+2, s_pos-m_pos-2);
  
  double day,min,sec;
  
  day = parse_double(d_string);
  min = parse_double(m_string);
  sec = parse_double(s_string);
  
  double res = sec/3600.0f+min/60.0f+day;

  res *= multiplier;
  
  return res;
}

bool exif_reader::is_gps_data_present()
{
  return ImageInfo.GpsInfoPresent == 1;
}
		
double exif_reader::get_gps_latitude()
{
  return convert_date_to_degree(ImageInfo.GpsLat);
}

double exif_reader::get_gps_longitude()
{
  return convert_date_to_degree(ImageInfo.GpsLong);
}

double exif_reader::get_gps_height()
{
  std::string str(ImageInfo.GpsAlt, 20);
  size_t m_pos = str.find('m');
  
  std::string height = str.substr(0,m_pos);
  
  return parse_double(height);
}

int exif_reader::get_image_orientation()
{
  return ImageInfo.Orientation;
}
  
void exif_reader::close_exif( )
{
  DiscardData();
}
