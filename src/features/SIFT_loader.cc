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

#include <iostream>
#include <fstream>
#include <cmath>
#include "SIFT_loader.hh"


SIFT_loader::SIFT_loader( )
{
  mNbFeatures = 0;
  mKeypoints.clear();
  mDescriptors.clear();
}

SIFT_loader::~SIFT_loader( )
{
  clear_data( );
}
    
void SIFT_loader::load_features( const char *filename, SIFT_FORMAT format )
{
  clear_data( );
  bool loaded = load_Lowe_features( filename );
  
  if( !loaded )
    std::cerr << "Could not load the features from the given file " << filename << std::endl;
}

bool SIFT_loader::save_features_lowe( const char *filename )
{
  //now open the .sift file and read out the interest points
  // format:
  //
  // nb_keypoints size_descriptor (should be 128)
  // for each keypoint:
  //	y x scale orientation(radians) 
  //    descripor( size_descriptor many unsigned char values) (6 lines with 20 values, 1 with 8 values)
  //the coordinates of the interest points are stored in a coordinate system in which the origin corresponds to the upper left of the image
  
  std::ofstream outstream( filename, std::ios::out );
  
  if ( !outstream.is_open() )
    return false;
  
  // save the number of keypoints and the size of the descriptors
  outstream << mNbFeatures << " 128" << std::endl;
  
  // load the keypoints and their descriptors
  for( uint32_t i=0; i<mNbFeatures; ++i )
  {
    outstream << mKeypoints[i].y << " " << mKeypoints[i].x << " " << mKeypoints[i].scale << " " << mKeypoints[i].orientation << std::endl;
    
    outstream << (int) mDescriptors[i][0] << " ";
    for( int j=1; j<128; ++j )
    {
      if( j%20 == 0 )
		outstream << std::endl;
	  outstream << (int) mDescriptors[i][j] << " ";
    }
    outstream << std::endl;
  }
  
  outstream.close();
  
  return true;
}


void SIFT_loader::clear_data( )
{
  for( uint32_t i=0; i<mNbFeatures; ++i )
  {
    if( mDescriptors[i] != 0 )
      delete [] mDescriptors[i];
    mDescriptors[i] = 0;
  }
  mKeypoints.clear();
  mDescriptors.clear();
}
    
uint32_t SIFT_loader::get_nb_features( )
{
  return mNbFeatures;
}

std::vector< unsigned char* >& SIFT_loader::get_descriptors( )
{
  return mDescriptors;
}

std::vector< SIFT_keypoint >& SIFT_loader::get_keypoints( )
{
  return mKeypoints;
}
    

bool SIFT_loader::load_Lowe_features( const char *filename )
{
  //now open the .sift file and read out the interest points
  // format:
  //
  // nb_keypoints size_descriptor (should be 128)
  // for each keypoint:
  //	y x scale orientation(radians) 
  //    descripor( size_descriptor many char values) (6 lines with 20 values, 1 with 8 values)
  //the coordinates of the interest points are stored in a coordinate system in which the origin corresponds to the upper left of the image
  
  std::ifstream instream( filename, std::ios::in );
  
  if ( !instream.is_open() )
    return false;
  
  // read the number of keypoints and the size of the descriptors
  uint32_t size_descriptor;
  instream >> mNbFeatures >> size_descriptor;
  
  if( size_descriptor != 128 )
    return false;
  
  // load the keypoints and their descriptors
  mKeypoints.resize(mNbFeatures);
  mDescriptors.resize(mNbFeatures,0);
  
  double x,y,scale,orientation;
  unsigned int descriptor_element;
  
  for( uint32_t i=0; i<mNbFeatures; ++i )
  {
    mDescriptors[i] = new unsigned char[128];
    instream >> y >> x >> scale >> orientation;
    mKeypoints[i] = SIFT_keypoint( x, y, scale, orientation );
    
    // read the descriptor
    for( int j=0; j<128; ++j )
    {
      instream >> descriptor_element;
      mDescriptors[i][j] = (unsigned char) descriptor_element;
    }
  }
  
  instream.close();
  
  return true;
}


