/*===========================================================================*\
 *                                                                           *
 *                            ACG Localizer                                  *
 *      Copyright (C) 2011-2012 by Computer Graphics Group, RWTH Aachen      *
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

#ifndef PARSE_BUNDLER_HH
#define PARSE_BUNDLER_HH

/**
 * Parses a bundle.out file created by Bundler into a set of 
 * 3D points and the corresponding descriptors.
 *
 * author: Torsten Sattler (tsattler@cs.rwth-aachen.de)
 * date  : 11-08-2012
**/

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdint.h>
#include <cstdlib>
#include <string>
#include "../features/SIFT_keypoint.hh"
#include "../features/SIFT_loader.hh"
#include "bundler_camera.hh"

////////////////////////////
// Simple definition of a 3D point
// consisting of x,y,z coordinates
// and a point color.
////////////////////////////
class point3D
{
  public:
    float x,y,z;
    
    unsigned char r,g,b;
    
    point3D()
    {
      x = y = z = 0.0;
      r = g = b = 0;
    }
    
    point3D( float x_, float y_, float z_ )
    {
      x = x_;
      y = y_;
      z = z_;
      r = g = b = 0;
    }
    
    point3D( float x_, float y_, float z_, unsigned char r_, unsigned char g_, unsigned char b_ )
    {
      x = x_;
      y = y_;
      z = z_;
      r = r_;
      g = g_;
      b = b_;
    }
    
    point3D( const point3D &other )
    {
      x = other.x;
      y = other.y;
      z = other.z;
      r = other.r;
      g = other.g;
      b = other.b;
    }
    
    
    void operator=( const point3D &other )
    {
      x = other.x;
      y = other.y;
      z = other.z;
      r = other.r;
      g = other.g;
      b = other.b;
    }
};

////////////////////////////
// A view of a 3D point as defined by Bundler. It
// contains the index of the camera the point was seen it,
// the keypoint index in the .key file of that camera corresponding
// to the projection of the 3D point, the x- and y-positions of that keypoint 
// in a reference frame where the center of the coordinate system is the center
// of the cameras image, as well as the scale and orientation of that keypoint.
// (see http://phototour.cs.washington.edu/bundler/bundler-v0.4-manual.html)
////////////////////////////
struct view
{
  uint32_t camera; // index of the camera this feature was detected in
  uint32_t key; // index of this feature in the corresponding .key file containing the cameras SIFT features
  float x,y;
  float scale; // scale at which the feature was detected as taken from the original image
  float orientation; // orientation of that features as detected in the image
    
  view()
  {
    camera = key = 0;
    x = y = 0.0;
    scale = 0.0f;
    orientation = 0.0f;
  }
};


////////////////////////////
// Information about a reconstructed 3D point. The structure contains information about the 3D position
// of the point, its color, and a list of views in which the point can be seen and the corresponding SIFT descriptors.
// The descriptor for view view_list[i] is stored in descriptors[128*i]...descriptors[128*i+127]
////////////////////////////
class feature_3D_info
{
  public:
    point3D point;
    std::vector< view > view_list; 

    std::vector< unsigned char > descriptors;
    
    // constructor
    feature_3D_info( )
    {
      view_list.clear();
      descriptors.clear();
    }
    
    // copy constructor
    feature_3D_info( const feature_3D_info &other )
    {
      view_list.clear();
      descriptors.clear();
      view_list.resize( other.view_list.size() );
      descriptors.resize( other.descriptors.size() );
      
      for( uint32_t i=0; i< (uint32_t) other.view_list.size(); ++i )
      {
        view_list[i].camera = other.view_list[i].camera;
        view_list[i].key = other.view_list[i].key;
        view_list[i].x = other.view_list[i].x;
        view_list[i].y = other.view_list[i].y;
        view_list[i].scale = other.view_list[i].scale;
        view_list[i].orientation = other.view_list[i].orientation;
      }
      
      for( uint32_t i=0; i<(uint32_t)other.descriptors.size(); ++i )
        descriptors[i] = other.descriptors[i];
      
      point = other.point;
    }
    
    // destructor
    ~feature_3D_info( )
    {
      view_list.clear();
      descriptors.clear();
    }
    
    // reset all data
    void clear_data()
    {
      view_list.clear();
      descriptors.clear();
    }
};

////////////////////////////
// Class to parse the output files generated by Bundler.
// The class is able to further load the .key files of all cameras in
// the reconstruction given that the .key files are not zipped (Bundler
// generally uses gzip to compress the .key files). It can also
// load the file format generated by Bundle2Info.
////////////////////////////
class parse_bundler
{
  public:
    
    // standard constructor
    parse_bundler();
    
    // standard destructor
    ~parse_bundler( );
    
    // get the number of 3D points in the reconstruction
    uint32_t get_number_of_points( );

    // get the number of cameras in the reconstruction
    uint32_t get_number_of_cameras( );
    
    // get the actual 3D points from the loaded reconstruction
    void get_points( std::vector< point3D > &points );
    
    // get the feature information extracted from the reconstruction
    std::vector< feature_3D_info >& get_feature_infos( );
    
    // get the cameras included in the reconstruction
    std::vector< bundler_camera >& get_cameras( );
    
    
    // Parses the output text file generated by bundler.
    // Input parameters: The filename of the output file (usually bundle.out) and
    // the filename of the list of images (usually list.txt).
    bool parse_data( const char* bundle_out_filename_, const char* image_list_filename );
    
    // Load the information from a binary file constructed with Bundle2Info.
    // The format parameter specifies whether the binary file contains 
    // no camera information (format=0, as generated by Bundle2Info) or does 
    // contain camera information (format=1, for example the aachen.info file 
    // released with the Aachen dataset from the paper
    //
    // Torsten Sattler, Tobias Weyand, Bastian Leibe, Leif Kobbelt. 
    // Image Retrieval for Image-Based Localization Revisited.
    // BMVC 2012.
    bool load_from_binary( const char* filename, const int format );
    
    // clear the loaded data
    void clear();
    
  private:
    
    std::vector< bundler_camera > mCameras;
    
    std::vector< feature_3D_info > mFeatureInfos;
    
    uint32_t mNbPoints, mNbCameras;
    
};

#endif 
