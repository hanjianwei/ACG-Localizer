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


#include "parse_bundler.hh"


//------------------------------    

parse_bundler::parse_bundler()
{
  mNbPoints = 0;
  mNbCameras = 0;
}

//------------------------------    

parse_bundler::~parse_bundler( )
{
  for( uint32_t i=0; i<mNbPoints; ++i )
  {
    mFeatureInfos[i].view_list.clear();
    mFeatureInfos[i].descriptors.clear();
  }
  mFeatureInfos.clear();
  
}


//------------------------------    

uint32_t parse_bundler::get_number_of_points( )
{
  return mNbPoints;
}

//------------------------------ 

uint32_t parse_bundler::get_number_of_cameras( )
{
  return mNbCameras;
}

//------------------------------    

void parse_bundler::get_points( std::vector< point3D > &points )
{
  points.reserve( mNbPoints );
  points.clear();
  for( uint32_t i=0; i<mNbPoints; ++i )
    points.push_back( point3D( mFeatureInfos[i].point ) );
}

//------------------------------    

std::vector< feature_3D_info >& parse_bundler::get_feature_infos( )
{
  return mFeatureInfos;  
}

//------------------------------  

std::vector< bundler_camera >& parse_bundler::get_cameras( )
{
  return mCameras;
}

//------------------------------    

bool parse_bundler::parse_data( const char* bundle_out_filename_, const char* image_list_filename )
{
  ////
  // intialize the points and their views
  for( uint32_t i=0; i<mNbPoints; ++i )
  {
    mFeatureInfos[i].view_list.clear();
    mFeatureInfos[i].descriptors.clear();
  }
  mFeatureInfos.clear();
  
  mNbCameras = mNbPoints = 0;
  

  ////
  // load the Bundler file
  
  std::cout << " Parsing " << bundle_out_filename_ << std::endl;
  
  std::ifstream instream( bundle_out_filename_, std::ios::in );
  
  if ( !instream.is_open() )
  {
    std::cerr << " Could not open the file " << bundle_out_filename_ << std::endl;
    return false;
  }
  
  // read the first line (containing only some information about the Bundler version)
  char header[4096];
  instream.getline( header, 4086, '\n' );
  std::cout << "  header of the file: " << header << std::endl;
  
  // get the number of cameras and the number of 3D points
  instream >> mNbCameras >> mNbPoints;
  
  mFeatureInfos.resize( mNbPoints );
  mCameras.resize( mNbCameras );
  
  // load the camera data
  std::cout << " skipping " << mNbCameras << " cameras" << std::endl;

  float camera_data = 0.0f;
  for( uint32_t i=0; i<mNbCameras; ++i )
  {
    instream >> mCameras[i].focal_length >> mCameras[i].kappa_1 >> mCameras[i].kappa_2;
    for( int j=0; j<3; ++j )
      instream >> mCameras[i].rotation( j,0 ) >> mCameras[i].rotation( j,1 ) >> mCameras[i].rotation( j,2 );
    instream >> mCameras[i].translation[0] >> mCameras[i].translation[1] >> mCameras[i].translation[2];
    mCameras[i].id = i;
  }
  
  std::cout << "   done " << std::endl;
  
  // load the points ...
  std::cout << " starting to load " << mNbPoints << " 3D points" << std::endl;
  int r,g,b;
  
  // .. and store for each camera which points ( index of the point in mFeatureInfos and the cameras id in the 
  // view list) it does see. This is needed to later on load the descriptors from the .key files
  std::vector< std::vector< std::pair< uint32_t, uint32_t > > > cam_feature_infos( mNbCameras );
  for( uint32_t i=0; i<mNbCameras; ++i )
    cam_feature_infos[i].clear();
  
  // read the 3D points together with their view lists
  for( uint32_t i=0; i<mNbPoints; ++i )
  {
    // read the position
    instream >> mFeatureInfos[i].point.x >> mFeatureInfos[i].point.y >> mFeatureInfos[i].point.z;
    
    // read the color of the point
    instream >> r >> g >> b;
    mFeatureInfos[i].point.r = (unsigned char) r;
    mFeatureInfos[i].point.g = (unsigned char) g;
    mFeatureInfos[i].point.b = (unsigned char) b;
    
    
    // now, read the view list
    uint32_t view_list_length;
    instream >> view_list_length;
    
    mFeatureInfos[i].view_list.resize(view_list_length);
    mFeatureInfos[i].descriptors.resize( 128*view_list_length, 0 );
    
    for( uint32_t j=0; j<view_list_length; ++j )
    {
      instream >> mFeatureInfos[i].view_list[j].camera >> mFeatureInfos[i].view_list[j].key >> mFeatureInfos[i].view_list[j].x >> mFeatureInfos[i].view_list[j].y;
      cam_feature_infos[mFeatureInfos[i].view_list[j].camera].push_back( std::make_pair( i, j ) );
    }
  }
  
  instream.close();
  
  std::cout << "   done " << std::endl;
  
  ////
  // if no file containing the filenames of the images in the reconstruction is specified, the we are done with loading
  if( image_list_filename == 0 )
    return true;
  
  ////
  // read the file containing the list of images and extract the names of the keyfiles from it
  std::cout << " extracting names of the .key files from " << image_list_filename << std::endl;
  
  std::vector< std::string > keyfilenames;
  keyfilenames.resize( mNbCameras );
  {    
    std::ifstream instream2( image_list_filename, std::ios::in );
    
    if ( !instream2.is_open() )
    {
      std::cerr << " Could not open the file " << image_list_filename << std::endl;
      return false;
    }
    
    char buffer[8192];
    
    for( uint32_t i=0; i<mNbCameras; ++i )
    {
      instream2.getline( buffer, 8192, '\n' );
      std::string strbuffer( buffer );
      std::stringstream sstream( strbuffer );
      
      std::string tmp;
      sstream >> tmp;
      
      // replace the .jpg ending with .key
      
      tmp.replace(tmp.size()-3,3,"key");
      keyfilenames[i] = tmp;
    }
    
    instream2.close();
  }
  
  std::cout << "   done " << std::endl;
  
  ////
  // now load the key-files containing the SIFT keys one for one
  std::cout << " starting to load the keypoints from " << mNbCameras << " many images " << std::endl;
  
  // keep track on how many keypoints could not be found
  uint32_t missing_keypoints = 0;
  
  for( uint32_t i=0; i<mNbCameras; ++i )
  {
    // load the .key file for that camera
  
    SIFT_loader key_loader;
    key_loader.load_features( keyfilenames[i].c_str(), LOWE );
    
    std::vector< unsigned char* >& descriptors = key_loader.get_descriptors();
    std::vector< SIFT_keypoint >& keypoints = key_loader.get_keypoints();
        
    // go through the descriptors and store the ones we are interested in
    uint32_t nb_des = (uint32_t) cam_feature_infos[i].size();
    
    for( uint32_t j=0; j<nb_des; ++j )
    {
      uint32_t feature_id = cam_feature_infos[i][j].first;
      uint32_t view_list_id = cam_feature_infos[i][j].second;
      uint32_t feature_descriptor_index = view_list_id * 128;
      uint32_t feature_id_keyfile = mFeatureInfos[feature_id].view_list[view_list_id].key;
      // copy the descriptor
      
      if( feature_id_keyfile >= (uint32_t) keypoints.size() )
      {
        std::cerr << " Trying to load the descriptor for feature " << feature_id << " from camera " << i << " viewlist entry " << view_list_id << " feature id in the keyfile: " << feature_id_keyfile << " but only " << keypoints.size() << " in keyfile " << std::endl;
        
        ++missing_keypoints;
    
        for( uint32_t k=0; k<128; ++k )
          mFeatureInfos[feature_id].descriptors[feature_descriptor_index+k] = 0;
        
        // set scale and orientation of the keypoint in the view list
        mFeatureInfos[feature_id].view_list[view_list_id].scale = -1.0f;
        mFeatureInfos[feature_id].view_list[view_list_id].orientation = 0.0f;
      }
      else
      {
        for( uint32_t k=0; k<128; ++k )
          mFeatureInfos[feature_id].descriptors[feature_descriptor_index+k] = descriptors[feature_id_keyfile][k];
    
        // set scale and orientation of the keypoint in the view list
        mFeatureInfos[feature_id].view_list[view_list_id].scale = keypoints[feature_id_keyfile].scale;
        mFeatureInfos[feature_id].view_list[view_list_id].orientation = keypoints[feature_id_keyfile].orientation;
      }
    }
    
    // clean up the reserved space
    uint32_t nb_keys_file = (uint32_t) key_loader.get_nb_features();
    
    for( uint32_t j=0; j<nb_keys_file; ++j )
    {
      if( descriptors[j] != 0 )
        delete [] descriptors[j];
      descriptors[j] = 0;
    }
    
    descriptors.clear();
    keypoints.clear();
    
    std::cout << "   " << i+1 << " / " << mNbCameras << std::endl;
  }
  
  std::cout << " Could not find " << missing_keypoints << " many keypoints" << std::endl;
  
  keyfilenames.clear();
  std::cout << "   done " << std::endl;
  
  return true;
}


//------------------------------    

bool parse_bundler::load_from_binary( const char* filename, const int format )
{
  ////
  // initialize the data structure for the 3D points
  for( uint32_t i=0; i<mNbPoints; ++i )
  {
    mFeatureInfos[i].view_list.clear();
    mFeatureInfos[i].descriptors.clear();
  }
  mFeatureInfos.clear();
  
    
  mNbCameras = mNbPoints = 0;
  
  // open file for reading
  std::ifstream ifs( filename, std::ios::in | std::ios::binary );
  if ( !ifs )
  {
    std::cerr << "Cannot read file " << filename << std::endl;
    return false;
  }
  
  // read the number of cameras
  ifs.read(( char* ) &mNbCameras, sizeof( uint32_t ) );
  
  // depending on the format, cameras will be loaded
  if( format == 1 )
  {
    // load the cameras
    mCameras.resize(mNbCameras);
    for( uint32_t i=0; i<mNbCameras; ++i )
    {
      double focal_length, kappa_1, kappa_2;
      int32_t width, height;
      double *rotation = new double[9];
      double *translation = new double[3];
      ifs.read( (char* ) &focal_length, sizeof( double ) );
      ifs.read( (char* ) &kappa_1, sizeof( double ) );
      ifs.read( (char* ) &kappa_2, sizeof( double ) );
      ifs.read( (char* ) &width, sizeof( int32_t ) );
      ifs.read( (char* ) &height, sizeof( int32_t ) );
      ifs.read( (char* ) rotation, 9*sizeof( double ) );
      ifs.read( (char* ) translation, 3*sizeof( double ) );
      
      mCameras[i].focal_length = focal_length;
      mCameras[i].kappa_1 = kappa_1;
      mCameras[i].kappa_2 = kappa_2;
      mCameras[i].width = width;
      mCameras[i].height = height;
      for( int j=0; j<3; ++j )
      {
        for( int k=0; k<3; ++k )
          mCameras[i].rotation( j, k ) = rotation[3*j+k];
      }
      for( int j=0; j<3; ++j )
        mCameras[i].translation[j] = translation[j];
      
      delete [] rotation;
      delete [] translation;
    }
  }
  
  // load the points
  ifs.read(( char* ) &mNbPoints, sizeof( uint32_t ) );
  mFeatureInfos.resize(mNbPoints);
  
  for( uint32_t i=0; i<mNbPoints; ++i )
  {
    float *pos = new float[3];
    ifs.read( (char* ) pos, 3*sizeof( float ) );
    mFeatureInfos[i].point.x = pos[0];
    mFeatureInfos[i].point.y = pos[1];
    mFeatureInfos[i].point.z = pos[2];
    delete [] pos;
    pos = 0;

    uint32_t size_view_list=0;
    ifs.read(( char* ) &size_view_list, sizeof( uint32_t ) );
    mFeatureInfos[i].view_list.resize(size_view_list);
    mFeatureInfos[i].descriptors.resize( 128*size_view_list, 0 );
    
    unsigned char *desc = new unsigned char[128];
    for( uint32_t j=0; j<size_view_list; ++j )
    {
      float x,y,scale, orientation;
      uint32_t cam_id;
      ifs.read( (char* ) &cam_id, sizeof( uint32_t ) );
      ifs.read( (char* ) &x, sizeof( float ) );
      ifs.read( (char* ) &y, sizeof( float ) );
      ifs.read( (char* ) &scale, sizeof( float ) );
      ifs.read( (char* ) &orientation, sizeof( float ) );
      mFeatureInfos[i].view_list[j].camera = cam_id;
      mFeatureInfos[i].view_list[j].x = x;
      mFeatureInfos[i].view_list[j].y = y;
      mFeatureInfos[i].view_list[j].scale = scale;
      mFeatureInfos[i].view_list[j].orientation = orientation;
      
      ifs.read( (char* ) desc, 128*sizeof( unsigned char ) );
      
      // store the descriptor
      for( uint32_t k=0; k<128; ++k )
        mFeatureInfos[i].descriptors[ 128*j+k] = desc[k];
    }
    delete [] desc;
  }
  
  ifs.close();
  
  return true;
}




//------------------------------    


void parse_bundler::clear()
{
 
  mNbCameras = 0;
  
  mNbPoints = (uint32_t) mFeatureInfos.size();
  
  for( uint32_t i=0; i<mNbPoints; ++i )
  {
    mFeatureInfos[i].view_list.clear();
    mFeatureInfos[i].descriptors.clear();
  }  
  mFeatureInfos.clear();
  
  mCameras.clear();
  
}