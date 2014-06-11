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


#define __STDC_LIMIT_MACROS

// C++ includes
#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <stdint.h>
#include <string>
#include <algorithm>
#include <climits>
#include <float.h>
#include <cmath>
#include <sstream>
#include <stdio.h>

// includes for classes dealing with SIFT-features
#include "features/SIFT_loader.hh"
#include "features/visual_words_handler.hh"

// stopwatch
#include "timer.hh"

// math functionality
#include "math/projmatrix.hh"
#include "math/matrix3x3.hh"

// RANSAC
#include "RANSAC.hh"

// exif reader to get the width and height out 
// of the exif tags of an image
#include "exif_reader/exif_reader.hh"

// simple vector class for 3D points
#include <OpenMesh/Core/Geometry/VectorT.hh>

const uint64_t sift_dim = 128;

////
// Classes to handle the two nearest neighbors (nn) of a descriptor.
// There are three classes:
// 1. Normal 2 nearest neighbors for integer distances
// 2. 2 nearest neighbors for integer distances, making sure 
//    that the second nearest neighbor does not belong to the same 3D point
// 3. Normal 2 nearest neighbors for floating point distances
//
// We store the distances to the 2 nearest neighbors as well as the ids of the 
// corresponding 3D points and update the 2 nearest neighbors if needed.
// The stored distances are squared Euclidean distances.
////  

// for unsigned char descriptors
class nearest_neighbors
{
  public:

    // constructors
    nearest_neighbors()
    {
      nn_idx1 = nn_idx2 = UINT32_MAX;
      dist1 = dist2 = -1;
    }
    
    nearest_neighbors( uint32_t nn1, uint32_t nn2, int d1, int d2 ) : nn_idx1( nn1 ), nn_idx2( nn2 ), dist1( d1 ), dist2( d2 )
    {}
      
    nearest_neighbors( uint32_t nn1, int d1 ) : nn_idx1( nn1 ), nn_idx2( UINT32_MAX ), dist1( d1 ), dist2( -1 )
    {}
    
    nearest_neighbors( const nearest_neighbors &other )
    {
      if( &other != this )
      {
        nn_idx1 = other.nn_idx1;
        nn_idx2 = other.nn_idx2;
        dist1 = other.dist1;
        dist2 = other.dist2;
      }
    }
    
    // update the 2 nn with a new distance to a 3D points 
    void update( uint32_t point, int dist )
    {
      if( dist1 < 0 )
      {
        nn_idx1 = point;
        dist1 = dist;
      }
      else
      {
        if( dist < dist1 )
        {
          nn_idx2 = nn_idx1;
          dist2 = dist1;
          nn_idx1 = point;
          dist1 = dist;
        }
        else if( dist < dist2 || dist2 < 0 )
        {
          nn_idx2 = point;
          dist2 = dist;
        }
      }
    }
    
    float get_ratio()
    {
      return float(dist1) / float(dist2);
    }
    
    uint32_t nn_idx1, nn_idx2;
    int dist1, dist2;
};

// for the case that multiple descriptors of the same 3D point are mapped to the same visual word
class nearest_neighbors_multiple
{
  public:
    nearest_neighbors_multiple()
    {
      nn_idx1 = nn_idx2 = UINT32_MAX;
      dist1 = dist2 = -1;
    }
    
    nearest_neighbors_multiple( uint32_t nn1, uint32_t nn2, int d1, int d2 ) : nn_idx1( nn1 ), nn_idx2( nn2 ), dist1( d1 ), dist2( d2 )
    {}
      
    nearest_neighbors_multiple( uint32_t nn1, int d1 ) : nn_idx1( nn1 ), nn_idx2( UINT32_MAX ), dist1( d1 ), dist2( -1 )
    {}
    
    nearest_neighbors_multiple( const nearest_neighbors_multiple &other )
    {
      if( &other != this )
      {
        nn_idx1 = other.nn_idx1;
        nn_idx2 = other.nn_idx2;
        dist1 = other.dist1;
        dist2 = other.dist2;
      }
    }
    
    void update( uint32_t point, int dist )
    {
      if( dist1 < 0 )
      {
        nn_idx1 = point;
        dist1 = dist;
      }
      else
      {
        if( dist < dist1 )
        {
          if( nn_idx1 != point )
          {
            nn_idx2 = nn_idx1;
            dist2 = dist1;
          }
          nn_idx1 = point;
          dist1 = dist;
        }
        else if( (dist < dist2 || dist2 < 0) && (point != nn_idx1) )
        {
          nn_idx2 = point;
          dist2 = dist;
        }
      }
    }
    
    float get_ratio()
    {
      return float(dist1) / float(dist2);
    }
    
    uint32_t nn_idx1, nn_idx2;
    int dist1, dist2;
};

// for float descriptors
class nearest_neighbors_float
{
  public:
    nearest_neighbors_float()
    {
      nn_idx1 = nn_idx2 = UINT32_MAX;
      dist1 = dist2 = -1.0;
    }
    
    nearest_neighbors_float( uint32_t nn1, uint32_t nn2, float d1, float d2 ) : nn_idx1( nn1 ), nn_idx2( nn2 ), dist1( d1 ), dist2( d2 )
    {}
      
    nearest_neighbors_float( uint32_t nn1, float d1 ) : nn_idx1( nn1 ), nn_idx2( UINT32_MAX ), dist1( d1 ), dist2( -1 )
    {}
    
    nearest_neighbors_float( const nearest_neighbors_float &other )
    {
      if( &other != this )
      {
        nn_idx1 = other.nn_idx1;
        nn_idx2 = other.nn_idx2;
        dist1 = other.dist1;
        dist2 = other.dist2;
      }
    }
    
    void update( uint32_t point, float dist )
    {
      if( dist1 < 0 )
      {
        nn_idx1 = point;
        dist1 = dist;
      }
      else
      {
        if( dist < dist1 )
        {
          nn_idx2 = nn_idx1;
          dist2 = dist1;
          nn_idx1 = point;
          dist1 = dist;
        }
        else if( dist < dist2 || dist2 < 0 )
        {
          nn_idx2 = point;
          dist2 = dist;
        }
      }
    }
    
    float get_ratio()
    {
      return dist1 / dist2;
    }
    
    uint32_t nn_idx1, nn_idx2;
    float dist1, dist2;
};

////
// functions to compute the squared distances between two SIFT-vectors
// there are different ways how the SIFT-vectors are stored, for each there is 
// one function
////

// First descriptor is stored in an array, while the second descriptor is stored in a vector (concatenation of vector entries)
// The second descriptor begins at position index*128
inline int compute_squared_SIFT_dist( const unsigned char * const v1, std::vector< unsigned char > &v2, uint32_t index )
{
  uint64_t index_( index );
  index_ *= sift_dim;
  int dist = 0;
  int x = 0;
  for( uint64_t i=0; i<sift_dim; ++i )
  {
    x = int( v1[i] ) - int( v2[index_+i] );
    dist += x*x;
  }
  return dist;
}

// same in case that one descriptors consists of floating point values
inline float compute_squared_SIFT_dist_float( const unsigned char * const v1, std::vector< float > &v2, uint32_t index )
{
  size_t index_( index );
  index_ *= sift_dim;
  float dist = 0;
  float x = 0;
  for( int i=0; i<sift_dim; ++i )
  {
    x = float( v1[i] ) - v2[index_+i];
    dist += x*x;
  }
  return dist;
}


// function to sort (2D feature, visual word) point pairs for the prioritized search.
inline bool cmp_priorities( const std::pair< uint32_t, uint32_t >& a, const std::pair< uint32_t, uint32_t >& b )
{
  return ( a.second < b.second );
}

////
// constants
////

// minimal number of inliers required to accept an image as registered
uint32_t minimal_RANSAC_solution = 12;

// SIFT-ratio value. Since we store squared distances, we need the squared value 0.7^2 = 0.49
float nn_ratio = 0.49f; 

// the assumed minimal inlier ratio of RANSAC
float min_inlier = 0.2f;

// stop RANSAC if 60 seconds have passed
double ransac_max_time = 60.0; 

//---------------------------------------------------------------------------------------------------------------------------------------------------------------

////
// Actual localization method
////

int main (int argc, char **argv)
{
  if( argc < 10 )
  {
    std::cout << "____________________________________________________________________________________________________________________________" << std::endl;
    std::cout << " -                                                                                                                        - " << std::endl;
    std::cout << " -        Localization method. Implementation of the localization framework proposed in the ICCV 2011 paper               - " << std::endl;
    std::cout << " -          T. Sattler, B. Leibe, L. Kobbelt. Fast Image-Based Localization using Direct 2D-to-3D Matching.               - " << std::endl;
    std::cout << " -                               2011 by Torsten Sattler (tsattler@cs.rwth-aachen.de)                                     - " << std::endl;
    std::cout << " -                                                                                                                        - " << std::endl;
    std::cout << " - usage: acg_localizer list nb_trees nb_cluster clusters descriptors mode in_ratio max_corr results                      - " << std::endl;
    std::cout << " - Parameters:                                                                                                            - " << std::endl;
    std::cout << " -  list                                                                                                                  - " << std::endl;
    std::cout << " -     List containing the filenames of all the .key files that should be used as query. It is assumed that the           - " << std::endl;
    std::cout << " -     corresponding images have the same filename except of ending in .jpg.                                              - " << std::endl;
    std::cout << " -                                                                                                                        - " << std::endl;
    std::cout << " -  nb_tree                                                                                                               - " << std::endl;
    std::cout << " -     The number of trees to use for the random kd-forest that is employed to do the assignments of 2D features to       - " << std::endl;
    std::cout << " -     visual words. Set to 1 for all experiments in the paper.                                                           - " << std::endl;
    std::cout << " -                                                                                                                        - " << std::endl;
    std::cout << " -  nb_cluster                                                                                                            - " << std::endl;
    std::cout << " -     The number of clusters in the file containing the cluster centers. The number of cluster must be the same          - " << std::endl;
    std::cout << " -     as for the assignments computed by compute_desc_assignments.                                                       - " << std::endl;
    std::cout << " -                                                                                                                        - " << std::endl;
    std::cout << " -  clusters                                                                                                              - " << std::endl;
    std::cout << " -     The cluster centers (visual words), stored in a textfile consisting of nb_clusters * 128 floating point values.    - " << std::endl;
    std::cout << " -                                                                                                                        - " << std::endl;
    std::cout << " -  descriptors                                                                                                           - " << std::endl;
    std::cout << " -     The assignments assigning descriptors (and 3D points) to visual words, computed by compute_desc_assignments.       - " << std::endl;
    std::cout << " -                                                                                                                        - " << std::endl;
    std::cout << " -  mode                                                                                                                  - " << std::endl;
    std::cout << " -     The way the descriptors in the assignments file are stored (0 = unsigned char, 1 = float). Set mode to 2 to        - " << std::endl;
    std::cout << " -     allow multiple descriptors (unsigned char) of the same point to be mapped to the same visual word.                 - " << std::endl;
    std::cout << " -                                                                                                                        - " << std::endl;
    std::cout << " -  in_ratio                                                                                                              - " << std::endl;
    std::cout << " -     Assumed minimal inlier ratio for RANSAC.                                                                           - " << std::endl;
    std::cout << " -                                                                                                                        - " << std::endl;
    std::cout << " -  max_cor                                                                                                               - " << std::endl;
    std::cout << " -     If max_cor correspondences are found, the search for more correspondences is terminated.                           - " << std::endl;
    std::cout << " -     Set to 0 to disable early termination.                                                                             - " << std::endl;
    std::cout << " -                                                                                                                        - " << std::endl;
    std::cout << " -  results                                                                                                               - " << std::endl;
    std::cout << " -     The program will write the results of the localization into a text file of name \"results\". It has the following  - " << std::endl;
    std::cout << " -     format, where every line in the file belongs to one query image and has the format                                 - " << std::endl;
    std::cout << " -       #inliers #(correspondences found) (time needed to compute the visual words, in seconds) (time needed for linear  - " << std::endl;
    std::cout << " -       search, in seconds) (time needed for RANSAC, in seconds) (total time needed, in seconds)                         - " << std::endl;
    std::cout << " -                                                                                                                        - " << std::endl;
    std::cout << "____________________________________________________________________________________________________________________________" << std::endl;
    return -1;
  }
  
  ////
  // get the parameters
  std::string keylist( argv[1] );
  size_t nb_trees = (size_t) atoi( argv[2] );
  uint32_t nb_clusters = (uint32_t) atoi( argv[3] );
  std::string cluster_file( argv[4] );
  std::string vw_assignments( argv[5] );
  int mode = atoi( argv[6] );
  if( mode < 0 || mode > 2 )
  {
    std::cerr << " ERROR: unknown mode " << mode << std::endl;
    return -1;
  }
  
  
  min_inlier = atof( argv[7] );
  
  std::cout << " Assumed minimal inlier-ratio: " << min_inlier << std::endl;
  
  size_t max_cor_early_term = (size_t) atoi( argv[8] );
  
  if( max_cor_early_term > 0 )
    std::cout << " Early termination after finding " << max_cor_early_term << " correspondences " << std::endl;
  else
    std::cout << " No early termination in correspondence search" << std::endl;
  
  std::string results( argv[9] );
  
  ////
  // create and open the output file
  std::ofstream ofs_details( results.c_str(), std::ios::out );
    
  if( !ofs_details.is_open() )
  {
    std::cerr << " Could not write results to " << results << std::endl;
    return 1;
  }
  
  ////
  // load the visual words and their tree 
  visual_words_handler vw_handler;
  vw_handler.set_nb_trees( nb_trees );
  vw_handler.set_nb_visual_words( nb_clusters );
  vw_handler.set_branching( 10 );

  vw_handler.set_method(std::string("flann"));
  vw_handler.set_flann_type(std::string("randomkd"));
  if( !vw_handler.create_flann_search_index( cluster_file ) )
  {
    std::cout << " ERROR: Could not load the cluster centers from " << cluster_file << std::endl;;
    return -1;
  } 
  std::cout << "  done " << std::endl;
    
  
  ////
  // load the assignments for the visual words
  
  std::cout << "* Loading and parsing the assignments ... " << std::endl;
  
  // store the 3D positions of the 3D points
  std::vector< OpenMesh::Vec3f > points3D;
  
  // store the descriptors in a vector simply by concatenating their entries
  // depending on the mode, either unsigned char entries or floating point entries are used
  std::vector< unsigned char > all_descriptors;
  std::vector< float > all_descriptors_float;
  
  // for every visual word, store a vector of (3D point id, descriptor id) pairs, where the 3D point id
  // is the index of the 3D point in points3D and the descriptor id is the position of the first entry of the corresponding 
  // descriptor in all_descriptors / all_descriptors_float
  std::vector< std::vector< std::pair< uint32_t, uint32_t > > > vw_points_descriptors(nb_clusters);
  
  // store per visual word the number of (point, descriptor) pairs store in it
  std::vector< uint32_t > nb_points_per_vw(nb_clusters,0);
  
  // number of non-empty visual words, the number of 3D points and the total number of descriptors
  uint32_t nb_non_empty_vw, nb_3D_points, nb_descriptors;
  
  for( uint32_t i=0; i<nb_clusters; ++i )
    vw_points_descriptors[i].clear();
  
  // load the assignments from a file generated by compute_desc_assignments
  {
    std::ifstream ifs( vw_assignments.c_str(), std::ios::in | std::ios::binary  );
    
    if ( !ifs )
    {
      std::cerr << " ERROR: Cannot read the visual word assignments " << vw_assignments << std::endl;
      return -1;
    }
    
    uint32_t nb_clusts;
    ifs.read(( char* ) &nb_3D_points, sizeof( uint32_t ) );
    ifs.read(( char* ) &nb_clusts, sizeof( uint32_t ) );
    ifs.read(( char* ) &nb_non_empty_vw, sizeof( uint32_t ) );
    ifs.read(( char* ) &nb_descriptors, sizeof( uint32_t ) );

    if( nb_clusts != nb_clusters )
      std::cerr << " WARNING: Number of clusters differs! " << nb_clusts << " " << nb_clusters << std::endl;
    
    std::cout << "  Number of non-empty clusters: " << nb_non_empty_vw << " number of points : " << nb_3D_points << " number of descriptors: " << nb_descriptors << std::endl;
    
    // read the 3D points and their visibility polygons
    points3D.resize(nb_3D_points);
    if( mode == 0 || mode == 2 )
      all_descriptors.resize(128*nb_descriptors);
    else
      all_descriptors_float.resize(128*nb_descriptors);
    
    
    // load the points
    float *point_data = new float[3];
    for( uint32_t i=0; i<nb_3D_points; ++i )
    {
      ifs.read(( char* ) point_data, 3 * sizeof( float ) );
      for( int j=0; j<3; ++j )
        points3D[i][j] = point_data[j];
    }
    delete [] point_data;
       
    // load the descriptors
    int tmp_int;
    uint64_t index = 0;
    for( uint32_t i=0; i<nb_descriptors; ++i, index += sift_dim )
    {
      for( uint64_t j=0; j<sift_dim; ++j )
      {
        if( mode == 0 || mode == 2 )
          ifs.read(( char* ) &all_descriptors[index+j], sizeof( unsigned char ) );
        else
          ifs.read(( char* ) &all_descriptors_float[index+j], sizeof( float ) );
      }
    }
    
    // now we load the assignments of the pairs (point_id, descriptor_id) to the visual words
    for( uint32_t i=0; i<nb_non_empty_vw; ++i )
    {
      uint32_t id, nb_pairs;
      ifs.read(( char* ) &id, sizeof( uint32_t ) );
      ifs.read(( char* ) &nb_pairs, sizeof( uint32_t ) );
      vw_points_descriptors[id].resize( nb_pairs );
      nb_points_per_vw[id] = nb_pairs;
      for( uint32_t j=0; j<nb_pairs; ++j )
      {
        ifs.read(( char* ) &vw_points_descriptors[id][j].first, sizeof( uint32_t ) );
        ifs.read(( char* ) &vw_points_descriptors[id][j].second, sizeof( uint32_t ) );
      }
    }
    ifs.close();
    
    std::cout << "  done loading and parsing the assignments " << std::endl;
  }
  

  
  ////
  // now load all the filenames of the query images
  std::vector< std::string > key_filenames;
  key_filenames.clear();
  {
    std::ifstream ifs( keylist.c_str(), std::ios::in );
    std::string tmp_string;
    
    while( !ifs.eof() )
    {
      tmp_string = "";
      ifs >> tmp_string;
      if( !tmp_string.empty() )
        key_filenames.push_back(tmp_string);
    }
    ifs.close();
    std::cout << " done loading " << key_filenames.size() << " keyfile names " << std::endl;
  }
  
  uint32_t nb_keyfiles = key_filenames.size();
  
  ////
  // do the actual localization
  
  // first initialize some variables used to compute statistics about the number of registered images
  double avrg_reg_time = 0.0;
  double avrg_reject_time = 0.0;
  double avrg_vw_time = 0.0;
  double avrg_nb_features = 0.0;
  double avrg_cor_computation_time_accepted = 0.0;
  double avrg_cor_computation_time_rejected = 0.0;
  double N = 0.0;
  double N_reject = 0.0;
  double mean_inlier_ratio_accepted = 0.0;
  double mean_inlier_ratio_rejected = 0.0;
  double mean_nb_correspondences_accepted = 0.0;
  double mean_nb_correspondences_rejected = 0.0;
  double mean_nb_features_accepted = 0.0;
  double mean_nb_features_rejected = 0.0;
  double vw_time = 0.0;
  double corr_time = 0.0;
  double RANSAC_time = 0.0;
  double avrg_RANSAC_time_registered = 0.0;
  double avrg_RANSAC_time_rejected = 0.0;
  
  // the number of registered images
  uint32_t registered = 0;
  
  // store all assignments of 2D features to visual words in one large vector (resized if necessary)
  // preallocated for speed
  std::vector< uint32_t > computed_visual_words( 50000, 0 );
  
  
  for( uint32_t i=0; i<nb_keyfiles; ++i, N+=1.0 )
  {
    std::cout << std::endl << " --------- " << i+1 << " / " << nb_keyfiles << " --------- " << std::endl;
    
    // load the features
    SIFT_loader key_loader;
    key_loader.load_features( key_filenames[i].c_str(), LOWE );
    
    std::vector< unsigned char* >& descriptors = key_loader.get_descriptors();
    std::vector< SIFT_keypoint >& keypoints = key_loader.get_keypoints();
    
    uint32_t nb_loaded_keypoints = (uint32_t) keypoints.size();
    
   
    // center the keypoints around the center of the image
    // first we need to get the dimensions of the image which we obtain from its exif tag
    int img_width, img_height;
    std::string jpg_filename( key_filenames[i] );
    jpg_filename.replace( jpg_filename.size()-3,3,"jpg");
    exif_reader::open_exif( jpg_filename.c_str() );
    img_width = exif_reader::get_image_width();
    img_height = exif_reader::get_image_height();
    exif_reader::close_exif();
    
    for( uint32_t j=0; j<nb_loaded_keypoints; ++j )
    {
      keypoints[j].x -= (img_width-1.0)/2.0f; 
      keypoints[j].y = (img_height-1.0)/2.0f - keypoints[j].y; 
    }
    
    std::cout << " loaded " << nb_loaded_keypoints << " descriptors from " << key_filenames[i] << std::endl;
    
    ////
    // assign the descriptors to the visual words   
    Timer timer;
    timer.Init();
    timer.Start();
    
    if( computed_visual_words.size() < nb_loaded_keypoints )
      computed_visual_words.resize( nb_loaded_keypoints );

    
    std::set< size_t > unique_vw;
    unique_vw.clear();

    vw_handler.set_nb_paths( 10 );
    vw_handler.assign_visual_words_ucharv( descriptors, nb_loaded_keypoints, computed_visual_words );
    
    timer.Stop();
    
    for( size_t j=0; j<nb_loaded_keypoints; ++j )
      unique_vw.insert( computed_visual_words[j] );
    
    std::cout << " assigned visual words in " << timer.GetElapsedTimeAsString() << " to " << unique_vw.size() << " unique vw" << std::endl;
    
    avrg_vw_time = avrg_vw_time * N / (N+1.0) + timer.GetElapsedTime() / (N+1.0);
    avrg_nb_features = avrg_nb_features * N / (N+1.0) + double(nb_loaded_keypoints) / (N+1.0);
    vw_time = timer.GetElapsedTime();
    
    ////
    // establish 2D-3D correspondences by using the vw to compute pairwise nearest neighbors
    timer.Init();
    timer.Start();

    Timer all_timer;
    all_timer.Init();
    all_timer.Start();
    
    uint32_t max_corr = 0;
    
    // sort the 2D points in ascending order of the number of 3D points belonging to their respective visual words
    std::vector< std::pair< uint32_t, uint32_t > > priorities( nb_loaded_keypoints );
    for( size_t j=0; j<nb_loaded_keypoints; ++j )
    {
      priorities[j].first = j;
      priorities[j].second = nb_points_per_vw[computed_visual_words[j]];
    }
    
    std::sort( priorities.begin(), priorities.end(), cmp_priorities );


    // we store for each 3D point the corresponding 2D feature as well as the squared distance
    // this is needed in case that two 2D features are assigned to one 3D point, because
    // we only want to keep the correspondence to the 2D point with the most similar descriptor
    // i.e. the smallest Euclidean distance in descriptor space
    std::map< uint32_t, std::pair< uint32_t, int > > corr_3D_to_2D;
    corr_3D_to_2D.clear();
    
    std::map< uint32_t, nearest_neighbors >::iterator map_it_2D;
    std::map< uint32_t, nearest_neighbors_float >::iterator map_it_2D_float;
    std::map< uint32_t, nearest_neighbors_multiple >::iterator map_it_2D_multple;
    std::map< uint32_t, std::pair< uint32_t, int > >::iterator map_it_3D;

    // compute nearest neighbors
    // we do a single case distinction wether the database consists of unsigned char descriptors or floating point descriptors
    // while resulting in ugly code, it should improve run-time since we do not need to make the distinction multiple times
    uint32_t nb_comparisons = 0;
    uint32_t nb_considered_points = 0;
    if( mode == 0 )
    {
      for( size_t j=0; j<nb_loaded_keypoints; ++j )
      {
        uint32_t j_index = priorities[j].first;
        uint32_t assignment = uint32_t( computed_visual_words[j_index] );
        
        nearest_neighbors nn;
        
        if( vw_points_descriptors[assignment].size() > 0 )
        {
          // find nearest neighbor for 2D feature, update nearest neighbor information for 3D points if necessary
          size_t nb_poss_assignments = vw_points_descriptors[assignment].size();
          
          nb_comparisons += nb_poss_assignments;
          
          for( size_t k=0; k<nb_poss_assignments; ++k )
          {
            uint32_t point_id = vw_points_descriptors[assignment][k].first;
            uint32_t desc_id = vw_points_descriptors[assignment][k].second;
            
            int dist = compute_squared_SIFT_dist( descriptors[j_index], all_descriptors, desc_id );
            
            nn.update( point_id, dist );
          }
        }
        
        // check if we have found a correspondence
        if( nn.dist1 >= 0 )
        {
          if( nn.dist2 >= 0 )
          {
            if( nn.get_ratio() < nn_ratio )
            {
              // we found one, so we need check for mutual nearest neighbors
              map_it_3D = corr_3D_to_2D.find( nn.nn_idx1 );
        
              if( map_it_3D != corr_3D_to_2D.end() )
              {
                if( map_it_3D->second.second > nn.dist1 )
                {
                  map_it_3D->second.first = j_index;
                  map_it_3D->second.second = nn.dist1;
                }
              }
              else
              {
                corr_3D_to_2D.insert( std::make_pair( nn.nn_idx1, std::make_pair( j_index, nn.dist1 ) ) );
              }
            }
          }
        }
        
        // stop the search if enough correspondences are found
        if( max_cor_early_term > 0 && corr_3D_to_2D.size() >= max_cor_early_term )
        {
          nb_considered_points = j+1;
          break;
        }
      }
    }
    else if( mode == 1 )
    {
      // floating point
      for( size_t j=0; j<nb_loaded_keypoints; ++j )
      {
        uint32_t j_index = priorities[j].first;
        uint32_t assignment = uint32_t( computed_visual_words[j_index] );
        
        nearest_neighbors_float nn;
        
        if( vw_points_descriptors[assignment].size() > 0 )
        {
          // find nearest neighbor for 2D feature, update nearest neighbor information for 3D points if necessary
          size_t nb_poss_assignments = vw_points_descriptors[assignment].size();
          nb_comparisons += nb_poss_assignments;
          
          for( size_t k=0; k<nb_poss_assignments; ++k )
          {
            uint32_t point_id = vw_points_descriptors[assignment][k].first;
            uint32_t desc_id = vw_points_descriptors[assignment][k].second;
            
            float dist = compute_squared_SIFT_dist_float( descriptors[j_index], all_descriptors_float, desc_id );
            
            nn.update( point_id, dist );
          }
        }
        
        // check if we have found a correspondence
        if( nn.dist1 >= 0.0 )
        {
          if( nn.dist2 >= 0.0 )
          {
            if( nn.get_ratio() < nn_ratio )
            {
              // we found one, so we need check for mutual nearest neighbors
              map_it_3D = corr_3D_to_2D.find( nn.nn_idx1 );
        
              if( map_it_3D != corr_3D_to_2D.end() )
              {
                if( map_it_3D->second.second > nn.dist1 )
                {
                  map_it_3D->second.first = j_index;
                  map_it_3D->second.second = nn.dist1;
                }
              }
              else
              {
                corr_3D_to_2D.insert( std::make_pair( nn.nn_idx1, std::make_pair( j_index, nn.dist1 ) ) );
              }
            }
          }
        }
        
        if( max_cor_early_term > 0 && corr_3D_to_2D.size() >= max_cor_early_term )
        {
          nb_considered_points = j+1;
          break;
        }
      }
    }
    else if( mode == 2 )
    {
      for( size_t j=0; j<nb_loaded_keypoints; ++j )
      {
        uint32_t j_index = priorities[j].first;
        uint32_t assignment = uint32_t( computed_visual_words[j_index] );
        
        nearest_neighbors_multiple nn;
        
        if( vw_points_descriptors[assignment].size() > 0 )
        {
          // find nearest neighbor for 2D feature, update nearest neighbor information for 3D points if necessary
          size_t nb_poss_assignments = vw_points_descriptors[assignment].size();
          nb_comparisons += nb_poss_assignments;
          
          for( size_t k=0; k<nb_poss_assignments; ++k )
          {
            uint32_t point_id = vw_points_descriptors[assignment][k].first;
            uint32_t desc_id = vw_points_descriptors[assignment][k].second;
            
            
            int dist = compute_squared_SIFT_dist( descriptors[j_index], all_descriptors, desc_id );
            
            nn.update( point_id, dist );
          }
        }
        
        // check if we have found a correspondence
        if( nn.dist1 >= 0 )
        {
          if( nn.dist2 >= 0 )
          {
            if( nn.get_ratio() < nn_ratio )
            {
              // we found one, so we need check for mutual nearest neighbors
              map_it_3D = corr_3D_to_2D.find( nn.nn_idx1 );
        
              if( map_it_3D != corr_3D_to_2D.end() )
              {
                if( map_it_3D->second.second > nn.dist1 )
                {
                  map_it_3D->second.first = j_index;
                  map_it_3D->second.second = nn.dist1;
                }
              }
              else
              {
                corr_3D_to_2D.insert( std::make_pair( nn.nn_idx1, std::make_pair( j_index, nn.dist1 ) ) );
              }
            }
          }
        }
        
        if( max_cor_early_term > 0 && corr_3D_to_2D.size() >= max_cor_early_term )
        {
          nb_considered_points = j+1;
          break;
        }
      }
    }
    
    if( nb_considered_points == 0 )
      nb_considered_points = nb_loaded_keypoints;
    
    

    ////
    // compute and store the correspondences such that we can easily hand them over to RANSAC
    
    // the 2D and 3D positions of features and points are simply concatenated into 2 vectors
    std::vector< float > c2D, c3D;
    c2D.clear();
    c3D.clear();
    
    // furthermore, we want to store the ids of the 2D features and the 3D points3D
    // first the 2D, then the 3D point
    std::vector< std::pair< uint32_t, uint32_t > > final_correspondences;
    final_correspondences.clear();
    
    // get the correspondences
    for( map_it_3D = corr_3D_to_2D.begin(); map_it_3D != corr_3D_to_2D.end(); ++map_it_3D )
    {
      c2D.push_back(keypoints[map_it_3D->second.first].x);
      c2D.push_back(keypoints[map_it_3D->second.first].y);
      
      c3D.push_back( points3D[map_it_3D->first][0] );
      c3D.push_back( points3D[map_it_3D->first][1] );
      c3D.push_back( points3D[map_it_3D->first][2] );
      
      final_correspondences.push_back( std::make_pair( map_it_3D->second.first, map_it_3D->first ) );
    }
    
    timer.Stop();
    std::cout << " computed correspondences in " << timer.GetElapsedTimeAsString() << ", considering " << nb_considered_points << " features " << std::endl;
    corr_time = timer.GetElapsedTime();
    
    
    ////
    // do the pose verification using RANSAC
    
    uint32_t nb_corr = c2D.size() / 2;
    RANSAC::computation_type = P6pt;
    RANSAC::stop_after_n_secs = true;
    RANSAC::max_time = ransac_max_time;
    RANSAC::error = 10.0f; // for P6pt this is the SQUARED reprojection error in pixels
    RANSAC ransac_solver;
    
    
    std::cout << " applying RANSAC on " << nb_corr << " correspondences " << std::endl;
    timer.Init();
    timer.Start();
    ransac_solver.apply_RANSAC( c2D, c3D, nb_corr, std::max( float( minimal_RANSAC_solution ) / float( nb_corr ), min_inlier ) ); 
    timer.Stop();
    RANSAC_time = timer.GetElapsedTime();
    
    all_timer.Stop();
    
    // output the solution:
    std::cout << "#### found solution ####" << std::endl;
    std::cout << " needed time: " << all_timer.GetElapsedTimeAsString() << std::endl;

    // get the solution from RANSAC
    std::vector< uint32_t > inlier;
    inlier.assign( ransac_solver.get_inliers().begin(), ransac_solver.get_inliers().end()  );
    
    // get the computed projection matrix
    Util::Math::ProjMatrix proj_matrix = ransac_solver.get_projection_matrix();
    
    // decompose the projection matrix
    Util::Math::Matrix3x3 Rot, K;
    proj_matrix.decompose( K, Rot );
    proj_matrix.computeInverse();
    proj_matrix.computeCenter();
    std::cout << " camera calibration: " << K << std::endl;
    std::cout << " camera rotation: " << Rot << std::endl;
    std::cout << " camera position: " << proj_matrix.m_center << std::endl;
        
    ofs_details << inlier.size() << " " << nb_corr << " " << vw_time << " " << corr_time << " " << RANSAC_time << " " << all_timer.GetElapsedTime() << std::endl;
    
    std::cout << "#########################" << std::endl;
    
    // determine whether the image was registered or not
    // also update the statistics about timing, ...
    if( inlier.size() >= minimal_RANSAC_solution )
    {
      double N_reg = registered;
      avrg_reg_time = avrg_reg_time * N_reg / (N_reg+1.0) + all_timer.GetElapsedTime() / (N_reg+1.0);
      mean_inlier_ratio_accepted = mean_inlier_ratio_accepted * N_reg / (N_reg+1.0) + double(inlier.size()) / (double( nb_corr ) * (N_reg+1.0));
      mean_nb_correspondences_accepted = mean_nb_correspondences_accepted * N_reg / (N_reg+1.0) + double(nb_corr) / (N_reg+1.0);
      mean_nb_features_accepted = mean_nb_features_accepted * N_reg / (N_reg+1.0) + double(nb_loaded_keypoints) / (N_reg+1.0);
      avrg_cor_computation_time_accepted = avrg_cor_computation_time_accepted * N_reg /(N_reg+1.0) + corr_time / (N_reg+1.0);
      avrg_RANSAC_time_registered = avrg_RANSAC_time_registered * N_reg / (N_reg+1.0) + RANSAC_time / (N_reg+1.0);
      ++registered;
    }
    else
    {
      avrg_reject_time = avrg_reject_time * N_reject / (N_reject+1.0) + all_timer.GetElapsedTime() / (N_reject+1.0);
      mean_inlier_ratio_rejected = mean_inlier_ratio_rejected * N_reject / (N_reject+1.0) + double(inlier.size()) / (double( nb_corr ) * (N_reject+1.0));
      mean_nb_correspondences_rejected = mean_nb_correspondences_rejected * N_reject / (N_reject+1.0) + double(nb_corr) / (N_reject+1.0);
      mean_nb_features_rejected = mean_nb_features_rejected * N_reject / (N_reject+1.0) + double(nb_loaded_keypoints) / (N_reject+1.0);
      avrg_cor_computation_time_rejected = avrg_cor_computation_time_rejected * N_reject /(N_reject+1.0) + corr_time / (N_reject+1.0);
      avrg_RANSAC_time_rejected = avrg_RANSAC_time_rejected * N_reject / (N_reject+1.0) + RANSAC_time / (N_reject+1.0);
      N_reject += 1.0;
    }
    
    // clean up
    for( uint32_t j=0; j<nb_loaded_keypoints; ++j )
    {
      if( descriptors[j] != 0 )
        delete [] descriptors[j];
      descriptors[j] = 0;
    }

    descriptors.clear();
    keypoints.clear();
    inlier.clear();
    
    // display statistics
    std::cout << std::endl << std::endl << " registered so far: " << registered << " / " << i+1 << std::endl;
    std::cout << " average time needed to compute the correspondences: registered: " << avrg_cor_computation_time_accepted << " rejected: " << avrg_cor_computation_time_rejected << std::endl;
    std::cout << "avrg. registration time: " << avrg_reg_time << " ( " << registered << " , avrg. inlier-ratio: " << mean_inlier_ratio_accepted << ", avrg. nb correspondences : " << mean_nb_correspondences_accepted << " ) avrg. rejection time: " << avrg_reject_time << " ( " << N_reject << ", avrg. inlier-ratio : " << mean_inlier_ratio_rejected << " avrg. nb correspondences : " << mean_nb_correspondences_rejected << " ) " << std::endl << std::endl;
  
  }
  
  ofs_details.close();
  
  ////
  // final statistics
  
  std::cout << std::endl << "#############################" << std::endl;
  std::cout << " total number registered: " << registered << " / " << nb_keyfiles << std::endl;
  std::cout << " average time for computing the vw assignments                                : " << avrg_vw_time << " s for " << avrg_nb_features << " features (on average)" << std::endl;
  std::cout << " average time for succesfully registering image (w/o time for vw assignments) : " << avrg_reg_time << " s " << std::endl;
  std::cout << " average time for rejecting an query image (w/o time for vw assignments)      : " << avrg_reject_time << " s " << std::endl;
  std::cout << " average inlier-ratio (registered)                                            : " << mean_inlier_ratio_accepted << std::endl;
  std::cout << " average inlier-ratio (rejected)                                              : " << mean_inlier_ratio_rejected << std::endl;
  std::cout << " average nb correspondences (registered)                                      : " << mean_nb_correspondences_accepted << std::endl;
  std::cout << " average nb correspondences (rejected)                                        : " << mean_nb_correspondences_rejected << std::endl;
  std::cout << " average nb features        (registered)                                      : " << mean_nb_features_accepted << std::endl;
  std::cout << " average nb features        (rejected)                                        : " << mean_nb_features_rejected << std::endl;
  std::cout << " avrg. time to compute the correspondences (registered)                       : " << avrg_cor_computation_time_accepted << std::endl;
  std::cout << " avrg. time to compute the correspondences (rejected)                         : " << avrg_cor_computation_time_rejected << std::endl;
  std::cout << " avrg. time for RANSAC (registered)                                           : " << avrg_RANSAC_time_registered << std::endl;
  std::cout << " avrg. time for RANSAC (rejected)                                             : " << avrg_RANSAC_time_rejected << std::endl;
  std::cout << " minimum inlier-ratio for RANSAC                                              : " << min_inlier << std::endl;
  std::cout << " stop after n correspondences                                                 : " << max_cor_early_term << std::endl;
  std::cout << " model consists of                                                            : " << nb_descriptors << " ";
  if (mode == 0 || mode == 2)
    std::cout << "unsigned char descriptors " << std::endl;
  else
    std::cout << "float descriptors " << std::endl;
  
  std::cout << "#############################" << std::endl;
  
 
  return 0;
}

