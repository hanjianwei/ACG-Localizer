/*===========================================================================*\
 *                                                                           *
 *                            ACG Localizer                                  *
 *      Copyright (C) 2012 by Computer Graphics Group, RWTH Aachen           *
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
#include <list>
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
#include <time.h>
#include <stdlib.h>
#include <queue>

// includes for classes dealing with SIFT-features
#include "features/SIFT_loader.hh"
#include "features/visual_words_handler.hh"

// stopwatch
#include "timer.hh"

// math functionality
#include "math/projmatrix.hh"
#include "math/matrix3x3.hh"

// tools to parse a bundler reconstruction
#include "sfm/parse_bundler.hh"
#include "sfm/bundler_camera.hh"

// RANSAC
#include "RANSAC.hh"

// exif reader to get the width and height out 
// of the exif tags of an image
#include "exif_reader/exif_reader.hh"

// ANN Libary, used to perform search in 3D
#include <ANN/ANN.h>

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
// structures
////  

// structure to indicate a feature / point correspondence search is needed
// this is used for (feature, visual word) pairs as well as for 3D points that should be matched 
// against the 2D features in the image
struct match_struct 
{
  // id of the 2D feature, the 3D point
  uint32_t feature_id;
  
  // cost of the matching, i.e., the number of descriptors in a vw (2D-to-3D)
  // or the number of SIFT distances that have to be calculated for 3D-to-2D matching
  uint32_t matching_cost;
  
  // 2D-to-3D matching (true) or 3D-to-2D matching (false)
  bool matching_type;
  
  match_struct() : feature_id(0), matching_cost(0), matching_type(true) {}
  
  match_struct( uint32_t f, uint32_t c, bool t ) : feature_id(f), matching_cost(c), matching_type(t) {} 
  
  match_struct( const match_struct &other ) : feature_id(other.feature_id), matching_cost(other.matching_cost), matching_type(other.matching_type) {}
};

////
// functions
////

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

// generic comparison function, using < to compare the second entry of two pairs
template <typename first_type, typename second_type >
inline bool cmp_second_entry_less( const std::pair< first_type, second_type >& a, const std::pair< first_type, second_type >& b )
{
  return (a.second < b.second);
}


// comparison function for the new prioritization function
inline bool cmp_priorities( const match_struct &a, const match_struct &b )
{
  return (a.matching_cost < b.matching_cost);
}


// check whether two sets have a common element or not
// sets are assumed to be stored in ascending order
// run time is in O( m + n ), where n and m are the sizes of the two sets
bool set_intersection_test( const std::set< uint32_t > &a, const std::set< uint32_t > &b )
{
  std::set< uint32_t >::const_iterator it_a = a.begin();
  std::set< uint32_t >::const_iterator it_b = b.begin();
  
  std::set< uint32_t >::const_iterator it_a_end = a.end();
  std::set< uint32_t >::const_iterator it_b_end = b.end();
  
  while( ( it_a != it_a_end ) && ( it_b != it_b_end ) )
  {
    if( *it_a < *it_b )
      ++it_a;
    else if( *it_b < *it_a )
      ++it_b;
    else // equal, that is the set has a common element
      return true;
  }
  
  return false;
}

////
// constants
////

// minimal number of inliers required to accept an image as registered
uint32_t minimal_RANSAC_solution = 12;

// SIFT-ratio value for 2D-to-3D matching. Since we store squared distances, we need the squared value 0.7^2 = 0.49
float nn_ratio_2D_to_3D = 0.49f; 

// SIFT-ratio value for 2D-to-3D matching. Since we store squared distances, we need the squared value 0.6^2 = 0.36
float nn_ratio_3D_to_2D = 0.36f; 

// the assumed minimal inlier ratio of RANSAC
float min_inlier = 0.2f;

// stop RANSAC if 60 seconds have passed
double ransac_max_time = 60.0; 

// the number of nearest neighbors to search for in 3D
int N_3D = 200;

// specify whether to use the RANSAC pre-filter (1) or not (0)
int ransac_filter = 1;

// compute a set cover from the images or use the original images
bool use_image_set_cover = true;

// the number of closest cameras to consider for clustering cameras
uint32_t consider_K_nearest_cams = 10;

// filter 3D points for visibiltiy before 3D-to-2D matching
bool filter_points = false;

// The number of correspondences to find before the search is terminated
size_t max_cor_early_term = 100;

//---------------------------------------------------------------------------------------------------------------------------------------------------------------

int main (int argc, char **argv)
{
  
  if( argc < 12 )
  {
    std::cout << "____________________________________________________________________________________________________________________________" << std::endl;
    std::cout << " -                                                                                                                        - " << std::endl;
    std::cout << " -        Localization method. Implementation of the localization framework proposed in the ECCV 2011 paper               - " << std::endl;
    std::cout << " -          T. Sattler, B. Leibe, L. Kobbelt. Improving Image-Based Localization by Active Correspondence Search.         - " << std::endl;
    std::cout << " -                               2012 by Torsten Sattler (tsattler@cs.rwth-aachen.de)                                     - " << std::endl;
    std::cout << " -                                                                                                                        - " << std::endl;
    std::cout << " - usage: acg_localizer_active_search list bundle_file nb_cluster clusters descriptors prioritization_strategy results    - " << std::endl;
    std::cout << " -                        N_3D ransac_pre_filter filter_points image_set_cover nb_cams_set_cover                          - " << std::endl;
    std::cout << " - Parameters:                                                                                                            - " << std::endl;
    std::cout << " -  list                                                                                                                  - " << std::endl;
    std::cout << " -     List containing the filenames of all the .key files that should be used as query. It is assumed that the           - " << std::endl;
    std::cout << " -     corresponding images have the same filename except of ending in .jpg.                                              - " << std::endl;
    std::cout << " -                                                                                                                        - " << std::endl;
    std::cout << " -  bundle_file                                                                                                           - " << std::endl;
    std::cout << " -     The bundle.out file generated by Bundler.                                                                          - " << std::endl;
    std::cout << " -                                                                                                                        - " << std::endl;
    std::cout << " -  nb_cluster                                                                                                            - " << std::endl;
    std::cout << " -     The number of clusters in the file containing the cluster centers. The number of cluster must be the same          - " << std::endl;
    std::cout << " -     as for the assignments computed by compute_desc_assignments.                                                       - " << std::endl;
    std::cout << " -                                                                                                                        - " << std::endl;
    std::cout << " -  clusters                                                                                                              - " << std::endl;
    std::cout << " -     The cluster centers (visual words), stored in a textfile consisting of nb_clusters * 128 floating point values.    - " << std::endl;
    std::cout << " -                                                                                                                        - " << std::endl;
    std::cout << " -  descriptors                                                                                                           - " << std::endl;
    std::cout << " -     The assignments assigning descriptors (and 3D points) to visual words. It is assumed that all descriptors          - " << std::endl;
    std::cout << " -     are stored as unsigned char values and that every visual word contains at most one descriptor for every 3D point,  - " << std::endl;
    std::cout << " -     i.e., the integer mean strategy is used to represent 3D points.                                                    - " << std::endl;
    std::cout << " -                                                                                                                        - " << std::endl;
    std::cout << " -  prioritization_strategy                                                                                               - " << std::endl;
    std::cout << " -     The strategy to combine 2D-to-3D and 3D-to-2D matching (see the paper for details):                                - " << std::endl;
    std::cout << " -      0 - Combined                                                                                                      - " << std::endl;
    std::cout << " -      1 - Direct                                                                                                        - " << std::endl;
    std::cout << " -      2 - Afterwards                                                                                                    - " << std::endl;
    std::cout << " -                                                                                                                        - " << std::endl;
    std::cout << " -  results                                                                                                               - " << std::endl;
    std::cout << " -     The program will write the results of the localization into a text file of name \"results\". It has the following  - " << std::endl;
    std::cout << " -     format, where every line in the file belongs to one query image and has the format                                 - " << std::endl;
    std::cout << " -       #inliers #(correspondences found) (time needed to compute the visual words, in seconds) (time needed for linear  - " << std::endl;
    std::cout << " -       search, in seconds) (time needed for RANSAC, in seconds) (total time needed, in seconds)                         - " << std::endl;
    std::cout << " -                                                                                                                        - " << std::endl;
    std::cout << " -  N_3D                                                                                                                  - " << std::endl;
    std::cout << " -     Specifies the number of nearest neighbors in 3D that are used as candidates for 3D-to-2D matching                  - " << std::endl;
    std::cout << " -                                                                                                                        - " << std::endl;
    std::cout << " -  ransac_pre_filter                                                                                                     - " << std::endl;
    std::cout << " -     Set to 0 to use not filter, set to 1 to use the RANSAC Pre-Filter described in the paper.                          - " << std::endl;
    std::cout << " -                                                                                                                        - " << std::endl;
    std::cout << " -  filter_points                                                                                                         - " << std::endl;
    std::cout << " -     Set to 1 to use the images in the reconstruction (or the set cover, see below) to filter the 3D points obtained    - " << std::endl;
    std::cout << " -     by nn-search in 3D: A point is accepted if it is visible with the 3D point triggering the search in at least one   - " << std::endl;
    std::cout << " -     image. Set to 0 to disable the filtering.                                                                          - " << std::endl;
    std::cout << " -                                                                                                                        - " << std::endl;
    std::cout << " -  image_set_cover                                                                                                       - " << std::endl;
    std::cout << " -     Set to 0 to use the original images from the reconstruction or to 1 to first compute a set cover that for every    - " << std::endl;
    std::cout << " -     image computes the nearest K images (in 3D) and selects those cameras for which the viewing direction differs      - " << std::endl;
    std::cout << " -     by less than 60° to the current camera to belong to the same set of cameras. Afterwards we compute a set cover     - " << std::endl;
    std::cout << " -     for those sets. You can specify the value for K with the parameter nb_cams_set_cover (default: 10).                - " << std::endl;
    std::cout << " -                                                                                                                        - " << std::endl;
    std::cout << "____________________________________________________________________________________________________________________________" << std::endl;
    return 1;
  }
  
 
  // rocketcluster
  std::string keylist( argv[1] );
  std::string bundle_file( argv[2] );
  uint32_t nb_clusters = (uint32_t) atoi( argv[3] );
  std::string cluster_file( argv[4] );
  std::string vw_assignments( argv[5] );
  int prioritization_strategy = atoi( argv[6] );
  if( prioritization_strategy != 0 && prioritization_strategy != 1 && prioritization_strategy != 2 )
  {
    std::cerr << " ERROR: Unknown prioritization strategy " << prioritization_strategy << ", aborting" << std::endl;
    return -1;
  }
  std::string outfile(argv[7] );
  
  std::cout << " Assumed minimal inlier-ratio: " << min_inlier << std::endl;
  std::cout << " Early termination after finding " << max_cor_early_term << " correspondences " << std::endl;
    
  
  N_3D = (uint32_t) atoi( argv[8] );
  std::cout << " Query expansion includes the next " << N_3D << " features " << std::endl;
  
  ransac_filter = atoi( argv[9] );
  std::cout << " RANSAC prefilter: " << ransac_filter << std::endl;
  
  filter_points = (bool) atoi( argv[10] );
  if( filter_points )
    std::cout << " Filtering 3D points before 3D-to-2D matching " << std::endl;
  else
    std::cout << " NOT Filtering 3D points before 3D-to-2D matching " << std::endl;
  
  use_image_set_cover = (bool) atoi( argv[11] );

  if( argc >= 13 )
    consider_K_nearest_cams = (uint32_t) atoi( argv[12] );
  
  if( use_image_set_cover )
    std::cout << " Using set covers instead of normal images with " << consider_K_nearest_cams << " nearest cams" << std::endl;
  else
    std::cout << " Using the original images " << std::endl;
  
  ////
  // load the visual words and their tree 
  visual_words_handler vw_handler;
  vw_handler.set_nb_trees( 1 );
  vw_handler.set_nb_visual_words( nb_clusters );
  vw_handler.set_branching( 10 );
  
  // store for every visual word at the finest level the id of its parents at levels 2 and 3
  uint32_t *parents_at_level_2 = new uint32_t[ nb_clusters ];
  uint32_t *parents_at_level_3 = new uint32_t[ nb_clusters ];
  
  for( uint32_t i=0; i<nb_clusters; ++i )
  {
    parents_at_level_2[i] = parents_at_level_3[i] = nb_clusters;
  }
  
  vw_handler.set_method(std::string("flann"));
  vw_handler.set_flann_type(std::string("hkmeans"));
    
  if( !vw_handler.create_flann_search_index( cluster_file ) )
  {
    std::cout << " ERROR: Could not load the cluster centers from " << cluster_file << std::endl;;
    return -1;
  }
    
  {
    // get the ids of the parents that we need
    int *parent_ids = new int[nb_clusters];
    for( uint32_t i=0; i<nb_clusters; ++i )
      parent_ids[i] = -1;
    
    vw_handler.get_parents_at_level_L( 2, parent_ids );
    
    // copy that to the real data structures
    uint32_t counter = 0;
    for( uint32_t i=0; i<nb_clusters; ++i )
    {
      if( parent_ids[i] == -1 )
      {
        std::cout << " uhoh, this should not happen!" << std::endl;
        ++counter;
      }
      else
      {
        parents_at_level_2[i] = (uint32_t) parent_ids[i];
        if( parents_at_level_2[i] < 0 || parents_at_level_2[i] >= 100 )
          std::cout << " WARNING: OUT OF RANGE: " << parents_at_level_2[i] << " from " << parent_ids[i] << std::endl;
      }
    }
    if( counter == nb_clusters )
    {
      std::cerr << " Some ERROR in getting the parents from level 2, stopping here " << std::endl;
      return -1;
    }
    
    for( uint32_t i=0; i<nb_clusters; ++i )
      parent_ids[i] = -1;
    
    vw_handler.get_parents_at_level_L( 3, parent_ids );
    
    // copy that
    counter = 0;
    for( uint32_t i=0; i<nb_clusters; ++i )
    {
      if( parent_ids[i] == -1 )
      {
        std::cout << " uhoh, this should not happen!" << std::endl;
        ++counter;
      }
      else
      {
        parents_at_level_3[i] = (uint32_t) parent_ids[i];
        if( parents_at_level_3[i] < 0 || parents_at_level_3[i] >= 1000 )
          std::cout << " WARNING: OUT OF RANGE: " << parents_at_level_3[i] << " from " << parent_ids[i] << std::endl;
      }
    }
    if( counter == nb_clusters )
    {
      std::cerr << " Some ERROR in getting the parents from level 3, stopping here " << std::endl;
      return -1;
    }
  }
 
  std::cout << "  done " << std::endl;
    
  
  ////
  // load the assignments of the 3D points to the visual words
  
  std::cout << "* Loading and parsing the assignments ... " << std::endl;
  
  ANNcoord** points3D = 0;
  
  // store for every visual word a list of (point id, descriptor id) pairs
  std::vector< std::vector< std::pair< uint32_t, uint32_t > > > vw_points_descriptors(nb_clusters);
  
  // store all descriptors in a single vector
  std::vector< unsigned char > all_descriptors;
  
  // for visual word, remember the number of points assigned to it
  std::vector< uint32_t > nb_points_per_vw(nb_clusters,0);
  
  
  uint32_t nb_non_empty_vw, nb_3D_points, nb_descriptors;
  
  // for every point, store its descriptor ids and the ids of the visual words this descriptors belong to
  std::vector< std::vector< uint32_t > > desc_per_point;
  std::vector< std::vector< uint32_t > > vws_per_point;
  
  for( uint32_t i=0; i<nb_clusters; ++i )
    vw_points_descriptors[i].clear();
  
  {
    std::ifstream ifs( vw_assignments.c_str(), std::ios::in | std::ios::binary  );
    
    if ( !ifs )
    {
      std::cerr << " ERROR: Cannot read the visual word assignments " << vw_assignments << std::endl;;
      return 1;
    }
    
    uint32_t nb_clusts;
    ifs.read(( char* ) &nb_3D_points, sizeof( uint32_t ) );
    ifs.read(( char* ) &nb_clusts, sizeof( uint32_t ) );
    ifs.read(( char* ) &nb_non_empty_vw, sizeof( uint32_t ) );
    ifs.read(( char* ) &nb_descriptors, sizeof( uint32_t ) );
    if( nb_clusts != nb_clusters )
    {
      std::cerr << " WARNING: Number of clusters differs! " << nb_clusts << " " << nb_clusters << std::endl;
    }
    
    std::cout << "  Number of non-empty clusters: " << nb_non_empty_vw << " number of points : " << nb_3D_points << " number of descriptors: " << nb_descriptors << std::endl;
    
    // read the 3D points and their visibility polygons
    points3D = new ANNcoord*[nb_3D_points];
    all_descriptors.resize(128*nb_descriptors);
    
    
    // load the points
    float *point_data = new float[3];
    for( uint32_t i=0; i<nb_3D_points; ++i )
    {
      points3D[i] = new ANNcoord[3];
      ifs.read(( char* ) point_data, 3 * sizeof( float ) );
      for( int j=0; j<3; ++j )
        points3D[i][j] = (ANNcoord) point_data[j];
    }
    delete [] point_data;
     
    desc_per_point.resize( nb_3D_points );
    vws_per_point.resize( nb_3D_points );
    for( uint32_t i=0; i<nb_3D_points; ++i )
    {
      desc_per_point[i].clear();
      vws_per_point[i].clear();
    }
    
    // load the descriptors
    int tmp_int;
    for( uint32_t i=0; i<nb_descriptors; ++i )
    {
      for( uint32_t j=0; j<128; ++j )
        ifs.read(( char* ) &all_descriptors[128*i+j], sizeof( unsigned char ) );
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
        
        // for every 3D point, remember the indices of its descriptors
        desc_per_point[ vw_points_descriptors[id][j].first ].push_back( vw_points_descriptors[id][j].second );
        vws_per_point[ vw_points_descriptors[id][j].first ].push_back( id );
      }
    }
    
    ifs.close();
    
    
    
    std::cout << "  done loading and parsing the assignments " << std::endl;
  }
  
  ////
  // load information about the 3D points: 
  // the connected component it belongs to and the ids of the images it is observed in
  // we obtain all these information by parsing a Bundler file
   
  // for every 3D point, store the id of its connected component and the ids of the images that see the point 
  std::vector< uint32_t > connected_component_id_per_point( nb_3D_points, 0 );
  std::vector< std::set< uint32_t > > images_per_point( nb_3D_points );
  
  for( uint32_t i=0; i<nb_3D_points; ++i )
    images_per_point[i].clear();
  
  {
    // parse the reconstruction
    parse_bundler parser;
    if( !parser.parse_data( bundle_file.c_str(), 0 ) )
    {
      std::cerr << " ERROR: Could not parse the bundler file " << bundle_file << std::endl;
      return -1;
    }
    uint32_t nb_cameras = parser.get_number_of_cameras();
    uint32_t nb_points_bundler = parser.get_number_of_points();
    std::vector< feature_3D_info >& feature_infos = parser.get_feature_infos();
    
    if( nb_points_bundler != nb_3D_points )
    {
      std::cerr << " ERROR: The number of points in the binary file ( " << nb_3D_points << " ) and in the reconstruction ( " << nb_points_bundler << " ) differ!" << std::endl;
      return -1;
    }
    
    // compute the connected components
    std::cout << "  Computing connected components " << std::endl;
    // for every camera, get the corresponding keypoint ids
    std::vector< std::vector< uint32_t > > cam_keys( nb_cameras );
    for( uint32_t i=0; i<nb_cameras; ++i )
      cam_keys[i].clear();


    for( uint32_t i=0; i<nb_points_bundler; ++i )
    {
      for( size_t j=0; j<feature_infos[i].view_list.size(); ++j )
        cam_keys[ feature_infos[i].view_list[j].camera ].push_back( i );
    }
      
    std::vector< int > cam_ccs( nb_cameras, -1 );
    std::vector< int > point_ccs( nb_points_bundler, -1 );
    int nb_ccs = 0;
    
    std::queue< uint32_t > remaining_points;
    
    for( uint32_t i=0; i<nb_points_bundler; ++i )
      remaining_points.push(i);
      
    while( !remaining_points.empty() )
    {
      uint32_t cur_point = remaining_points.front();
      
      remaining_points.pop();
      
      if( point_ccs[ cur_point] != -1 )
        continue;
      
      std::queue< uint32_t > recursive_points;
        
      // create a new connected component
      point_ccs[cur_point] = nb_ccs;
      ++nb_ccs;
      
      for( size_t j=0; j<feature_infos[cur_point].view_list.size(); ++j )
      {
        uint32_t cam_id = feature_infos[cur_point].view_list[j].camera;
        if( cam_ccs[cam_id] != -1 && cam_ccs[cam_id] != point_ccs[cur_point] )
          std::cout << " ERROR: ambigous ids for camera! " <<std::endl;
        
        if( cam_ccs[ cam_id ] == -1 )
        {
          cam_ccs[ cam_id ] = point_ccs[ cur_point ];
          for( size_t k=0; k<cam_keys[ cam_id ].size(); ++k )
          {
            if( point_ccs[ cam_keys[ cam_id ][k] ] != -1 && point_ccs[ cam_keys[ cam_id ][k] ] != point_ccs[ cur_point ] )
          std::cout << " ERROR: ambigous ids for point " <<std::endl;
            
            if( point_ccs[ cam_keys[ cam_id ][k] ] == -1 )
          recursive_points.push( cam_keys[ cam_id ][k] );
          }
        }
      }
      
      while( !recursive_points.empty() )
      {
        uint32_t c = recursive_points.front();
        recursive_points.pop();
        
        if( point_ccs[c] != -1 && point_ccs[c] != point_ccs[ cur_point ] )
          std::cout << " ERROR: ambigous ids for point " <<std::endl;
        
        if( point_ccs[c] == -1 )
        {
          // create a new connected component
          point_ccs[c] = point_ccs[cur_point];
          
          for( size_t j=0; j<feature_infos[c].view_list.size(); ++j )
          {
            uint32_t cam_id = feature_infos[c].view_list[j].camera;
            if( cam_ccs[cam_id] != -1 && cam_ccs[cam_id] != point_ccs[cur_point] )
              std::cout << " ERROR: ambigous ids for camera! " <<std::endl;
            
            if( cam_ccs[ cam_id ] == -1 )
            {
              cam_ccs[ cam_id ] = point_ccs[ cur_point ];
              for( size_t k=0; k<cam_keys[ cam_id ].size(); ++k )
              {
                if( point_ccs[ cam_keys[ cam_id ][k] ] != -1 && point_ccs[ cam_keys[ cam_id ][k] ] != point_ccs[ cur_point ] )
                  std::cout << " ERROR: ambigous ids for point " <<std::endl;
                
                if( point_ccs[ cam_keys[ cam_id ][k] ] == -1 )
                  recursive_points.push( cam_keys[ cam_id ][k] );
              }
            }
          }
        }
      }
    }
    std::cout << "   Found " << nb_ccs << " connected components " << std::endl;
      
    // check if all non-empty images and all points belong to a connected component
    for( uint32_t i=0; i<nb_points_bundler; ++i )
    {
      if( point_ccs[i] == -1 || point_ccs[i] >= nb_ccs )
      {
        std::cout << " ERROR : Point " << i << " has connected component id " << point_ccs[i] << std::endl;
        return -1;
      }
    }
      
    for( uint32_t i=0; i<nb_cameras; ++i )
    {
      if( (cam_ccs[i] == -1 || cam_ccs[i] >= nb_ccs) && (cam_keys[i].size() > 0 ) )
      {
        std::cout << " ERROR : Non-empty camera " << i << " has connected component id " << cam_ccs[i] << std::endl;
        return -1;
      }
      
      // check if any of the points visible in this camera has a different conncected component id
      for( size_t k=0; k<cam_keys[ i ].size(); ++k )
      {
        if( point_ccs[ cam_keys[i][k] ] != cam_ccs[i] )
        {
          std::cerr << " ERROR : Found point " << cam_keys[i][k] << " in camera " << i << " that belongs to component " << point_ccs[ cam_keys[i][k] ] << " while the camera belongs to component " << cam_ccs[i] << std::endl;
          return -1;
        }
      }
      
      if( use_image_set_cover )
        std::sort( cam_keys[i].begin(), cam_keys[i].end() );
      else
        cam_keys[i].clear();
    }
      
    ////
    // if we want to represent the set of images with a smaller set, we now compute it
    
    // recall for every image by which other images it is covered
    std::vector< std::set< uint32_t > > image_covered_by;
    image_covered_by.clear();
    std::vector< std::set< uint32_t > > images_covered_by_image;
    images_covered_by_image.clear(); 
    
    if( use_image_set_cover )
    {
      image_covered_by.resize( nb_cameras );
      
      images_covered_by_image.resize( nb_cameras );
      
      std::cout << "  Computing the set cover for all images, each image covers itself and (at most) the " << consider_K_nearest_cams << " images that have the largest number of 3D points in common with it " << std::endl;
      
      for( uint32_t i=0; i<nb_cameras; ++i )
      {
        images_covered_by_image[i].clear();
        images_covered_by_image[i].insert(i);
      }
      
      ////
      // now add cameras that are close in 3D and have a similar viewing direction (angle between directions beneath 60°)
      std::vector< OpenMesh::Vec3f > camera_positions( nb_cameras );
      std::vector< OpenMesh::Vec3f > cam_viewing_dirs( nb_cameras );
      
      std::vector< bundler_camera > &bundle_cams = parser.get_cameras();
      
      for( uint32_t i=0; i<nb_cameras; ++i )
      {
        camera_positions[i] = bundle_cams[i].get_cam_position_f();
        cam_viewing_dirs[i] = bundle_cams[i].get_cam_global_vec_f( OpenMesh::Vec3d( 0.0, 0.0, -1.0 ) );
      }
      
      for( uint32_t i=0; i<nb_cameras; ++i )
      {
        std::vector< std::pair< uint32_t, float > > cameras_distances( nb_cameras );
        
        for( uint32_t j=0; j<nb_cameras; ++j )
        {
          cameras_distances[j].first = j;
          cameras_distances[j].second = ( camera_positions[i] - camera_positions[j] ).length();
        }
        
        std::sort( cameras_distances.begin(), cameras_distances.end(), cmp_second_entry_less< uint32_t, float > );
        
        // now pick the 10 nearest cameras looking in a similar direction out of the nearest 20
        // cameras from the same connected component
        uint32_t counter = 0;
        uint32_t found = 0;
        for( uint32_t j=0; j<nb_cameras && counter<consider_K_nearest_cams /*&& found<5*/; ++j )
        {
          // take care that the camera we are looking at is not camera i but is in the 
          // same connected component!
          if( cameras_distances[j].first == i || cam_ccs[i] != cam_ccs[ cameras_distances[j].first ] )
            continue;
          
          // remember we only want to look at the 20 nearest cameras and pick at most 10 from them
          ++counter;
          
          if( ( cam_viewing_dirs[i] | cam_viewing_dirs[ cameras_distances[j].first ] ) >= 0.5f )
          {
            images_covered_by_image[i].insert( cameras_distances[j].first );
            ++found;
          }
        }
        
      }
      
      camera_positions.clear();
      cam_viewing_dirs.clear();
      
      ////
      // now we compute the set cover
      
      // for every image, track how many new images it can cover
      std::vector< std::pair< uint32_t, uint32_t > > nb_new_images_covered( nb_cameras );
      
      for( uint32_t i=0; i<nb_cameras; ++i )
      {
        image_covered_by[i].clear();
        nb_new_images_covered[i].first = i;
        nb_new_images_covered[i].second = images_covered_by_image[i].size();
      }
      
      // map image ids to a smaller range
      std::vector< int > new_image_ids( nb_cameras, -1 );
      int size_set_cover = 0;
      
      while( !nb_new_images_covered.empty() )
      {
        std::sort( nb_new_images_covered.begin(), nb_new_images_covered.end(), cmp_second_entry_less< uint32_t, uint32_t > );
        if( nb_new_images_covered.back().second == 0 )
          break;
        
        uint32_t cam_id_ = nb_new_images_covered.back().first;
        
        // add new image to set cover
        new_image_ids[ cam_id_ ] = size_set_cover;
        
        // markt its images as covered
        for( std::set< uint32_t >::const_iterator it = images_covered_by_image[ cam_id_ ].begin(); it != images_covered_by_image[ cam_id_ ].end(); ++it )
          image_covered_by[ *it ].insert( (uint32_t) size_set_cover );
        
        ++size_set_cover;
        
        // pop first element
        nb_new_images_covered.pop_back();
        
        // recompute the nb of new images each image can cover
        for( std::vector< std::pair< uint32_t, uint32_t > >::iterator it = nb_new_images_covered.begin(); it != nb_new_images_covered.end(); ++it )
        {
          it->second = 0;
          for( std::set< uint32_t >::const_iterator it2 = images_covered_by_image[ it->first ].begin(); it2 != images_covered_by_image[ it->first ].end(); ++it2 )
          {
            if( image_covered_by[ *it2 ].empty() )
          it->second += 1;
          }
        }
      }
      
      nb_new_images_covered.clear();
      
      std::cout << "   Set cover contains " << size_set_cover << " cameras out of " << nb_cameras << std::endl;
    }
    else
    {
      image_covered_by.clear();
      images_covered_by_image.clear(); 
    }
      
    // copy the information
    for( uint32_t i=0; i<nb_points_bundler; ++i )
      connected_component_id_per_point[i] = point_ccs[i];
    
    cam_ccs.clear();
    point_ccs.clear();
    
    
    std::cout << "  Reading the images per 3D point " << std::endl;
    for( uint32_t i=0; i<nb_points_bundler; ++i )
    {
      images_per_point[i].clear();
      
      for( size_t j=0; j<feature_infos[i].view_list.size(); ++j )
      {
        if( use_image_set_cover )
        {
          uint32_t cam_id_ = feature_infos[i].view_list[j].camera;
          for( std::set< uint32_t >::const_iterator it = image_covered_by[ cam_id_ ].begin(); it != image_covered_by[ cam_id_ ].end(); ++it )
            images_per_point[i].insert( *it );
        }
        else
          images_per_point[i].insert( feature_infos[i].view_list[j].camera );
      }
      
      if( images_per_point[i].size() == 0 )
        std::cout << " WARNING: Point " << i << " is visible in no image!" << std::endl;
    }
    
    // clean up
    if( use_image_set_cover )
    {
      for( uint32_t i=0; i<nb_cameras; ++i )
      {
        image_covered_by[i].clear();
        images_covered_by_image[i].clear();
      }
      image_covered_by.clear();
      images_covered_by_image.clear();
    }
      
    parser.clear();
    feature_infos.clear();
  }
  
  
  ////
  // create the kd-trees for the 3D points to enable search in 3D, one for each connected component
  
  // get the number of connected components
  uint32_t nb_connected_components = *std::max_element( connected_component_id_per_point.begin(), connected_component_id_per_point.end() );
  nb_connected_components += 1;
  
  std::cout << " * creating kd-trees for 3D points, one for each of the " << nb_connected_components << " connected components " << std::endl;
  
  // for every connected component, get the number of points in it
  std::vector< uint32_t > nb_points_per_component( nb_connected_components, 0 );
  for( std::vector< uint32_t >::const_iterator it = connected_component_id_per_point.begin(); it != connected_component_id_per_point.end(); ++it )
    nb_points_per_component[ *it ] += 1;
  
  // for every connected component, get pointers to its 3D points
  float **points_per_component[ nb_connected_components ];
  uint32_t *indices_per_component[ nb_connected_components ];
  
  // store pointers to the appropriate points
  {
    std::vector< uint32_t > cc_point_counter( nb_connected_components, 0 );
    
    for( uint32_t i=0; i<nb_connected_components; ++i )
    {
      points_per_component[i] = new float*[ nb_points_per_component[i] ];
      indices_per_component[i] = new uint32_t[ nb_points_per_component[i] ];
    }
    
    for( uint32_t i=0; i<nb_3D_points; ++i )
    {
      uint32_t cc_id = connected_component_id_per_point[i];
      points_per_component[ cc_id ][ cc_point_counter[ cc_id ] ] = points3D[i];
      indices_per_component[ cc_id ][ cc_point_counter[ cc_id ] ] = i;
      cc_point_counter[ cc_id ] += 1;
    }
  }
  
  // create the trees
  ANNkd_tree *kd_trees[ nb_connected_components ];
  
  for( uint32_t i=0; i<nb_connected_components; ++i )
  {
    int nb_points_3D_int = (int) nb_points_per_component[i];
    kd_trees[i] = new ANNkd_tree( points_per_component[i], nb_points_3D_int, 3 );
  }
  
  int nb_points_3D_int = (int) nb_3D_points;
  ANNkd_tree *kd_tree = new ANNkd_tree( points3D, nb_points_3D_int, 3 );
  
  // and the search structures
  annMaxPtsVisit( 0 );
  
  ANNidxArray indices = new ANNidx[ N_3D ];
  ANNdistArray distances = new ANNdist[ N_3D ];
  
  std::cout << "  done " << std::endl;
 
  
  ////
  // now load the filenames of the query images
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
      {
        key_filenames.push_back(tmp_string);
      }
    }
    ifs.close();
    std::cout << " done loading " << key_filenames.size() << " keyfile names " << std::endl;
  }
  
  uint32_t nb_keyfiles = key_filenames.size();
  
  ////
  // do the actual localization
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
  double mean_assignment_time_QE = 0.0;
  double mean_QE_time = 0.0;
  
  
  
  
  std::ofstream ofs( outfile.c_str(), std::ios::out );
    
  if( !ofs.is_open() )
  {
    std::cerr << " Could not write results to " << outfile << std::endl;
    return 1;
  }
  
  
  uint32_t registered = 0;
  
  std::vector< uint32_t > computed_visual_words( 50000, 0 );
  std::vector< uint32_t > computed_visual_words_low_dim( 50000,0 );

  
  ////
  // for every visual word on level 2 or 3 (so at max 1000 visual words)
  // store the ids of the 2D features that are mapped to that visual word
  std::vector< std::vector< uint32_t > > features_per_vw( 1000 );
  for( int i=0; i<1000; ++i )
    features_per_vw[i].resize(100);
  
  
  for( uint32_t i=0; i<nb_keyfiles; ++i, N+=1.0 )
  {
    std::cout << std::endl << " --------- " << i+1 << " / " << nb_keyfiles << " --------- " << std::endl;
    
    SIFT_loader key_loader;
    std::cout << key_filenames[i] << std::endl;
    key_loader.load_features( key_filenames[i].c_str(), LOWE );
    
    std::vector< unsigned char* >& descriptors = key_loader.get_descriptors();
    std::vector< SIFT_keypoint >& keypoints = key_loader.get_keypoints();
    
    uint32_t nb_loaded_keypoints = (uint32_t) keypoints.size();
    
   
    // center the keypoints around the center of the image
    // first we need to get the dimensions of the image
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
    
    std::cout << " loaded " << nb_loaded_keypoints << " descriptors" << std::endl;
    
    // assign the descriptors to the visual words   
    
    Timer timer;
    timer.Init();
    timer.Start();
    
   
    if( computed_visual_words.size() < nb_loaded_keypoints )
    {
      computed_visual_words.resize( nb_loaded_keypoints );
      computed_visual_words_low_dim.resize( nb_loaded_keypoints );
    }
    
    
    std::set< size_t > unique_vw;
    unique_vw.clear();
    
    vw_handler.set_nb_paths( 1 );
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
    
    // first, compute for every feature in the image its visual word in the coarser vocabulary
    // compute the lower dimensions
    int low_dim_choosen = 0;
    int nb_low_dim = 100;
    
    if( nb_loaded_keypoints > 5000 )
    {
      low_dim_choosen = 1;
      nb_low_dim = 1000;
    }
    
    for( int j=0; j<nb_low_dim; ++j )
      features_per_vw[j].clear();
    
    
      
    if( low_dim_choosen == 0 )
    {
      for( size_t j=0; j<nb_loaded_keypoints; ++j )
        features_per_vw[ parents_at_level_2[ computed_visual_words[j] ] ].push_back( j );
    }
    else
    {
      for( size_t j=0; j<nb_loaded_keypoints; ++j )
        features_per_vw[ parents_at_level_3[ computed_visual_words[j] ] ].push_back( j );
    }
      
   
    // store the correspondences for RANSAC
    std::vector< float > c2D, c3D;
    c2D.clear();
    c3D.clear();
    
    std::vector< std::pair< uint32_t, uint32_t > > final_correspondences; // first the 2D, then the 3D point
    final_correspondences.clear();
    
    uint32_t max_corr = 0;
    uint32_t no_nn = 0;
    uint32_t no_sn_neighbor = 0;
    uint32_t failed_ratio = 0;
    
    
    // compute the priorities
    
    std::list< match_struct > priorities( nb_loaded_keypoints );
    std::list< match_struct >::iterator priorities_it = priorities.begin();
    
    for( uint32_t j=0; j<nb_loaded_keypoints; ++j, ++priorities_it )
    {
      priorities_it->feature_id = j;
      priorities_it->matching_cost = nb_points_per_vw[computed_visual_words[j]];
      priorities_it->matching_type = true;
    }
    
    // keep track which 2D features are used for correspondences
    std::vector< bool > feature_in_correspondence( nb_loaded_keypoints, false );
    
    // keep track which 2D features are used for 3D-to-2D correspondences and 
    // to which point they match. This is needed to be able to override 3D-to-2D
    // correspondences with new 3D-to-2D correspondences.
    // To signal that a feature is not part of a 3D-to-2D correspondence, we use
    // UINT32_MAX
    std::vector< uint32_t > point_per_feature( nb_loaded_keypoints, UINT32_MAX );
    
    // keep track which 3D points have been used in the query expansion
    std::set< uint32_t > used_3D_points;
    used_3D_points.clear();
    
    priorities.sort( cmp_priorities );

    
    // similarly, we store for each 3D point the corresponding 2D feature as well as the squared distance
    std::map< uint32_t, std::pair< uint32_t, int > > corr_3D_to_2D;
    corr_3D_to_2D.clear();
    
    std::map< uint32_t, nearest_neighbors >::iterator map_it_2D;
    std::map< uint32_t, nearest_neighbors_float >::iterator map_it_2D_float;
    std::map< uint32_t, nearest_neighbors_multiple >::iterator map_it_2D_multple;
    std::map< uint32_t, std::pair< uint32_t, int > >::iterator map_it_3D;

    // compute nearest neighbors using 2D-to-3D and 3D-to-2D matching

    uint32_t nb_considered_points = 0;
    uint32_t nb_considered_points_counter = 0;
    if( prioritization_strategy == 0 || prioritization_strategy == 1 )
    {
      for( priorities_it = priorities.begin(); priorities_it != priorities.end(); ++priorities_it )
      {
//         std::cout << priorities_it->feature_id << " " << priorities_it->matching_type << " " << priorities_it->matching_cost << std::endl;
        // check the matching type, and handle the different matching directions accordingly
        if( priorities_it->matching_type )
        {
          
          ////
          // 2D-to-3D matching, similar to the ICCV 2011 paper
        
          uint32_t j_index = priorities_it->feature_id;
          
          if( feature_in_correspondence[ j_index ] )
            continue;
          
          uint32_t assignment = uint32_t( computed_visual_words[j_index] );
          
          ++nb_considered_points_counter;
          
          nearest_neighbors nn;
          
          if( priorities_it->matching_cost > 0 )
          {
            // find nearest neighbor for 2D feature, update nearest neighbor information for 3D points if necessary
            size_t nb_poss_assignments = vw_points_descriptors[assignment].size();
            
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
              if( nn.get_ratio() < nn_ratio_2D_to_3D )
              {
                // we found one, so we need check for mutual nearest neighbors
                map_it_3D = corr_3D_to_2D.find( nn.nn_idx1 );
              
                if( map_it_3D != corr_3D_to_2D.end() )
                {
                  // a correspondence to the same 3D point already exists
                  // so we have to check whether we have to update it or not
                  if( map_it_3D->second.second > nn.dist1 )
                  {
                    feature_in_correspondence[ map_it_3D->second.first ] = false;
                    
                    map_it_3D->second.first = j_index;
                    map_it_3D->second.second = nn.dist1;
                    
                    feature_in_correspondence[j_index ] = true;
                  }
                }
                else
                {
                  corr_3D_to_2D.insert( std::make_pair( nn.nn_idx1, std::make_pair( j_index, nn.dist1 ) ) );
                  feature_in_correspondence[ j_index ] = true;
                  used_3D_points.insert( nn.nn_idx1 );
                  
                  // avoid query expansion if we are not going to use it anyways
                  if( corr_3D_to_2D.size() >= max_cor_early_term )
                  {
                    nb_considered_points = nb_considered_points_counter;
                    break;
                  }
                  
                  //// START ACTIVE SEARCH
                  // we only have to do the nn query in 3D space if we have 
                  // found a new correspondence, not updated an old one (should be happening seldomly anyways)
                  
                  // find the nearest neighbors in 3D
                  // ANN will throw an exception if we search for more points than contained in the connected component, so 
                  // we have to adjust the number of points we search for
                  int N3D_ = std::min( N_3D, (int) nb_points_per_component[ connected_component_id_per_point[ nn.nn_idx1 ] ] );
                  kd_trees[ connected_component_id_per_point[ nn.nn_idx1 ] ]->annkSearch( points3D[nn.nn_idx1], N3D_, indices, distances );
                  
                  ////
                  // find new matching possibilities and insert them into the correct position 
                  // in the list
                  std::list< match_struct > new_matching_possibilities;
                  new_matching_possibilities.clear();
                  

                  for( int kk=0; kk<N3D_; ++kk )
                  {
                    if( indices[kk] < 0 )
                      break;
                                      
                    uint32_t candidate_point = indices_per_component[ connected_component_id_per_point[ nn.nn_idx1 ] ][ indices[kk] ];
                    
                    // check if we have already used or visited this 3D point (or promise to visit it later on
                    if( used_3D_points.find( candidate_point ) == used_3D_points.end() )
                    {
                      // visibility filter
                      if( filter_points && ( !set_intersection_test( images_per_point[candidate_point], images_per_point[ nn.nn_idx1 ] ) ) )
                        continue;
                      
                      // promise that we will (eventually) look at this 3D point
                      used_3D_points.insert(candidate_point);
                                            
                      ////
                      // compute the matching cost of this particular 3D point
                      // therefore, we have to identify to which visual words (on a higher level)
                      // the descriptors of this point are assigned how many times
                      
                      match_struct new_match( candidate_point, 0, false );
                      
                      if( low_dim_choosen == 0 )
                      {
                        for(std::vector< uint32_t >::const_iterator it_vws = vws_per_point[ candidate_point].begin(); it_vws != vws_per_point[ candidate_point].end(); ++it_vws)
                          new_match.matching_cost += features_per_vw[ parents_at_level_2[ *it_vws ] ].size();
                      }
                      else
                      {
                        for(std::vector< uint32_t >::const_iterator it_vws = vws_per_point[ candidate_point].begin(); it_vws != vws_per_point[ candidate_point].end(); ++it_vws)
                          new_match.matching_cost += features_per_vw[ parents_at_level_3[ *it_vws ] ].size();
                      }
                      
                      // push back the new matching possibilities
                      new_matching_possibilities.push_back( new_match );
                    }
                  }
                  
                  ////
                  // insert the list of new matching possibilities into our prioritization search structure
                  // to do so efficiently, we sort the new list and use insertion sort to update the prioritization scheme
                  
                  // sort
                  new_matching_possibilities.sort( cmp_priorities );
                  
                  // depending on the prioritization_strategy we either insert all new matching possibilities directly after the next 
                  // entry in the current prioritization list (prioritization_strategy == 1) or use insertion sort to insert it in a
                  // sorted fashion (prioritization_strategy == 0)
                  
                  if( prioritization_strategy == 0 )
                  {
                    // insertion sort: Notice that we start with the NEXT entry in the current prioritization list
                    // since we are nearly done with the old one
                    std::list< match_struct >::iterator insertion_it = priorities_it;
                    ++insertion_it;
                    
                    std::list< match_struct >::iterator to_insert_it = new_matching_possibilities.begin();
                    
                    while( to_insert_it != new_matching_possibilities.end() )
                    {
                      if( insertion_it == priorities.end() )
                      {
                        priorities.insert( insertion_it, to_insert_it, new_matching_possibilities.end() );
                        break;
                      }
                      
                      if( insertion_it->matching_cost > to_insert_it->matching_cost )
                      {
                        priorities.insert( insertion_it, *to_insert_it );
                        ++to_insert_it;
                      }
                      else
                        ++insertion_it;
                    }
                  }
                  else
                  {
                    std::list< match_struct >::iterator insertion_it = priorities_it;
                    ++insertion_it;
                    
                    priorities.insert( insertion_it, new_matching_possibilities.begin(), new_matching_possibilities.end() );
                  }
                  
                  new_matching_possibilities.clear();
                  
                  //// STOP ACTIVE SEARCH
                }
              }
            }
          }
        }
        else
        {
          ////
          // 3D-to-2D correspondence, handle it accordingly
          uint32_t candidate_point = priorities_it->feature_id;
          
          // check if we have already found a correspondence for that 3D point
          // if so, we don't need to find a new one since it was found during
          // 2D-to-3D matching which we trust more
          if( corr_3D_to_2D.find( candidate_point ) != corr_3D_to_2D.end() )
            continue;
            
          
          uint32_t nb_desc_for_point = (uint32_t) desc_per_point[ candidate_point ].size();
          std::vector< uint32_t > low_dim_vw_ids( nb_desc_for_point );
          
              
          // map all descriptors of that point to the lower dimensional descriptors
          std::set< uint32_t > low_dim_vw;
          low_dim_vw.clear();
          
          uint32_t counter = 0;
          if( low_dim_choosen == 0 )
          {
            for( std::vector< uint32_t >::const_iterator it_vws = vws_per_point[ candidate_point].begin(); it_vws != vws_per_point[ candidate_point].end(); ++it_vws, ++counter )
            {
              low_dim_vw_ids[counter] = parents_at_level_2[ *it_vws ];
              low_dim_vw.insert( parents_at_level_2[ *it_vws ] );
            }
          }
          else
          {
            for( std::vector< uint32_t >::const_iterator it_vws = vws_per_point[ candidate_point].begin(); it_vws != vws_per_point[ candidate_point].end(); ++it_vws, ++counter )
            {
              low_dim_vw_ids[counter] = parents_at_level_3[ *it_vws ];
              low_dim_vw.insert( parents_at_level_3[ *it_vws ] );
            }
          }

          
          // try to find a new correspondence
          nearest_neighbors_multiple nn_exp;
          uint64_t descriptor_index(0);
          for( std::set< uint32_t >::const_iterator activated_vw = low_dim_vw.begin(); activated_vw != low_dim_vw.end(); ++activated_vw )
          {
            counter = 0;
            for( std::vector< uint32_t >::const_iterator it_desc = desc_per_point[ candidate_point].begin(); it_desc != desc_per_point[ candidate_point].end(); ++it_desc, ++counter )
            {
              if( low_dim_vw_ids[counter] != *activated_vw )
                continue;
              
              descriptor_index = (*it_desc) * sift_dim;
              
              for( std::vector< uint32_t >::const_iterator feature_it = features_per_vw[ *activated_vw ].begin(); feature_it != features_per_vw[ * activated_vw ].end(); ++feature_it )
              {
                int dist = 0;
                int x = 0;
                for( uint64_t jj=0; jj<128; ++jj )
                {
                  x = ( (int) descriptors[*feature_it][jj] ) - ( (int) all_descriptors[ descriptor_index + jj ] );
                  dist += x*x;
                }
                
                nn_exp.update( *feature_it, dist );
              }
            }
          }
          
          if( nn_exp.dist1 >= 0 )
          {
            if( nn_exp.dist2 >= 0 )
            {
              if( nn_exp.get_ratio() < nn_ratio_3D_to_2D )
              {
                // we don't want to overwrite 2D-to-3D correspondence, only 3D-to-2D correspondence 
                if( !feature_in_correspondence[ nn_exp.nn_idx1 ] )
                {
                  // no existing correspondence
                  corr_3D_to_2D.insert( std::make_pair( candidate_point, std::make_pair( nn_exp.nn_idx1, nn_exp.dist1 ) ) );
                  feature_in_correspondence[ nn_exp.nn_idx1 ] = true;
                  point_per_feature[ nn_exp.nn_idx1 ] = candidate_point;
                }
                else if( point_per_feature[ nn_exp.nn_idx1 ] != UINT32_MAX )
                {
                  // only 3D-to-2D correspondence
                  // overwrite if the absolute SIFT distance is smaller 
                  
                  // first, we to get an iterator to the corresponding 3D point
                  map_it_3D = corr_3D_to_2D.find( point_per_feature[ nn_exp.nn_idx1 ] );
                  
                  // this has to exist, otherwise we would not have labeled it as having a correspondence
                  if( map_it_3D->second.second > nn_exp.dist1 )
                  {
                    // update the correspondence
                    
                    // which means we first have to remove the old one
                    corr_3D_to_2D.erase( map_it_3D );
                    
                    point_per_feature[ nn_exp.nn_idx1 ] = candidate_point;
                    
                    corr_3D_to_2D.insert( std::make_pair( candidate_point, std::make_pair( nn_exp.nn_idx1, nn_exp.dist1 ) ) );
                  }
                }
              }
            }
          }
        }
        
//         std::cout << " "  << corr_3D_to_2D.size() << std::endl;
        
        if( corr_3D_to_2D.size() >= max_cor_early_term )
        {
          nb_considered_points = nb_considered_points_counter;
          break;
        }
      }
    }
    else if( prioritization_strategy == 2 )
    {
      ////
      // first perform 2D-to-3D matching, then 3D-to-2D matching until enough correspondences are found
      for( priorities_it = priorities.begin(); priorities_it != priorities.end(); ++priorities_it )
      {
        ////
        // 2D-to-3D matching, similar to ICCV 2011 version
          
        uint32_t j_index = priorities_it->feature_id;
        
        // testing whether a correspondence is already part of a correspondence
        // is not necessary for pure 2D-to-3D matching, but we will later need
        // the information for 3D-to-2D matching
        
        uint32_t assignment = uint32_t( computed_visual_words[j_index] );
        
        ++nb_considered_points_counter;
        
        nearest_neighbors nn;


        if( priorities_it->matching_cost > 0 )
        {
          // find nearest neighbor for 2D feature, update nearest neighbor information for 3D points if necessary
          size_t nb_poss_assignments = vw_points_descriptors[assignment].size();
          
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
            if( nn.get_ratio() < nn_ratio_2D_to_3D )
            {
              // we found one, so we need check for mutual nearest neighbors
              map_it_3D = corr_3D_to_2D.find( nn.nn_idx1 );
        
              if( map_it_3D != corr_3D_to_2D.end() )
              {
                // a correspondence to the same 3D point already exists
                // so we have to check whether we have to update it or not
                if( map_it_3D->second.second > nn.dist1 )
                {
                  feature_in_correspondence[ map_it_3D->second.first ] = false;
                  
                  map_it_3D->second.first = j_index;
                  map_it_3D->second.second = nn.dist1;
                  
                  feature_in_correspondence[j_index ] = true;
                }
              }
              else
              {
                corr_3D_to_2D.insert( std::make_pair( nn.nn_idx1, std::make_pair( j_index, nn.dist1 ) ) );
                feature_in_correspondence[ j_index ] = true;
                used_3D_points.insert( nn.nn_idx1 );
              }
            }
          }
        }

        if( corr_3D_to_2D.size() >= max_cor_early_term )
        {
          nb_considered_points = nb_considered_points_counter;
          break;
        }
      }
      
      ////
      // check whether we have to compute correspondences from 3D-to-2D
      if( corr_3D_to_2D.size() < max_cor_early_term )
      {
        ////
        // perform the nn search in 3D to get new potential matches
        
        std::list< match_struct > point_to_image_matches;
        point_to_image_matches.clear();
        
        for( map_it_3D = corr_3D_to_2D.begin(); map_it_3D != corr_3D_to_2D.end(); ++map_it_3D )
        {
          //// START ACTIVE SEARCH
          int N3D_ = std::min( N_3D, (int) nb_points_per_component[ connected_component_id_per_point[ map_it_3D->first ] ] );
          kd_trees[ connected_component_id_per_point[ map_it_3D->first ] ]->annkSearch( points3D[map_it_3D->first], N3D_, indices, distances );
          
          ////
          // find new matching possibilities and insert them into the correct position 
          // in the list

          for( int kk=0; kk<N3D_; ++kk )
          {
            if( indices[kk] < 0 )
              break;
            
            uint32_t candidate_point = indices_per_component[ connected_component_id_per_point[ map_it_3D->first ] ][ indices[kk] ];
            
            // check if we have already used or visited this 3D point (or promise to visit it later on
            if( used_3D_points.find( candidate_point ) == used_3D_points.end() )
            {
              // visibility filter
              if( filter_points && ( !set_intersection_test( images_per_point[candidate_point], images_per_point[ map_it_3D->first ] ) ) )
                continue;
              
              // promise that we will (eventually) look at this 3D point
              used_3D_points.insert(candidate_point);
              
              ////
              // compute the matching cost of this particular 3D point
              // therefore, we have to identify to which visual words (on a higher level)
              // the descriptors of this point are assigned how many times
              
              match_struct new_match( candidate_point, 0, false );
              
              if( low_dim_choosen == 0 )
              {
                for( std::vector< uint32_t >::const_iterator it_vws = vws_per_point[ candidate_point].begin(); it_vws != vws_per_point[ candidate_point].end(); ++it_vws )
                  new_match.matching_cost += features_per_vw[ parents_at_level_2[ *it_vws ] ].size();
              }
              else
              {
                for( std::vector< uint32_t >::const_iterator it_vws = vws_per_point[ candidate_point].begin(); it_vws != vws_per_point[ candidate_point].end(); ++it_vws )
                  new_match.matching_cost += features_per_vw[ parents_at_level_3[ *it_vws ] ].size();
              }
              
              // push back the new matching possibilities
              point_to_image_matches.push_back( new_match );
            }
          }
        }
      
        // sort the list in ascending number of search cost
        point_to_image_matches.sort( cmp_priorities );
        //
        //// STOP ACTIVE SEARCH

        ////
        // do the 3D-to-2D matching
        for( priorities_it = point_to_image_matches.begin(); priorities_it != point_to_image_matches.end(); ++priorities_it )
        {
          ////
          // 3D-to-2D correspondence, handle it accordingly
          uint32_t candidate_point = priorities_it->feature_id;
          
          // check if we have already found a correspondence for that 3D point
          // if so, we don't need to find a new one since it was found during
          // 2D-to-3D matching which we trust more
          if( corr_3D_to_2D.find( candidate_point ) != corr_3D_to_2D.end() )
            continue;
            
          
          uint32_t nb_desc_for_point = (uint32_t) desc_per_point[ candidate_point ].size();
          std::vector< uint32_t > low_dim_vw_ids( nb_desc_for_point );
          
              
          // map all descriptors of that point to the lower dimensional descriptors
          std::set< uint32_t > low_dim_vw;
          low_dim_vw.clear();
          
          uint32_t counter = 0;

          if( low_dim_choosen == 0 )
          {
            for( std::vector< uint32_t >::const_iterator it_vws = vws_per_point[ candidate_point].begin(); it_vws != vws_per_point[ candidate_point].end(); ++it_vws, ++counter )
            {
              low_dim_vw_ids[counter] = parents_at_level_2[ *it_vws ];
              low_dim_vw.insert( parents_at_level_2[ *it_vws ] );
            }
          }
          else
          {
            for( std::vector< uint32_t >::const_iterator it_vws = vws_per_point[ candidate_point].begin(); it_vws != vws_per_point[ candidate_point].end(); ++it_vws, ++counter )
            {
              low_dim_vw_ids[counter] = parents_at_level_3[ *it_vws ];
              low_dim_vw.insert( parents_at_level_3[ *it_vws ] );
            }
          }
            
          // try to find a new correspondence
          nearest_neighbors_multiple nn_exp;
          for( std::set< uint32_t >::const_iterator activated_vw = low_dim_vw.begin(); activated_vw != low_dim_vw.end(); ++activated_vw )
          {
            counter = 0;
            for( std::vector< uint32_t >::const_iterator it_desc = desc_per_point[ candidate_point].begin(); it_desc != desc_per_point[ candidate_point].end(); ++it_desc, ++counter )
            {
              if( low_dim_vw_ids[counter] != *activated_vw )
                continue;
              
              for( std::vector< uint32_t >::const_iterator feature_it = features_per_vw[ *activated_vw ].begin(); feature_it != features_per_vw[ * activated_vw ].end(); ++feature_it )
              {
                int dist = 0;
                int x = 0;
                for( uint32_t jj=0; jj<128; ++jj )
                {
                  x = ( (int) descriptors[*feature_it][jj] ) - ( (int) all_descriptors[ *it_desc + jj ] );
                  dist += x*x;
                }
                
                nn_exp.update( *feature_it, dist );
              }
            }
          }
            
          if( nn_exp.dist1 >= 0 )
          {
            if( nn_exp.dist2 >= 0 )
            {
              if( nn_exp.get_ratio() < nn_ratio_3D_to_2D )
              {
                // we don't want to overwrite 2D-to-3D correspondence, only 3D-to-2D correspondence 
                if( !feature_in_correspondence[ nn_exp.nn_idx1 ] )
                {
                  // no existing correspondence
                  corr_3D_to_2D.insert( std::make_pair( candidate_point, std::make_pair( nn_exp.nn_idx1, nn_exp.dist1 ) ) );
                  feature_in_correspondence[ nn_exp.nn_idx1 ] = true;
                  point_per_feature[ nn_exp.nn_idx1 ] = candidate_point;
                }
                else if( point_per_feature[ nn_exp.nn_idx1 ] != UINT32_MAX )
                {
                  // only 3D-to-2D correspondence
                  // overwrite if the absolute SIFT distance is smaller 
                  
                  // first, we to get an iterator to the corresponding 3D point
                  map_it_3D = corr_3D_to_2D.find( point_per_feature[ nn_exp.nn_idx1 ] );
                  
                  // this has to exist, otherwise we would not have labeled it as having a correspondence
                  if( map_it_3D->second.second > nn_exp.dist1 )
                  {
                    // update the correspondence
                    
                    // which means we first have to remove the old one
                    corr_3D_to_2D.erase( map_it_3D );
                    
                    point_per_feature[ nn_exp.nn_idx1 ] = candidate_point;
                    
                    corr_3D_to_2D.insert( std::make_pair( candidate_point, std::make_pair( nn_exp.nn_idx1, nn_exp.dist1 ) ) );
                  }
                  
                }
              }
            }
          }
            
          if( corr_3D_to_2D.size() >= max_cor_early_term )
          {
            nb_considered_points = nb_considered_points_counter;
            break;
          }
        }
      }
    }

    
    if( nb_considered_points == 0 )
      nb_considered_points = nb_loaded_keypoints;
    
    
    ////
    // Establish the correspondences needed for RANSAC-based pose estimation
    
    if( ransac_filter == 1 )  
    {
      ////
      // Apply the RANSAC Pre-filter
      uint32_t max_set_size = 0;
     
      uint32_t nb_found_corr = (uint32_t) corr_3D_to_2D.size();
      
      // map found points to index range 0...N (actually the other direction)
      std::vector< uint32_t > index_to_point( nb_found_corr, 0 );
      
      // go through all points found as correspondences, establish edges in the form of images
      std::map< uint32_t, std::list< uint32_t > > image_edges;
      image_edges.clear();
      
      uint32_t point_counter = 0;
      
      std::map< uint32_t, std::list< uint32_t > >::iterator it;
      
      for( map_it_3D = corr_3D_to_2D.begin(); map_it_3D != corr_3D_to_2D.end(); ++map_it_3D, ++point_counter )
      {
        
        index_to_point[ point_counter ] = map_it_3D->first;
        
        for( std::set< uint32_t >::const_iterator it_images_point = images_per_point[ map_it_3D->first ].begin(); it_images_point != images_per_point[ map_it_3D->first ].end(); ++it_images_point )
        {
          it = image_edges.find( *it_images_point );
          
          if( it != image_edges.end() )
            it->second.push_back( point_counter );
          else
          {
            std::list< uint32_t > new_edge( 1, point_counter );
            image_edges.insert( std::make_pair< uint32_t, std::list< uint32_t > >( *it_images_point, new_edge ) );
          }
        }
      }
      
      ////
      // now find connected components, for every point, keep track which points belong to which connected component
      // also keep track of which images have been used
      
      std::vector< int > cc_per_corr( nb_found_corr, -1 );
      int current_cc = -1;
      int max_cc = -1;
      uint32_t size_current_cc = 0;
      
      std::set< uint32_t >::const_iterator img_it;
      
      for( point_counter = 0; point_counter < nb_found_corr; ++point_counter )
      {
        if( cc_per_corr[ point_counter ] < 0 )
        {
          // start new cc
          ++current_cc;
          size_current_cc = 0;
        
          std::queue< uint32_t > point_queue;
          point_queue.push( point_counter );
          
          // breadth first search for remaining points in connected component
          while( !point_queue.empty() )
          {
            uint32_t curr_point_id = point_queue.front();
            point_queue.pop();
            
            if( cc_per_corr[ curr_point_id ] < 0 )
            {
              cc_per_corr[ curr_point_id ] = current_cc;
              ++size_current_cc;
              
              // add all points in images visible by this point
              for( std::set< uint32_t >::const_iterator it_images_point = images_per_point[ index_to_point[ curr_point_id ] ].begin(); it_images_point != images_per_point[ index_to_point[ curr_point_id ] ].end(); ++it_images_point )
              {
                it = image_edges.find( *it_images_point );
                
                for( std::list< uint32_t >::const_iterator p_it = it->second.begin(); p_it != it->second.end(); ++p_it )
                {
                  if( cc_per_corr[ *p_it ] < 0 )
                    point_queue.push( *p_it );
                }
                
                // clear the image, we do not the this multi-edge anymore
                it->second.clear();
              }
            }
          }
          
          if( size_current_cc > max_set_size )
          {
            max_set_size = size_current_cc;
            max_cc = current_cc;
          }
         
        }
      }
      
      
      ////
      // now generate the correspondences
      
      c2D.reserve( 2*max_set_size );
      c3D.reserve( 3*max_set_size );
      
      point_counter = 0;
      for( map_it_3D = corr_3D_to_2D.begin(); map_it_3D != corr_3D_to_2D.end(); ++map_it_3D, ++point_counter )
      {
        if( cc_per_corr[ point_counter ] == max_cc )
        {
          c2D.push_back(keypoints[map_it_3D->second.first].x);
          c2D.push_back(keypoints[map_it_3D->second.first].y);
          
          c3D.push_back( points3D[map_it_3D->first][0] );
          c3D.push_back( points3D[map_it_3D->first][1] );
          c3D.push_back( points3D[map_it_3D->first][2] );
          
          final_correspondences.push_back( std::make_pair( map_it_3D->second.first, map_it_3D->first ) );
        }
      }
 
    }
    else
    {
      // normal correspondence computation
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
    }
      
    timer.Stop();
    std::cout << " computed correspondences in " << timer.GetElapsedTimeAsString() << ", considering " << nb_considered_points << " features " << " ( " << double(nb_considered_points) / double(nb_loaded_keypoints) * 100.0 << " % ) " << std::endl;
    corr_time = timer.GetElapsedTime();
    
   
    ////
    // do the pose verification using RANSAC
      
    RANSAC::computation_type = P6pt;
    RANSAC::stop_after_n_secs = true;
    RANSAC::max_time = ransac_max_time;
    RANSAC::error = 10.0f; // for P6pt this is the SQUARED reprojection error in pixels
    RANSAC ransac_solver;
    
    uint32_t nb_corr = c2D.size() / 2;
    
  
    std::cout << " applying RANSAC on " << nb_corr << " correspondences out of " << corr_3D_to_2D.size() << std::endl;
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
    
    Util::Math::ProjMatrix proj_matrix = ransac_solver.get_projection_matrix();
    Util::Math::ProjMatrix proj_matrix2 = ransac_solver.get_projection_matrix();
    
    // decompose the projection matrix
    Util::Math::Matrix3x3 Rot, K;
    proj_matrix.decompose( K, Rot );
    proj_matrix.computeInverse();
    proj_matrix.computeCenter();
    std::cout << " camera calibration: " << K << std::endl;
    std::cout << " camera rotation: " << Rot << std::endl;
    std::cout << " camera position: " << proj_matrix.m_center << std::endl;
    
    // write details to the details file. per line:
    // #inlier #correspondences time for vw time for corr time for ransac
    ofs << inlier.size() << " " << nb_corr << " " << vw_time << " " << corr_time << " " << RANSAC_time << " " << all_timer.GetElapsedTime() << std::endl;
    
    std::cout << "#########################" << std::endl;
    
    ////
    // Update statistics
    
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
    
    std::cout << std::endl << std::endl << " registered so far: " << registered << " / " << i+1 << std::endl;
    std::cout << " average time needed to compute the correspondences: registered: " << avrg_cor_computation_time_accepted << " rejected: " << avrg_cor_computation_time_rejected << std::endl;
    std::cout << "avrg. registration time: " << avrg_reg_time << " ( " << registered << " , avrg. inlier-ratio: " << mean_inlier_ratio_accepted << ", avrg. nb correspondences : " << mean_nb_correspondences_accepted << " ) avrg. rejection time: " << avrg_reject_time << " ( " << N_reject << ", avrg. inlier-ratio : " << mean_inlier_ratio_rejected << " avrg. nb correspondences : " << mean_nb_correspondences_rejected << " ) " << std::endl << std::endl;
  }
  
  ////
  // Overall statistics
  
  std::cout << std::endl << "#############################" << std::endl;
  std::cout << " total number registered: " << registered << " / " << nb_keyfiles << std::endl;
  std::cout << " average time for computing the vw assignments          : " << avrg_vw_time << " s for " << avrg_nb_features << " features (on average)" << std::endl;
  std::cout << " average time for succesfully registering image         : " << avrg_reg_time << " s " << std::endl;
  std::cout << " average time for rejecting an query image              : " << avrg_reject_time << " s " << std::endl;
  std::cout << " average inlier-ratio (registered)                      : " << mean_inlier_ratio_accepted << std::endl;
  std::cout << " average inlier-ratio (rejected)                        : " << mean_inlier_ratio_rejected << std::endl;
  std::cout << " average nb correspondences (registered)                : " << mean_nb_correspondences_accepted << std::endl;
  std::cout << " average nb correspondences (rejected)                  : " << mean_nb_correspondences_rejected << std::endl;
  std::cout << " average nb features        (registered)                : " << mean_nb_features_accepted << std::endl;
  std::cout << " average nb features        (rejected)                  : " << mean_nb_features_rejected << std::endl;
  std::cout << " avrg. time to compute the correspondences (registered) : " << avrg_cor_computation_time_accepted << std::endl;
  std::cout << " avrg. time to compute the correspondences (rejected)   : " << avrg_cor_computation_time_rejected << std::endl;
  std::cout << " avrg. time for RANSAC (registered)                     : " << avrg_RANSAC_time_registered << std::endl;
  std::cout << " avrg. time for RANSAC (rejected)                       : " << avrg_RANSAC_time_rejected << std::endl;
  std::cout << " minimum inlier-ratio for RANSAC                        : " << min_inlier << std::endl;
  std::cout << " stop after n correspondences                           : " << max_cor_early_term << std::endl;
  std::cout << " model consists of                                      : " << nb_descriptors << " ";
  std::cout << "unsigned char descriptors " << std::endl;
  if( ransac_filter > 0 )
    std::cout << " using the RANSAC prefiler " << ransac_filter << std::endl;
  if( use_image_set_cover )
    std::cout << " using the set cover with the " << consider_K_nearest_cams << " nearest cams " << std::endl;
  else
    std::cout << " Using the original images " << std::endl;
  if( prioritization_strategy == 0 )
    std::cout << " Combined Prioritization " << std::endl;
  else if( prioritization_strategy == 1 )
    std::cout << " Direct Prioritization " << std::endl;
  else if( prioritization_strategy == 2 )
    std::cout << " Afterwards Prioritization " << std::endl;
  if( filter_points )
    std::cout << " Filtering 3D points before 3D-to-2D matching " << std::endl;
  else
    std::cout << " NOT Filtering 3D points before 3D-to-2D matching " << std::endl;
  std::cout << "#############################" << std::endl;
  
  ofs.close();

  // delete kd-trees
  for( uint32_t i=0; i<nb_connected_components; ++i )
  {
    for( uint32_t j=0; j<nb_points_per_component[i]; ++j )
      points_per_component[i][j] = 0;
    
    delete [] indices_per_component[i];
    indices_per_component[i] = 0;
    
    delete [] points_per_component[i];
    points_per_component[i] = 0;
  }
  
  for( uint32_t i=0; i<nb_3D_points; ++i )
  {
    if( points3D[i] != 0 )
      delete [] points3D[i];
    points3D[i] = 0;
  }
  delete [] points3D;
  points3D = 0;
  
  for( uint32_t i=0; i<nb_connected_components; ++i )
  {
    delete kd_trees[i];
    kd_trees[i] = 0;
  }
  
  delete kd_tree;
  kd_tree = 0;
  
  delete [] indices;
  delete [] distances;
  
  delete [] parents_at_level_2;
  parents_at_level_2 = 0;
  
  delete [] parents_at_level_3;
  parents_at_level_3 = 0;
  
  annClose();
  

  for( uint32_t i=0; i<nb_3D_points; ++i )
    images_per_point[i].clear();
  images_per_point.clear();
  
  return 0;
}

