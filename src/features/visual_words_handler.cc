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

#include "visual_words_handler.hh"



visual_words_handler::visual_words_handler( )
{
  mMethod = 2;
  mNbTrees = 8;
  mNbPath = 1;
  mFlannIndexType = 2;
  mBranching = 10;
  mTargetPrecision = 0.9f;
  mBuildWeight = 0.1f;
  mMemoryWeight = 0.0f;
  mSampleFraction = 0.1f;
  mNbNearestNeighbors = 2;
  mNbVisualWords = 100000;
  
  mFlannIndex = 0;
  
  mMaxDescriptors = 30000;
  
  mClusterCentersFlann.data = new float[mNbVisualWords * 128];

  mClusterCentersFlann.rows = mNbVisualWords;
  mClusterCentersFlann.cols = 128;
  
}

//---------------------------------------------------

visual_words_handler::~visual_words_handler( )
{
  if( mClusterCentersFlann.data != 0 )
    mClusterCentersFlann.free();
  if( mFlannAssignments.data != 0 )
    mFlannAssignments.free();
  mFlannAssignmentsKNN.free();
  if( mFlannDistances.data != 0 )
    mFlannDistances.free();
  if( mFlannDistancesKNN.data != 0 )
    mFlannDistancesKNN.free();
  if( mFlannFeatures.data != 0 )
    mFlannFeatures.free();
  
  if( mFlannIndex != 0 )
    delete mFlannIndex;
  mFlannIndex = 0;
  
}

//---------------------------------------------------

void visual_words_handler::set_method( const std::string &method )
{
  if ( method == "flann" )
    mMethod = 2;
  else if ( method == "linear" )
    mMethod = 0;
}

//---------------------------------------------------

void visual_words_handler::set_nb_trees( int nb_trees )
{
  mNbTrees = nb_trees;
}


//---------------------------------------------------

void visual_words_handler::set_nb_paths( int nb_path )
{
  mNbPath = nb_path;
}

//---------------------------------------------------

void visual_words_handler::set_flann_type( const std::string &type )
{
  if( type == "auto" )
    mFlannIndexType = 0;
  else if( type == "hkmeans" )
    mFlannIndexType = 1;
  else if( type == "randomkd" )
    mFlannIndexType = 2;
}

//---------------------------------------------------

void visual_words_handler::set_nb_visual_words( uint32_t nb_vw )
{
  if( nb_vw != mNbVisualWords )
  {
    mClusterCentersFlann.free();
    mClusterCentersFlann.data = 0; 
    mClusterCentersFlann.data = new float[nb_vw *128 ];
    mClusterCentersFlann.rows = (size_t) nb_vw;
  }
  mNbVisualWords = nb_vw;
}

//---------------------------------------------------

void visual_words_handler::set_branching(int nb_branches)
{
  mBranching = nb_branches;
}


//---------------------------------------------------

void visual_words_handler::set_target_precision( float precision )
{
  mTargetPrecision = precision;
}
    
//---------------------------------------------------

void visual_words_handler::set_build_weight( float weight )
{
  mBuildWeight = weight;
}
    
//---------------------------------------------------

void visual_words_handler::set_memory_weight( float weight )
{
  mMemoryWeight = weight;
}
    
//---------------------------------------------------

void visual_words_handler::set_sample_size( float rel_size )
{
  mSampleFraction = rel_size;
}

//---------------------------------------------------

void visual_words_handler::set_nb_nearest_neighbors( int nb_neighbors )
{
//   // check if we have to resize the matrices in which flann stores the assignments 
//   // and the distances
//   if( nb_neighbors != mNbNearestNeighbors )
//   {
//     mNbNearestNeighbors = nb_neighbors;
//     mFlannAssignmentsKNN.free();
//     mFlannAssignmentsKNN.data = new int[ mNbNearestNeighbors * mMaxDescriptors ];
//     mFlannAssignmentsKNN.rows = mMaxDescriptors;
//     mFlannAssignmentsKNN.cols = mNbNearestNeighbors;
//     
//     mFlannDistancesKNN.free();
//     mFlannDistancesKNN.data = new float[ mNbNearestNeighbors * mMaxDescriptors ];
//     mFlannDistancesKNN.rows = mMaxDescriptors;
//     mFlannDistancesKNN.cols = mNbNearestNeighbors;
//   }
//   else
    mNbNearestNeighbors = nb_neighbors;
  
}

//---------------------------------------------------

void visual_words_handler::set_cluster_centers( std::vector< float > &cluster_centers )
{
  size_t nb_clusts = cluster_centers.size() / 128;
  if( nb_clusts != (size_t) mNbVisualWords )
  {
    // resize
    mClusterCentersFlann.free();
    mClusterCentersFlann.data = 0;
    mClusterCentersFlann.data = new float[nb_clusts *128 ];
    mClusterCentersFlann.rows = (size_t) nb_clusts;
    mNbVisualWords = (uint32_t) nb_clusts;
  }
  
  size_t index = 0;
  for( size_t i=0; i<nb_clusts; ++i )
  {
    for( size_t j=0; j<128; ++j )
    {
      mClusterCentersFlann[i][j] = cluster_centers[index+j];
    }
    index += 128;
  }
}

//---------------------------------------------------

void visual_words_handler::set_cluster_centers( std::string &cluster_file )
{
  std::ifstream ifs( cluster_file.c_str(), std::ios::in );
  if ( !ifs )
  {
      std::cerr << "[Visual_Words_Handler]: ERROR: Cannot read the clusters from " << cluster_file << std::endl;;
      return;
  }
  
  for( uint32_t i=0; i<mNbVisualWords; ++i )
  {
    float tmp_flt;
    for( uint32_t j=0; j<128; ++j )
    {
      ifs >> tmp_flt;
      (mClusterCentersFlann.data)[ 128*i + j ] = tmp_flt;
    }
  }
  
  ifs.close();
}

//---------------------------------------------------

std::vector<float> visual_words_handler::get_cluster_centers_flann()
{
  std::vector< float > cluster_centers_( mNbVisualWords * uint32_t(128), 0.0f );
  uint32_t index = 0;
  for( uint32_t i=0; i<mNbVisualWords; ++i )
  {
    for( uint32_t j=0; j<128; ++j )
      cluster_centers_[index+j] = (mClusterCentersFlann.data)[ index + j ];
    index += (uint32_t) 128;
  }

  return cluster_centers_;
}

//---------------------------------------------------

void visual_words_handler::get_cluster_centers_flann( std::vector< float > &centers_ )
{
  centers_.resize( mNbVisualWords * uint32_t(128), 0.0f );
  uint32_t index = 0;
  for( uint32_t i=0; i<mNbVisualWords; ++i )
  {
    for( uint32_t j=0; j<128; ++j )
      centers_[index+j] = (mClusterCentersFlann.data)[ index + j ];
    index += (uint32_t) 128;
  }
}

//---------------------------------------------------

bool visual_words_handler::load_trees_flann( std::string &cluster_file, std::string &tree_file )
{

  if( mFlannIndex != 0 )
    delete [] mFlannIndex;
  mFlannIndex = 0;
  
  // load the clusters from a text file
  {
    std::ifstream ifs( cluster_file.c_str(), std::ios::in );
    if ( !ifs )
    {
      std::cerr << "[Visual_Words_Handler]: ERROR: Cannot read the clusters from " << cluster_file << std::endl;
      return false;
    }
    
    for( uint32_t i=0; i<mNbVisualWords; ++i )
    {
      float tmp_flt;
      for( uint32_t j=0; j<128; ++j )
      {
		ifs >> tmp_flt;
		mClusterCentersFlann[i][j] = tmp_flt;
      }
    }
    
    ifs.close();
  }
  
  // check if the file for the tree exists
  bool exists = false;
  {
    std::ifstream in( tree_file.c_str() );
    exists = in.is_open();
    in.close();
  }
  
  if( exists )
  {
    // load the tree
    std::cout << "[Visual_Words_Handler]: Loading tree... " << std::endl;
    mFlannIndex = new flann::Index< flann::L2< float > >( mClusterCentersFlann, flann::SavedIndexParams( tree_file ) );    
    std::cout << "[Visual_Words_Handler]: Tree loaded " << std::endl;
  }
  else
  {
    // create the tree and save it
    std::cout << "[Visual_Words_Handler]: Building the tree... " << std::endl;
    
    
    // NOTE: Make sure that we allways generate the same tree for a set of features / cluster centers and are not 
    // afflicted by any calls to rand or srand
    srand(1);
    
    switch( mFlannIndexType )
    {
      case 0:
		mFlannIndex = new flann::Index< flann::L2< float > >( mClusterCentersFlann, flann::AutotunedIndexParams( mTargetPrecision, mBuildWeight, mMemoryWeight, mSampleFraction ) );
		break;
      case 1:
		mFlannIndex = new flann::Index< flann::L2< float > >( mClusterCentersFlann, flann::KMeansIndexParams( mBranching ) );
		break;
      case 2:
		mFlannIndex = new flann::Index< flann::L2< float > >( mClusterCentersFlann, flann::KDTreeIndexParams( mNbTrees ) );
		break;
	
    }
    mFlannIndex->buildIndex();
    mFlannIndex->save( tree_file );
    std::cout << "[Visual_Words_Handler]: Tree built. " << std::endl;
  }
  
  initialize();
  std::cout << "[Visual_Words_Handler]: Arrays initialized. " << std::endl;
  return true;
}

//---------------------------------------------------

void visual_words_handler::rebuild_flann_index()
{
  if( mFlannIndex != 0 )
    delete [] mFlannIndex;
  mFlannIndex = 0;
  
  std::cout << "[Visual_Words_Handler]: Building the tree... " << std::endl;

  
  // NOTE: Make sure that we allways generate the same tree for a set of features / cluster centers and are not 
  // afflicted by any calls to rand or srand
  srand(1);
  
  if( mMethod == 2 )
  {
    switch( mFlannIndexType )
    {
      case 0:
		mFlannIndex = new flann::Index< flann::L2< float > >( mClusterCentersFlann, flann::AutotunedIndexParams( mTargetPrecision, mBuildWeight, mMemoryWeight, mSampleFraction ) );
		break;
      case 1:
		mFlannIndex = new flann::Index< flann::L2< float > >( mClusterCentersFlann, flann::KMeansIndexParams( mBranching ) );
		break;
      case 2:
		mFlannIndex = new flann::Index< flann::L2< float > >( mClusterCentersFlann, flann::KDTreeIndexParams( mNbTrees ) );
		break;
	
    }
  }
  mFlannIndex->buildIndex();
  
  std::cout << "[Visual_Words_Handler]: Tree built. " << std::endl;
  initialize();
  std::cout << "[Visual_Words_Handler]: Arrays initialized. " << std::endl;
}

//---------------------------------------------------

bool visual_words_handler::create_flann_search_index( std::string &cluster_file )
{
  if( mFlannIndex != 0 )
    delete [] mFlannIndex;
  mFlannIndex = 0;
  
  // load the clusters from a text file
  {
    std::ifstream ifs( cluster_file.c_str(), std::ios::in );
    if ( !ifs )
    {
	std::cerr << "[Visual_Words_Handler]: ERROR: Cannot read the clusters from " << cluster_file << std::endl;;
	return false;
    }
    
    for( uint32_t i=0; i<mNbVisualWords; ++i )
    {
      float tmp_flt;
      for( uint32_t j=0; j<128; ++j )
      {
		ifs >> tmp_flt;
		mClusterCentersFlann[i][j] = tmp_flt;
      }
    }
    
    ifs.close();
  }

  // create the tree and save it
  std::cout << "[Visual_Words_Handler]: Building the tree for " <<  mNbVisualWords << " points ... " << std::endl;

  
  // NOTE: Make sure that we allways generate the same tree for a set of features / cluster centers and are not 
  // afflicted by any calls to rand or srand
  srand(1);
  
  if( mMethod == 2 )
  {
    switch( mFlannIndexType )
    {
      case 0:
		mFlannIndex = new flann::Index< flann::L2< float > >( mClusterCentersFlann, flann::AutotunedIndexParams( mTargetPrecision, mBuildWeight, mMemoryWeight, mSampleFraction ) );
		break;
      case 1:
		mFlannIndex = new flann::Index< flann::L2< float > >( mClusterCentersFlann, flann::KMeansIndexParams( mBranching ) );
		break;
      case 2:
		mFlannIndex = new flann::Index< flann::L2< float > >( mClusterCentersFlann, flann::KDTreeIndexParams( mNbTrees ) );
		break;
    }
  }
  
  mFlannIndex->buildIndex();
  std::cout << "[Visual_Words_Handler]: Tree built. " << std::endl;
  
  initialize();
  std::cout << "[Visual_Words_Handler]: Arrays initialized. " << std::endl;
  return true;
}

//---------------------------------------------------

bool visual_words_handler::create_flann_search_index( std::vector< float > &cluster_centers )
{
  if( mFlannIndex != 0 )
    delete [] mFlannIndex;
  mFlannIndex = 0;
  
  // copy the clusters
  size_t nb_clusts = cluster_centers.size() / 128;
  if( nb_clusts != (size_t) mNbVisualWords )
  {
    // resize
    mClusterCentersFlann.free();
    mClusterCentersFlann.data = new float[nb_clusts *128 ];
    mClusterCentersFlann.rows = (size_t) nb_clusts;
    mNbVisualWords = (uint32_t) nb_clusts;
  }
  
  size_t index = 0;
  for( size_t i=0; i<nb_clusts; ++i )
  {
//     size_t index = i *128;
    for( size_t j=0; j<128; ++j )
    {
	  mClusterCentersFlann[i][j] = cluster_centers[index+j];
    }
    index += 128;
  }
                                                 
  // create the tree and save it
  std::cout << "[Visual_Words_Handler]: Building the tree... " << std::endl;
  
  // NOTE: Make sure that we allways generate the same tree for a set of features / cluster centers and are not 
  // afflicted by any calls to rand or srand
  srand(1);
  
  if( mMethod == 2 )
  {
    switch( mFlannIndexType )
    {
      case 0:
		mFlannIndex = new flann::Index< flann::L2< float > >( mClusterCentersFlann, flann::AutotunedIndexParams( mTargetPrecision, mBuildWeight, mMemoryWeight, mSampleFraction ) );
		break;
      case 1:
		mFlannIndex = new flann::Index< flann::L2< float > >( mClusterCentersFlann, flann::KMeansIndexParams( mBranching ) );
		break;
      case 2:
		mFlannIndex = new flann::Index< flann::L2< float > >( mClusterCentersFlann, flann::KDTreeIndexParams( mNbTrees ) );
		break;
    }
  }

  mFlannIndex->buildIndex();
  std::cout << "[Visual_Words_Handler]: Tree built. " << std::endl;
  
  initialize();
  std::cout << "[Visual_Words_Handler]: Arrays initialized. " << std::endl;
  return true;
}

//---------------------------------------------------

void visual_words_handler::get_parents_at_level_L( int L, int* &ids )
{
  if( mMethod == 2 && mFlannIndexType == 1 && mFlannIndex != 0 )
  {
    if( dynamic_cast< flann::KMeansIndex< flann::L2< float > >* >( mFlannIndex->getIndex() ) != 0 )
    {
      flann::KMeansIndex< flann::L2< float > >* tree_ = dynamic_cast< flann::KMeansIndex< flann::L2< float > >* >( mFlannIndex->getIndex() );
      tree_->getClusterCentersOnLevelL( L, ids );
      tree_ = 0;
    }
    else
      std::cout << "[Visual_Words_Handler]: ERROR: Could not get the parents at the higher levels. " << std::endl;
  }
}

//---------------------------------------------------
    
bool visual_words_handler::assign_visual_words_uchar( std::vector< unsigned char > &descriptors, uint32_t nb_descriptors, std::vector< uint32_t > &assignments )
{
  if( assignments.size() < (size_t) nb_descriptors )
    assignments.resize(nb_descriptors);
    
  if( mMethod == 0 )
  {
    // do a linear search
    float vec[128];
    uint32_t index = 0;
    for( uint32_t i=0; i<nb_descriptors; ++i )
    {
      // copy the descriptor
//       uint32_t index = i*128;
      for( uint32_t j=0; j<128; ++j )
	vec[j] = (float) descriptors[index+j];
      
      // find the nearest neighbor
      float max_dist = 1e20f;
      float dist = 0.0f;
      float x;
      uint32_t index2 = 0;
      for( uint32_t j=0; j<mNbVisualWords; ++j )
      {
	dist = 0.0f;
// 	index2 = j*128;
	for( uint32_t k=0; k<128; ++k )
	{
	  x = vec[k] - mClusterCentersFlann.data[index2+k];
	  dist += x*x;
	}
	
	if( dist < max_dist )
	{
	  max_dist = dist;
	  assignments[i] = j;
	}
	index2 += 128;
      }
      index += 128;
    }
    
    return true;
  }
  else if( mMethod == 2 )
  {
    // FLANN
//     std::cout << " FLANN : " << nb_descriptors << " paths : " << mNbPath << " index type : " << mFlannIndexType << std::endl;
    if( mFlannIndex == 0 )
      return false;
    // copy the descriptors
    resize( nb_descriptors );
    uint32_t index = 0;
    for( uint32_t i=0; i<nb_descriptors; ++i )
    {
      for( uint32_t j=0; j<128; ++j )
      {
		mFlannFeatures[i][j] = (float) descriptors[index+j];
      }
      index += 128;
    }
    
    // compute the assignments
    if( mFlannIndexType == 0 )
    {
      flann::SearchParams params( FLANN_CHECKS_AUTOTUNED );
      mFlannIndex->knnSearch( mFlannFeatures, mFlannAssignments, mFlannDistances, 1, params );
    }
    else
    {
      flann::SearchParams params( mNbPath );
      mFlannIndex->knnSearch( mFlannFeatures, mFlannAssignments, mFlannDistances, 1, params );
    }
    
    // copy the assignments
    for( uint32_t i=0; i<nb_descriptors; ++i )
    {
      assignments[i] = (uint32_t) mFlannAssignments.data[i];
    }
    
    return true;
  }
  
  return false;
}



//---------------------------------------------------
    
bool visual_words_handler::assign_visual_words_ucharv( std::vector< unsigned char* > &descriptors, uint32_t nb_descriptors, std::vector< uint32_t > &assignments )
{
  if( assignments.size() < (size_t) nb_descriptors )
    assignments.resize(nb_descriptors);
  
  if( mMethod == 0 )
  {
    // do a linear search
    float vec[128];
    for( uint32_t i=0; i<nb_descriptors; ++i )
    {
      // copy the descriptor
      for( uint32_t j=0; j<128; ++j )
	vec[j] = (float) descriptors[i][j];
      
      // find the nearest neighbor
      float max_dist = 1e20f;
      float dist = 0.0f;
      float x;
      uint32_t index2 = 0;
      for( uint32_t j=0; j<mNbVisualWords; ++j )
      {
	dist = 0.0f;
// 	index2 = j*128;
	for( uint32_t k=0; k<128; ++k )
	{
	  x = vec[k] - mClusterCentersFlann.data[index2+k];
	  dist += x*x;
	}
	
	if( dist < max_dist )
	{
	  max_dist = dist;
	  assignments[i] = j;
	}
	index2 += 128;
      }
//       std::cout << " Linear: " << i << " -> " << assignments[i] << " dist :  " << max_dist << std::endl;
    }
    
    return true;
  }
  else if( mMethod == 2 )
  {
    // FLANN
    if( mFlannIndex == 0 )
      return false;
    
    // copy the descriptors
    resize( nb_descriptors );
    for( uint32_t i=0; i<nb_descriptors; ++i )
    {
      for( uint32_t j=0; j<128; ++j )
      {
		mFlannFeatures[i][j] = (float) descriptors[i][j];
      }
    }
    
    // compute the assignments
    if( mFlannIndexType == 0 )
    {
      flann::SearchParams params( FLANN_CHECKS_AUTOTUNED );
      mFlannIndex->knnSearch( mFlannFeatures, mFlannAssignments, mFlannDistances, 1, params );
    }
    else
    {
      flann::SearchParams params( mNbPath );
      mFlannIndex->knnSearch( mFlannFeatures, mFlannAssignments, mFlannDistances, 1, params );
    }
    
    // copy the assignments
    for( uint32_t i=0; i<nb_descriptors; ++i )
    {
      assignments[i] = (uint32_t) mFlannAssignments.data[i];
    }
    
    return true;
  }
  
  return false;
}

//---------------------------------------------------

bool visual_words_handler::assign_visual_words_float( std::vector< float > &descriptors, uint32_t nb_descriptors, std::vector< uint32_t > &assignments )
{
  if( assignments.size() < (size_t) nb_descriptors )
    assignments.resize(nb_descriptors);
  
  if( mMethod == 0 )
  {
    // do a linear search
    float vec[128];
    uint32_t index = 0;
    for( uint32_t i=0; i<nb_descriptors; ++i )
    {
      // copy the descriptor
//       uint32_t index = i*128;
      for( uint32_t j=0; j<128; ++j )
	vec[j] = descriptors[index+j];
      
      // find the nearest neighbor
      float max_dist = 1e20f;
      float dist = 0.0f;
      float x;
      uint32_t index2 = 0;
      for( uint32_t j=0; j<mNbVisualWords; ++j )
      {
	dist = 0.0f;
// 	index2 = j*128;
	for( uint32_t k=0; k<128; ++k )
	{
	  x = vec[k] - mClusterCentersFlann.data[index2+k];
	  dist += x*x;
	}
	
	if( dist < max_dist )
	{
	  max_dist = dist;
	  assignments[i] = j;
	}
	index2 += 128;
      }
      
      index += 128;
    }
    
    return true;
  }
  else if( mMethod == 2 )
  {
    // FLANN
    
    if( mFlannIndex == 0 )
      return false;
    
    // copy the descriptors
    resize( nb_descriptors );
    uint32_t index = 0;
    for( uint32_t i=0; i<nb_descriptors; ++i )
    {
      for( uint32_t j=0; j<128; ++j )
      {
		mFlannFeatures[i][j] = descriptors[index+j];
      }
      index += 128;
    }
    
    // compute the assignments
    if( mFlannIndexType == 0 )
    {
      flann::SearchParams params( FLANN_CHECKS_AUTOTUNED );
      mFlannIndex->knnSearch( mFlannFeatures, mFlannAssignments, mFlannDistances, 1, params );
    }
    else
    {
      flann::SearchParams params( mNbPath );
      mFlannIndex->knnSearch( mFlannFeatures, mFlannAssignments, mFlannDistances, 1, params );
    }
    
    // copy the assignments
    for( uint32_t i=0; i<nb_descriptors; ++i )
      assignments[i] = (uint32_t) mFlannAssignments.data[i];
    
    return true;
  }
  
  
  return false;
}

//---------------------------------------------------
    
bool visual_words_handler::assign_visual_words_uchar( std::vector< unsigned char > &descriptors, uint32_t nb_descriptors, std::vector< uint32_t > &assignments, std::vector< float > &distances )
{
  if( assignments.size() < (size_t) nb_descriptors )
    assignments.resize(nb_descriptors);
    
  if( mMethod == 0 )
  {
    // do a linear search
    float vec[128];
    uint32_t index = 0;
    for( uint32_t i=0; i<nb_descriptors; ++i )
    {
      // copy the descriptor
//       uint32_t index = i*128;
      for( uint32_t j=0; j<128; ++j )
	vec[j] = (float) descriptors[index+j];
      
      // find the nearest neighbor
      float max_dist = 1e20f;
      float dist = 0.0f;
      float x;
      uint32_t index2 = 0;
      for( uint32_t j=0; j<mNbVisualWords; ++j )
      {
	dist = 0.0f;
// 	index2 = j*128;
	for( uint32_t k=0; k<128; ++k )
	{
	  x = vec[k] - mClusterCentersFlann.data[index2+k];
	  dist += x*x;
	}
	
	if( dist < max_dist )
	{
	  max_dist = dist;
	  assignments[i] = j;
	}
	index2 += 128;
      }
      distances[i] = max_dist;
      index += 128;
    }
    
    return true;
  }
  else if( mMethod == 2 )
  {
    // FLANN
//     std::cout << " FLANN : " << nb_descriptors << " paths : " << mNbPath << " index type : " << mFlannIndexType << std::endl;
    if( mFlannIndex == 0 )
      return false;
    // copy the descriptors
    resize( nb_descriptors );
    uint32_t index = 0;
    for( uint32_t i=0; i<nb_descriptors; ++i )
    {
      for( uint32_t j=0; j<128; ++j )
      {
		mFlannFeatures[i][j] = (float) descriptors[index+j];
      }
      index += 128;
//       std::cout << std::endl;
    }
    
    // compute the assignments
    if( mFlannIndexType == 0 )
    {
      flann::SearchParams params( FLANN_CHECKS_AUTOTUNED );
      mFlannIndex->knnSearch( mFlannFeatures, mFlannAssignments, mFlannDistances, 1, params );
    }
    else
    {
      flann::SearchParams params( mNbPath );
      mFlannIndex->knnSearch( mFlannFeatures, mFlannAssignments, mFlannDistances, 1, params );
    }
    
    // copy the assignments
    for( uint32_t i=0; i<nb_descriptors; ++i )
    {
      assignments[i] = (uint32_t) mFlannAssignments.data[i];
      distances[i] = mFlannDistances.data[i];
    }
    
    return true;
  }

  
  return false;
}



//---------------------------------------------------
    
bool visual_words_handler::k_nn_search_flann_uchar( std::vector< unsigned char > &descriptors, uint32_t nb_descriptors, std::vector< uint32_t > &assignments, std::vector< float > &distances )
{
  if( mFlannIndex == 0 )
      return false;
  
  // check if we need to resize
  if( assignments.size() < (size_t) nb_descriptors )
    assignments.resize(nb_descriptors);
  

  // FLANN
  // copy the descriptors
  resize( nb_descriptors );
  uint32_t index = 0;
  for( uint32_t i=0; i<nb_descriptors; ++i )
  {
    for( uint32_t j=0; j<128; ++j )
    {
      mFlannFeatures[i][j] = (float) descriptors[index+j];
    }
    index += 128;
  }
  
  // compute the assignments
  if( mFlannIndexType == 0 )
  {
    flann::SearchParams params( FLANN_CHECKS_AUTOTUNED );
    mFlannIndex->knnSearch( mFlannFeatures, mFlannAssignmentsKNN, mFlannDistancesKNN, mNbNearestNeighbors, params );
  }
  else
  {
    flann::SearchParams params( mNbPath );
    mFlannIndex->knnSearch( mFlannFeatures, mFlannAssignmentsKNN, mFlannDistancesKNN, mNbNearestNeighbors, params );
  }
  
  // copy the assignments
  uint32_t k = (uint32_t) mNbNearestNeighbors;
  for( uint32_t i=0; i<nb_descriptors; ++i )
  {
    for( uint32_t j=0; j<k; ++j )
    {
      assignments[k*i+j] = (uint32_t) mFlannAssignmentsKNN.data[k*i+j];
      distances[k*i+j] = (uint32_t) mFlannDistancesKNN.data[k*i+j];
    }
  }
  
  return true;
}



//---------------------------------------------------
    
bool visual_words_handler::k_nn_search_flann_ucharv( std::vector< unsigned char* > &descriptors, uint32_t nb_descriptors, std::vector< uint32_t > &assignments, std::vector< float > &distances )
{
  if( mFlannIndex == 0 )
      return false;
  
  if( assignments.size() < (size_t) nb_descriptors )
    assignments.resize(nb_descriptors);
    
  // FLANN
  // copy the descriptors
  resize( nb_descriptors );
  for( uint32_t i=0; i<nb_descriptors; ++i )
  {
    for( uint32_t j=0; j<128; ++j )
    {
      mFlannFeatures[i][j] = (float) descriptors[i][j];
    }
  }
  
  // compute the assignments
  if( mFlannIndexType == 0 )
  {
    flann::SearchParams params( FLANN_CHECKS_AUTOTUNED );
    mFlannIndex->knnSearch( mFlannFeatures, mFlannAssignmentsKNN, mFlannDistancesKNN, mNbNearestNeighbors, params );
  }
  else
  {
    flann::SearchParams params( mNbPath );
    mFlannIndex->knnSearch( mFlannFeatures, mFlannAssignmentsKNN, mFlannDistancesKNN, mNbNearestNeighbors, params );
  }
  
  // copy the assignments
  uint32_t k = (uint32_t) mNbNearestNeighbors;
  for( uint32_t i=0; i<nb_descriptors; ++i )
  {
    for( uint32_t j=0; j<k; ++j )
    {
      assignments[k*i+j] = (uint32_t) mFlannAssignmentsKNN.data[k*i+j];
      distances[k*i+j] = (uint32_t) mFlannDistancesKNN.data[k*i+j];
    }
  }
  
  return true;
}

//---------------------------------------------------

bool visual_words_handler::k_nn_search_flann_float( std::vector< float > &descriptors, uint32_t nb_descriptors, std::vector< uint32_t > &assignments, std::vector< float > &distances )
{
  if( mFlannIndex == 0 )
    return false;
  
  if( assignments.size() < (size_t) nb_descriptors )
    assignments.resize(nb_descriptors);
  
  
  // FLANN
  // copy the descriptors
  resize( nb_descriptors );
  uint32_t index = 0;
  for( uint32_t i=0; i<nb_descriptors; ++i )
  {
    for( uint32_t j=0; j<128; ++j )
    {
      mFlannFeatures[i][j] = descriptors[index+j];
    }
    index += 128;
  }
  
  // compute the assignments
  if( mFlannIndexType == 0 )
  {
    flann::SearchParams params( FLANN_CHECKS_AUTOTUNED );
    mFlannIndex->knnSearch( mFlannFeatures, mFlannAssignmentsKNN, mFlannDistancesKNN, mNbNearestNeighbors, params );
  }
  else
  {
    flann::SearchParams params( mNbPath );
    mFlannIndex->knnSearch( mFlannFeatures, mFlannAssignmentsKNN, mFlannDistancesKNN, mNbNearestNeighbors, params );
  }
  
  // copy the assignments
  uint32_t k = (uint32_t) mNbNearestNeighbors;
  for( uint32_t i=0; i<nb_descriptors; ++i )
  {
    for( uint32_t j=0; j<k; ++j )
    {
      assignments[k*i+j] = (uint32_t) mFlannAssignmentsKNN.data[k*i+j];
      distances[k*i+j] = (uint32_t) mFlannDistancesKNN.data[k*i+j];
    }
  }
  
  return true;
}

//---------------------------------------------------

void visual_words_handler::load_from_file( std::string &filename, uint32_t nb_assignments, std::vector< uint32_t > & assignments )
{
  if( assignments.size() < (size_t) nb_assignments )
    assignments.resize(nb_assignments,0);
  
  std::cout << "[Visual_Words_Handler]: Loading the visual words from " << filename << std::endl;
  std::ifstream ifs( filename.c_str(), std::ios::in );
  
  if ( !ifs )
  {
    std::cerr << "[Visual_Words_Handler]: ERROR: Cannot read the assignments from " << filename << std::endl;;
    return;
  }
      
  for( uint32_t j=0; j<nb_assignments; ++j )
  {
    ifs >> assignments[j];
  }

  ifs.close();
  
  std::cout << "[Visual_Words_Handler]: Done loading the visual words. " << std::endl;
}

//---------------------------------------------------

/** loads assignments from a file. The difference to the function above is that the filename parameter is a prefix of the real filename. 
    The rest of the filename is computed based on the parameters specified.
**/
void visual_words_handler::load_from_file_prefix( std::string &filename, uint32_t nb_assignments, std::vector< uint32_t > & assignments )
{
  if( assignments.size() < (size_t) nb_assignments )
    assignments.resize(nb_assignments,0);
  
  std::string vw_file( filename );
  if( mNbVisualWords == 1000000 )
    vw_file.append("1M.");
  else if( mNbVisualWords == 100000 )
    vw_file.append("100k.");
  if( mNbTrees == 256 )
    vw_file.append("256.");
  else if( mNbTrees == 128 )
    vw_file.append("128.");
  else if( mNbTrees == 64 )
    vw_file.append("64.");
  else if( mNbTrees == 16 )
    vw_file.append("16.");
  else if( mNbTrees == 8 )
    vw_file.append("8.");
  else if( mNbTrees == 4 )
    vw_file.append("4.");
  else if( mNbTrees == 2 )
    vw_file.append("2.");
  vw_file.append("vw");
  
  std::cout << "[Visual_Words_Handler]: loading the visual words from " << vw_file << std::endl;
  
  std::ifstream ifs( vw_file.c_str(), std::ios::in );
  
  if ( !ifs )
  {
    std::cerr << "[Visual_Words_Handler]: ERROR: Cannot read the assignments from " << vw_file << std::endl;;
    return;
  }
      
  for( uint32_t j=0; j<nb_assignments; ++j )
  {
    ifs >> assignments[j];
  }

  ifs.close();
  
  std::cout << "[Visual_Words_Handler]: Done loading the visual words. " << std::endl;
}

//---------------------------------------------------

void visual_words_handler::assign_and_save_uchar( std::string &filename, std::vector< unsigned char > &descriptors, uint32_t nb_descriptors )
{
  std::vector< uint32_t > assignments( nb_descriptors, mNbVisualWords+1 );
  assign_visual_words_uchar( descriptors, nb_descriptors, assignments );
  
  std::ofstream ofs( filename.c_str(), std::ios::out );
  
  if ( !ofs )
  {
    std::cerr << "[Visual_Words_Handler]: ERROR: Cannot write the assignments to " << filename << std::endl;;
    return;
  }
  
  for( uint32_t i=0; i<nb_descriptors; ++i )
    ofs << assignments[i] << std::endl;
  
  ofs.close();
}

//---------------------------------------------------

void visual_words_handler::assign_and_save_float( std::string &filename, std::vector< float > &descriptors, uint32_t nb_descriptors )
{
  std::vector< uint32_t > assignments( nb_descriptors, mNbVisualWords+1 );
  assign_visual_words_float( descriptors, nb_descriptors, assignments );
  
  std::ofstream ofs( filename.c_str(), std::ios::out );
  
  if ( !ofs )
  {
    std::cerr << "[Visual_Words_Handler]: ERROR: Cannot write the assignments to " << filename << std::endl;;
    return;
  }
  
  for( uint32_t i=0; i<nb_descriptors; ++i )
    ofs << assignments[i] << std::endl;
  
  ofs.close();
}

//---------------------------------------------------

void visual_words_handler::assign_and_save_prefix_uchar( std::string &filename, std::vector< unsigned char > &descriptors, uint32_t nb_descriptors )
{
  std::vector< uint32_t > assignments( nb_descriptors, mNbVisualWords+1 );
  assign_visual_words_uchar( descriptors, nb_descriptors, assignments );
  
  std::string vw_file( filename );
  if( mNbVisualWords == 1000000 )
    vw_file.append("1M.");
  else if( mNbVisualWords == 100000 )
    vw_file.append("100k.");
  if( mNbTrees == 256 )
    vw_file.append("256.");
  else if( mNbTrees == 128 )
    vw_file.append("128.");
  else if( mNbTrees == 64 )
    vw_file.append("64.");
  else if( mNbTrees == 16 )
    vw_file.append("16.");
  else if( mNbTrees == 8 )
    vw_file.append("8.");
  else if( mNbTrees == 4 )
    vw_file.append("4.");
  else if( mNbTrees == 2 )
    vw_file.append("2.");
  vw_file.append("vw");
  
  std::ofstream ofs( vw_file.c_str(), std::ios::out );
  
  if ( !ofs )
  {
    std::cerr << "[Visual_Words_Handler]: ERROR: Cannot write the assignments to " << vw_file << std::endl;;
    return;
  }
  
  for( uint32_t i=0; i<nb_descriptors; ++i )
    ofs << assignments[i] << std::endl;
  
  ofs.close();
}

//---------------------------------------------------

void visual_words_handler::assign_and_save_prefix_float( std::string &filename, std::vector< float > &descriptors, uint32_t nb_descriptors )
{
  std::vector< uint32_t > assignments( nb_descriptors, mNbVisualWords+1 );
  assign_visual_words_float( descriptors, nb_descriptors, assignments );
  
  std::string vw_file( filename );
  if( mNbVisualWords == 1000000 )
    vw_file.append("1M.");
  else if( mNbVisualWords == 100000 )
    vw_file.append("100k.");
  if( mNbTrees == 256 )
    vw_file.append("256.");
  else if( mNbTrees == 128 )
    vw_file.append("128.");
  else if( mNbTrees == 64 )
    vw_file.append("64.");
  else if( mNbTrees == 16 )
    vw_file.append("16.");
  else if( mNbTrees == 8 )
    vw_file.append("8.");
  else if( mNbTrees == 4 )
    vw_file.append("4.");
  else if( mNbTrees == 2 )
    vw_file.append("2.");
  vw_file.append("vw");
  
  std::ofstream ofs( vw_file.c_str(), std::ios::out );
  
  if ( !ofs )
  {
    std::cerr << "[Visual_Words_Handler]: ERROR: Cannot write the assignments to " << vw_file << std::endl;;
    return;
  }
  
  for( uint32_t i=0; i<nb_descriptors; ++i )
    ofs << assignments[i] << std::endl;
  
  ofs.close();
}


//---------------------------------------------------

void visual_words_handler::resize( uint32_t nb_descriptors )
{
  if( mMaxDescriptors < (size_t) nb_descriptors )
  {
    mMaxDescriptors = nb_descriptors;
	
    // resize
    if( mMethod == 2 )
    {
      mFlannAssignments.free();
      mFlannAssignments.data = new int[ mMaxDescriptors ];
      mFlannAssignments.rows = mMaxDescriptors;
      mFlannAssignments.cols = 1;
      
      mFlannAssignmentsKNN.free();
      mFlannAssignmentsKNN.data = new int[ mNbNearestNeighbors * mMaxDescriptors ];
      mFlannAssignmentsKNN.rows = mMaxDescriptors;
      mFlannAssignmentsKNN.cols = mNbNearestNeighbors;
      
      mFlannDistances.free();
      mFlannDistances.data = new float[ mMaxDescriptors ];
      mFlannDistances.rows = mMaxDescriptors;
      mFlannDistances.cols = 1;
      
      mFlannDistancesKNN.free();
      mFlannDistancesKNN.data = new float[ mNbNearestNeighbors * mMaxDescriptors ];
      mFlannDistancesKNN.rows = mMaxDescriptors;
      mFlannDistancesKNN.cols = mNbNearestNeighbors;

      mFlannFeatures.free();
      mFlannFeatures.data = new float[128*mMaxDescriptors];
      mFlannFeatures.rows = mMaxDescriptors;
    }
  }
  mFlannFeatures.rows = nb_descriptors;
}

//---------------------------------------------------

void visual_words_handler::initialize( )
{
    
  // resize

  if( mMethod == 2 )
  {
    if( mFlannAssignments.data != 0 )
      mFlannAssignments.free();
    mFlannAssignments.data = new int[ mMaxDescriptors ];
    mFlannAssignments.rows = mMaxDescriptors;
    mFlannAssignments.cols = 1;
    
    if( mFlannAssignmentsKNN.data != 0 )
      mFlannAssignmentsKNN.free();
    mFlannAssignmentsKNN.data = new int[ mNbNearestNeighbors * mMaxDescriptors ];
    mFlannAssignmentsKNN.rows = mMaxDescriptors;
    mFlannAssignmentsKNN.cols = mNbNearestNeighbors;
    
    if( mFlannDistances.data != 0 )
      mFlannDistances.free();
    mFlannDistances.data = new float[ mMaxDescriptors ];
    mFlannDistances.rows = mMaxDescriptors;
    mFlannDistances.cols = 1;
    
    if( mFlannDistancesKNN.data != 0 )
      mFlannDistancesKNN.free();
    mFlannDistancesKNN.data = new float[ mNbNearestNeighbors * mMaxDescriptors ];
    mFlannDistancesKNN.rows = mMaxDescriptors;
    mFlannDistancesKNN.cols = mNbNearestNeighbors;

    mFlannFeatures.free();
    mFlannFeatures.data = new float[size_t(128)*mMaxDescriptors];
    mFlannFeatures.rows = mMaxDescriptors;
    mFlannFeatures.cols = 128;
  }

}

//---------------------------------------------------
