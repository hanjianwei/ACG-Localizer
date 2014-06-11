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


#ifndef VISUAL_WORDS_HANDLER_HH
#define VISUAL_WORDS_HANDLER_HH

/**
 *    Class to assign features to visual words using FLANN.
 *    Can also be used for simple nearest neighbor search.
 *  
 *  author : Torsten Sattler (tsattler@cs.rwth-aachen.de)
 *  date : 10-03-2012
**/ 


#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <stdint.h>
#include <cmath>
#include <time.h>

#include <flann/flann.hpp>


class visual_words_handler
{
  public:
    //! constructor
    visual_words_handler( );
    
    //! destructor
    ~visual_words_handler( );
    
    //! set the method ("flann" or "linear") to use. Per default, flann 
    void set_method( const std::string &method );
    
    // start : functions to set the parameters
    //! specifiy the number of trees to use. Has to be called before computing a new index for flann.
    void set_nb_trees( int nb_trees );
    
    
    //! specifies the number of paths to check
    void set_nb_paths( int nb_path );
    
    //! specifies the type of the flann index (call before creating a new one): "hkmeans", "randomkd" (default), "auto".
    void set_flann_type( const std::string &type );
    
    //! set the number of visual words
    void set_nb_visual_words( uint32_t nb_vw );
    
    //! set the branching of the hierachical k-means tree (flann). Has to be called before computing a new index for flann.
    void set_branching( int nb_branches );
    
    //! set the target precision for autotuned search structure (flann). Has to be called before computing a new index for flann.
    void set_target_precision( float precision );
    
    //! set the build weight for autotuned search structure (flann). Has to be called before computing a new index for flann.
    void set_build_weight( float weight );
    
    //! set the memory weight for autotuned search structure (flann). Has to be called before computing a new index for flann.
    void set_memory_weight( float weight );
    
    //! set the relative sample size for autotuned search structure (flann). Has to be called before computing a new index for flann.
    void set_sample_size( float rel_size );
    
    /**
     * set the number of nearest neighbors to search for (per default 2) when doing k-nn search with flann. Has to be called before computing a new index for flann.
     * Similarly, the function sets the number of nearest neighbors to find using only set of dimensions when doing the visual word assignments.
    **/
    void set_nb_nearest_neighbors( int nb_neighbors );
    
   
    // end : functions to set the parameters
    
    //! set the cluster centers
    void set_cluster_centers( std::vector< float > &cluster_centers );
    void set_cluster_centers( std::string &cluster_file );
    
    //! get the cluster centers
    std::vector<float> get_cluster_centers_flann();
    void get_cluster_centers_flann( std::vector< float > &centers_ );
    
    /**
     * load the tree(s) for flann assignments. If neccessary (and specified), the tree is created. Returns true if the tree was loaded or created
    **/
    bool load_trees_flann( std::string &cluster_file, std::string &tree_file );
    
    //! (re-)builds the flann index
    void rebuild_flann_index();
    
    //! creates a new search index for flann, returns false if index could not be constructed
    bool create_flann_search_index( std::string &cluster_file );
    bool create_flann_search_index( std::vector< float > &cluster_centers );
    
    //! if using a vocabulary tree (flann & hkmeans), get the cluster center ids of the parents at level L for all leaves  
    void get_parents_at_level_L( int L, int* &ids );

    //! assign visual words using the method defined by set_method and stores them in assignments. Returns false if assignments could not be computed
    bool assign_visual_words_uchar( std::vector< unsigned char > &descriptors, uint32_t nb_descriptors, std::vector< uint32_t > &assignments );
    bool assign_visual_words_ucharv( std::vector< unsigned char* > &descriptors, uint32_t nb_descriptors, std::vector< uint32_t > &assignments );
    bool assign_visual_words_float( std::vector< float > &descriptors, uint32_t nb_descriptors, std::vector< uint32_t > &assignments );
    
    /** 
     * assign visual words using the method defined by set_method and stores them in assignments. Returns false if assignments could not be computed.
     * Also returns the squared distances to the cluster centers.
    **/
    bool assign_visual_words_uchar( std::vector< unsigned char > &descriptors, uint32_t nb_descriptors, std::vector< uint32_t > &assignments, std::vector< float > &distances );
    
    /**
     * do a k-nn search using flann. Works only if a flann-index has been created previously. Returns true if assignments could be computed, false otherwise.
     * Both the ids of the nearest neighbors as well as the squared Euclidean distances to the neighbors are stored in assignments respectively distances, where 
     * the first k entries belong to the first query point (k being set by set_nb_nearest_neighbors).
     **/
    bool k_nn_search_flann_uchar( std::vector< unsigned char > &descriptors, uint32_t nb_descriptors, std::vector< uint32_t > &assignments, std::vector< float > &distances );
    bool k_nn_search_flann_ucharv( std::vector< unsigned char* > &descriptors, uint32_t nb_descriptors, std::vector< uint32_t > &assignments, std::vector< float > &distances );
    bool k_nn_search_flann_float( std::vector< float > &descriptors, uint32_t nb_descriptors, std::vector< uint32_t > &assignments, std::vector< float > &distances );
    
    //! load assignments from a file. Needs the number of assignments to expect. Assignments are stored in assignments
    void load_from_file( std::string &filename, uint32_t nb_assignments, std::vector< uint32_t > & assignments );
    
    /** loads assignments from a file. The difference to the function above is that the filename parameter is a prefix of the real filename. 
        The rest of the filename is computed based on the parameters specified.
    **/
    void load_from_file_prefix( std::string &filename, uint32_t nb_assignments, std::vector< uint32_t > & assignments );
    
    //! compute and save the assignments to a file specified by filename
    void assign_and_save_uchar( std::string &filename, std::vector< unsigned char > &descriptors, uint32_t nb_descriptors );
    void assign_and_save_float( std::string &filename, std::vector< float > &descriptors, uint32_t nb_descriptors );
    
    //! same as above, only that the filename is used as a prefix and the rest of the filename depends on the parameters specified.
    void assign_and_save_prefix_uchar( std::string &filename, std::vector< unsigned char > &descriptors, uint32_t nb_descriptors );
    void assign_and_save_prefix_float( std::string &filename, std::vector< float > &descriptors, uint32_t nb_descriptors );

  private:
    //! resize the datastructures to contain feature and assignment information
    void resize( uint32_t nb_descriptors );
    
    //! initialize the datastructures to contain feature and assignment information
    void initialize();
    
    
    //! the method to use: 0 (linear search), 2 (flann, default)
    int mMethod;
    
	//! the matrix containing the cluster centers (for flann)
	flann::Matrix< float > mClusterCentersFlann;
	
	//! the search index for flann
	flann::Index< flann::L2< float > > *mFlannIndex;
	
	//! FLANN:: datastructure for feature data
	flann::Matrix< float > mFlannFeatures;
	
    //! FLANN: datastructure for distances (k-nn search)
	flann::Matrix< float > mFlannDistancesKNN;
	
	//! FLANN: datastructure for distances (visual words)
	flann::Matrix< float > mFlannDistances;
    
    //! FLANN: datastructure for assignments (visual words)
    flann::Matrix< int > mFlannAssignments;
    
    //! FLANN: datastructure for assignments (k-nn search)
    flann::Matrix< int > mFlannAssignmentsKNN;
    
    //! number of trees (default 8)
    int mNbTrees;
    
    //! number of paths (default 1)
    int mNbPath;
    
    //! branching of the hierarchical k-means tree (default: 10)
    int mBranching;
    
    //! target precision for FLANNs AutoTunedIndex (default: 0.9)
    float mTargetPrecision;
    
    //! Build Weight for FLANNs AutoTunedIndex (default: 0.1)
    float mBuildWeight;
    
    //! Memory Weight for FLANNs AutoTunedIndex (default: 0)
    float mMemoryWeight;
    
    //! Sample Fraction for FLANNs AutoTunedIndex (default: 0.1)
    float mSampleFraction;
    
    //! The number of nearest neighbors to search for (default 2) by flann for k-nn search
    int mNbNearestNeighbors;
    
    //! flann index type
    int mFlannIndexType; // 0 : autotune, 1: hkmeans, 2: randomkd (default)
    
    //! number of visual words (default 100k)
    uint32_t mNbVisualWords; 
    
    //! current maximal number of descriptors the assignments structures can hold. Default: 50000
    size_t mMaxDescriptors;
      
};


#endif