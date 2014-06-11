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

#ifndef RANSAC_HH
#define RANSAC_HH

#include <iostream>
#include <vector>
#include <set>
#include <utility>
#include <cmath>
#include <algorithm>
#include <assert.h>
#include <limits>
#include <stdint.h>
#include "math/projmatrix.hh"
#include "solver/solverbase.hh"
#include "solver/solverproj.hh"
#include "math/pseudorandomnrgen.hh"
#include "timer.hh"
#include <ctime>


/**
 *    Implementation of the WaldSac (SPRT-RANSAC) variant of RANSAC.
 *    For details please see the relevant papers:
 *      Chum, O. and Matas. J. : Optimal Randomized RANSAC, PAMI, August 2008
 *    and
 *      Matas, J. and Chum, O.: Randomized RANSAC with Sequential Probability Ratio Test, ICCV 2005
 *    
 *    If you are using this implementation, you should cite one of these papers.
 *    
 *    The implementation includes also the local optimization proposed in
 *      Chum, O., Matas, J. and Obdrzalek, S.: Enhancing RANSAC by Generalized Model Optimization, ACCV 2004
 *    and 
 *      Chum, O., Matas, J. and Kittler, J.: Locally Optimized RANSAC, DAGM 2003
 *
 *  author : Torsten Sattler (tsattler@cs.rwth-aachen.de)
 *  date : 10-02-2011
**/ 

enum ransac_variant{
  SPRT_LO_RANSAC = 2
};

enum ransac_computation_type{
  P6pt = 4
};


/*
 * Implementation of SPRT-LO-RANSAC
 *
 * Author: Torsten Sattler <tsattler@cs.rwth-aachen.de>
 * Date: 10-03-2012
 */

class RANSAC
{
  public:
    RANSAC();
    
    ~RANSAC();
    
    /**
     * standard RANSAC with iterative refinement
     * parameters:
     *   c1 - vector containing the (x,y) coordinates of the found image points in the first image
     *   c2 - vector containing the (x,y) coordinates (or (x,y,z)-coordinates in case of 2D<->3D correspondences
     *        of the found image points in the second image (or the 3D point cloud)
     *   nb_correspondences - number of found correspondences, for 2D<->2D correspondences 
     *                        c1.size() = c2.size() = 2*nb_correspondences should hold
     *                        (c1.size() = 2*nb_correspondences and c2.size() = 3*nb_correspondences for 
     *                        2D<->3D correspondences).
     *   min_inlier_ratio - assumed minimal inlier-ratio
    **/
    void apply_RANSAC( const std::vector< float > &c1, const std::vector< float > &c2, uint32_t nb_correspondences, const float min_inlier_ratio = 0.2f );
    
    //! minimal size of the inlier set of an accepted solution
    void set_minimal_consensus_size( const uint32_t );
    
    //! error threshold for inlier<->outlier classification
    static double error;
    
    //! the type of spatial verification to be performerd
    static ransac_computation_type computation_type;
    
    //! the number of Local Optimization steps LO-RANSAC should take 
    static uint32_t nb_lo_steps;
    
    //! print information, measure time
    static bool silent, measure_time;
    
    //! stop RANSAC after finite amount of time?
    static double max_time; // maximal time after which RANSAC is stopped
    static bool stop_after_n_secs;
    
    //! the maximum number of samples taken in the local optimization step. 
    static uint32_t max_number_of_LO_samples;
    
    //! The t_M variable for the SPRT, specifying the cost of generating a hypothesis relative to evaluating a correspondence. Set to -1 (default) for an automatic, computation-type dependent assignment.
    static float t_M;
    
    //! get the projection matrix computed by RANSAC
    ProjMatrix& get_projection_matrix()
    {
      return projection_matrix;
    }
    
    //! set the projection matrix
    void set_projection_matrix( ProjMatrix& p )
    {
      projection_matrix = p;
    }
    
    /**
	 * Get the indices of the computed inliers. The i-th entry of the returned vector 
	 * is the index to the i-th inlier in c1, respectively c2.
	**/
    std::vector< uint32_t >& get_inliers()
    {
      return inlier_correspondences;
    }
    
    //! get the number of inliers computed by RANSAC
    uint32_t get_number_of_inliers()
    {
      return inlier_correspondences.size();
    }
    
      
    //! the RANSAC type used. Per default set to SPRT_LO_RANSAC
    static ransac_variant inner_RANSAC_type;
    
    //! get the number of steps needed by RANSAC (not including local optimization steps)
    uint32_t get_nb_ransac_steps()
    {
      return taken_samples;
    }
    
    //! get the time needed to perform RANSAC (in seconds)
    double get_elapsed_time()
    {
      return elapsed_time;
    }
    

    // variables for estimating the value of t_M (SPRT) through experiments
//     static double avrg_model_generation_time, nb_generated_models, avrg_model_evaluation_time, nb_evaluations;
       

  private: 
	//! Implementation of WaldSaC (SPRT Ransac), see O. Chum and J. Matas. PAMI, Vol. 30, No. 8. 2008
    void unified_SPRT_LO_RANSAC( const std::vector< float > &c1, const std::vector< float > &c2, const float min_inlier_ratio = 0.2f );
    
    
    void compute_final_correspondences( const std::vector< float > &c1, const std::vector< float > &c2, std::vector< uint32_t > &indices );  
   
    void initialize();
    
    void clear_solver();
    
    void add_correspondence( const std::vector< float > &c1, const std::vector< float > &c2, uint32_t index );
    
    bool evaluate_correspondece( const std::vector< float > &c1, const std::vector< float > &c2, uint32_t index );
    
    bool solve_system();
    
    void store_hypothesis( const bool get_epipoles = false );
    
    void set_hypothesis();
    
    // compute the maximal number of steps needed by RANSAC 
    uint32_t get_max_ransac_steps( float inlier_ratio );
    uint32_t get_max_ransac_steps( float inlier_ratio, int nb_samples_to_take );
    
    // compute the number of steps needed by SPRT RANSAC for the current test (this includes computing all the h_i)
    // the arguments are the current inlier ratio and the number of the current test
    uint32_t SPRT_get_max_sprt_ransac_steps( float epsilon, int l );
    uint32_t SPRT_get_max_sprt_ransac_steps( float epsilon, int l, int nb_samples_to_take );
    
    // SPRT: compute the decision threshold A, see Chum, Matas. Optimal Randomized RANSAC. PAMI, Vol. 30, No. 8. 2008
    float sprt_compute_A( float eps, float delta );
    
    uint32_t size_inlier_set;
    uint32_t taken_samples;
    double elapsed_time;
    double initialization_time;
    float inlier_ratio;
    float SPRT_m_s; // mean number of models per sample
    uint32_t random_sample;
    
    // the number of LO samples is the maximal number of correspondences choosen in the LO-step (not necessarily minimal)
    uint32_t nr_ransac_steps, number_preliminary_matches, number_of_samples, minimal_consensus_set_size;
    
    //! solvers for different geometric models
    Util::CorrSolver::SolverProj solverP6pt;
    
    //! timer for stopping calculation time of algorithms
    Timer RANSAC_timer;
    
    //! the estimated hypothetes
    ProjMatrix projection_matrix;
    Util::CorrSolver::Vector3D cam_orientation_tmp, cam_position_tmp;
    
    //! different kinds of correspondence sets computed through RANSAC
    std::vector< uint32_t > correspondence_indices;
    std::vector< uint32_t > inlier_correspondences;
    std::vector< uint32_t > inlier;
    std::vector< uint32_t > outlier;
    
    //! datastructures for storing SPRT-related values
    std::vector< float > epsilon_i, delta_i, A_i, h_i;
    std::vector< uint32_t > k_i;
    size_t nb_SPRT_tests;
    
	//! Random Number Generator
    Util::Math::PRNG_u32i random_number_gen;
    
    // multiplicator for look-up in correspondences array (2 for RANSAC, 3 for SCRAMSAC)
    uint32_t index_multiplicator;
    
    double internat_error;
};
#endif
