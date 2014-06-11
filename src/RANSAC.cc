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


#include "RANSAC.hh"
#include <utility>

// small typedef for the maximal value of uint32_t, used to 
// avoid having to define __STDC_LIMIT_MACROS 
#define MAX_VAL_UINT32_T 4294967295

// to speed up the computation, we avoid to recompute log(0.05f) but store it
#define LOG_5_PER -2.99573

// initialize static variables
ransac_computation_type RANSAC::computation_type = P6pt;
double RANSAC::error = 1.0;
bool RANSAC::silent = false;
bool RANSAC::measure_time = true;
double RANSAC::max_time = 10.0; // 10 seconds
bool RANSAC::stop_after_n_secs = false;
uint32_t RANSAC::max_number_of_LO_samples = 0;
float RANSAC::t_M = -1.0f;
uint32_t RANSAC::nb_lo_steps = 0;

ransac_variant RANSAC::inner_RANSAC_type = SPRT_LO_RANSAC;

// SPRT: Lookup table for starting values for h_i
// SPRT: depending on the values of epsilon (first dimension), 
// SPRT: epsilon_i (second dimension) and delta_i (third dimension). 
// SPRT: first the used values for epsilon (same for epsilon_i)
const int nb_epsilon_values = 8;
const float SPRT_eps_val[8] = { 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f };
// SPRT: values for delta_i. NOTE THAT WE ASSUME delta_i TO BE AT MOST 0.42.
// SPRT: For delta_i > 0.42 the computation of the h_i might lead to the wrong 
// SPRT: solution, i.e. h_i = 0 instead of h_i >1 
const int nb_delta_values = 10;
const float SPRT_delta_val[10] = {  0.01f, 0.02f, 0.07f, 0.12f, 0.17f, 0.22f, 0.27f, 0.32f, 0.37f, 0.42f };
// SPRT: Notice that only the case of epsilon >= epsilon_i > delta_i is meaningful.
// SPRT: Now follow the values for h_i. 0 indicates a non-meaningful parameter combination.
const float SPRT_h_i_val[8][8][10] = { 1.0f, 1.0f, 1.0f, 1.0015f, 1.0249f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
                       0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
                       0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
                       0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
                       0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
                       0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
                       0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
                       0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
                       1.6642f, 1.7299f, 2.1589f, 3.062f, 6.8894f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
                       1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0007f, 1.1308f, 0.0f, 0.0f, 0.0f, 
                       0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
                       0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
                       0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
                       0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
                       0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
                       0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
                       2.3958f, 2.511f, 3.3098f, 5.0339f, 12.3835f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
                       1.4658f, 1.4975f, 1.6723f, 1.9206f, 2.3428f, 3.2742f, 7.2705f, 0.0f, 0.0f, 0.0f, 
                       1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0001f, 1.1615f, 0.0f, 
                       0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
                       0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
                       0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
                       0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
                       0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
                       3.2527f, 3.4146f, 4.5761f, 7.1345f, 18.1024f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
                       1.9981f, 2.0543f, 2.3846f, 2.8683f, 3.6984f, 5.5358f, 13.4329f, 0.0f, 0.0f, 0.0f, 
                       1.3779f, 1.3972f, 1.4959f, 1.6163f, 1.7805f, 2.0283f, 2.4583f, 3.4147f, 7.5328f, 0.0f, 
                       1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0001f, 1.0001f, 
                       0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
                       0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
                       0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
                       0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
                       4.2999f, 4.515f, 6.0786f, 9.5662f, 24.5884f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
                       2.6433f, 2.7221f, 3.2052f, 3.932f, 5.1915f, 7.9902f, 20.0375f, 0.0f, 0.0f, 0.0f, 
                       1.8283f, 1.863f, 2.0522f, 2.2912f, 2.6209f, 3.121f, 3.9915f, 5.93f, 14.282f, 0.0f, 
                       1.3367f, 1.35f, 1.4156f, 1.4894f, 1.5799f, 1.6982f, 1.8634f, 2.1155f, 2.5557f, 3.5371f, 
                       1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 
                       0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
                       0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
                       0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
                       5.65f, 5.9326f, 7.9949f, 12.6205f, 32.6092f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
                       3.4734f, 3.5781f, 4.2326f, 5.2359f, 6.9898f, 10.9034f, 27.7798f, 0.0f, 0.0f, 0.0f, 
                       2.404f, 2.453f, 2.7335f, 3.0995f, 3.6105f, 4.3906f, 5.7528f, 8.7914f, 21.8925f, 0.0f, 
                       1.7615f, 1.7858f, 1.9137f, 2.0634f, 2.2497f, 2.4947f, 2.8383f, 3.364f, 4.2829f, 6.3331f, 
                       1.3251f, 1.3352f, 1.3838f, 1.4357f, 1.4956f, 1.5679f, 1.659f, 1.7796f, 1.9491f, 2.2087f, 
                       1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 
                       0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
                       0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
                       7.5527f, 7.9306f, 10.6887f, 16.8849f, 43.7002f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
                       4.6432f, 4.7833f, 5.6642f, 7.0274f, 9.4258f, 14.7974f, 38.0023f, 0.0f, 0.0f, 0.0f, 
                       3.2139f, 3.2803f, 3.6693f, 4.1888f, 4.9231f, 6.0516f, 8.0298f, 12.4514f, 31.5338f, 0.0f, 
                       2.356f, 2.3911f, 2.5854f, 2.8216f, 3.1203f, 3.5167f, 4.0755f, 4.9328f, 6.4339f, 9.7866f, 
                       1.7753f, 1.7941f, 1.8911f, 2.0f, 2.1278f, 2.2836f, 2.481f, 2.7431f, 3.1124f, 3.6786f, 
                       1.3458f, 1.3541f, 1.3935f, 1.4344f, 1.4796f, 1.5315f, 1.593f, 1.6682f, 1.7637f, 1.8903f, 
                       1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 
                       0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
                       10.8055f, 11.3461f, 15.2922f, 24.1588f, 62.5457f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
                       6.6429f, 6.8433f, 8.1047f, 10.0615f, 13.5147f, 21.2668f, 54.8035f, 0.0f, 0.0f, 0.0f, 
                       4.598f, 4.6932f, 5.2538f, 6.0104f, 7.0894f, 8.758f, 11.6945f, 18.2739f, 46.7041f, 0.0f, 
                       3.3708f, 3.4216f, 3.7094f, 4.0683f, 4.5298f, 5.1483f, 6.0258f, 7.3778f, 9.7513f, 15.0601f, 
                       2.5408f, 2.5694f, 2.726f, 2.9098f, 3.1304f, 3.4028f, 3.751f, 4.2155f, 4.8725f, 5.8821f, 
                       1.9284f, 1.9444f, 2.0277f, 2.1195f, 2.2239f, 2.3455f, 2.4908f, 2.6696f, 2.8972f, 3.2f, 
                       1.4386f, 1.4461f, 1.4822f, 1.5192f, 1.5591f, 1.6033f, 1.6535f, 1.7119f, 1.7814f, 1.8666f, 
                       1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f
                        };

RANSAC::RANSAC()
{
  initialize();
  nb_SPRT_tests = 100;
  epsilon_i.resize(nb_SPRT_tests,0.0f);
  delta_i.resize(nb_SPRT_tests,0.0f);
  A_i.resize(nb_SPRT_tests,0.0f);
  h_i.resize(nb_SPRT_tests,0.0f);
  k_i.resize(nb_SPRT_tests,0);
}

//-----------------------------------


RANSAC::~RANSAC()
{
  initialize();
  epsilon_i.clear();
  delta_i.clear();
  A_i.clear();
  h_i.clear();
  k_i.clear();
}

//-----------------------------------


void RANSAC::apply_RANSAC( const std::vector< float > &c1, const std::vector< float > &c2, uint32_t nb_correspondences, const float min_inlier_ratio )
{
  //initialize counting variables and various parameters
  initialize();
  
  //start measuring time
  if( measure_time )
  {
    RANSAC_timer.Init();
    RANSAC_timer.Start();
  }
  
  // use all correspondences
  correspondence_indices.resize(nb_correspondences);
  for( uint32_t i=0; i<nb_correspondences; ++i )
    correspondence_indices[i] = i;
  //apply RANSAC procedure
  
  if( !silent )
    std::cout << "[RANSAC] Applying RANSAC on " << nb_correspondences << " correspondences " << std::endl;
  
  switch( inner_RANSAC_type )
  {
    case SPRT_LO_RANSAC:
    {
      unified_SPRT_LO_RANSAC( c1, c2, min_inlier_ratio );
      break;
    }
    
    default:
      break;
  }
  
  if( measure_time )
  {
    RANSAC_timer.Stop();
    elapsed_time = RANSAC_timer.GetElapsedTime();
  }
  
  if( !silent )
  {
    std::cout << "[RANSAC] Applying RANSAC took " << RANSAC_timer.GetElapsedTimeAsString() << std::endl;
  }
  
  // explicitly compute the correspondences of the best hypothesis found in RANSAC
  compute_final_correspondences( c1, c2, correspondence_indices );  
}

//-----------------------------------


void RANSAC::set_minimal_consensus_size( const uint32_t size)
{
  minimal_consensus_set_size = size;
}

//-----------------------------------


void RANSAC::unified_SPRT_LO_RANSAC( const std::vector< float > &c1, const std::vector< float > &c2, const float min_inlier_ratio )
{
  //initialization
  size_inlier_set = 0;
  uint32_t nb_correspondences = (uint32_t) correspondence_indices.size();

  if( nb_correspondences >= number_of_samples )
  {
    //some more initialization
   
    uint32_t inlier_found = 0;
	    
	// compute the inlier ratio and the maximum number of steps that RANSAC has to take 
    inlier_ratio = std::max( min_inlier_ratio, ((float) number_of_samples) / ((float) nb_correspondences) );
    uint32_t max_steps = get_max_ransac_steps( inlier_ratio );
	
    if(!silent) 
	  std::cout << "[RANSAC] initial inlier ratio : " << inlier_ratio << " resulting in at most " << max_steps << " steps " << std::endl;

    size_inlier_set = number_of_samples-1;
    
    uint32_t new_found_inlier = 0;
    
    std::vector<uint32_t> randomly_choosen_corr_indices, LO_randomly_choosen_corr_indices;
    randomly_choosen_corr_indices.resize( number_of_samples, 0 );
    LO_randomly_choosen_corr_indices.resize( max_number_of_LO_samples, 0 );
    
    //actual SPRT RANSAC algorithm
    
    // initialize the 0-th SPRT test
    epsilon_i[0] = inlier_ratio;
    delta_i[0] = 0.01f;
    
    A_i[0] = sprt_compute_A( epsilon_i[0], delta_i[0] );
    int current_test = 0;
    int old_test = -1;
    
    // handle timing. Need to split into two parts due to nummerics
    double elapsed_time_total = 0.0;
    double elapsed_time_sec = 0.0;
    double old_time = 0.0;
    double sum_time = 0.0;
    int old_mins = 0;
    int elapsed_mins = 0;
    uint32_t nb_LO_samples;
    uint32_t counter = 0;
    Timer time_out_timer;
    time_out_timer.Init();
    time_out_timer.Start();
    
    float eval_1 = delta_i[0] / epsilon_i[0];
    float eval_0 = (1.0f-delta_i[0]) / (1.0f-epsilon_i[0]);
    float lambda = 1.0f;
    float delta_hat = delta_i[0];
    float old_ratio = inlier_ratio;
    float nb_rejected = 1.0f;
    bool bad_model = false;
	
    taken_samples = 0;
    k_i[0] = 0;
	
    while( true )
    {
    
      for( ; taken_samples < max_steps/*0.01*/;)
      {
		if( stop_after_n_secs )
		{
		  elapsed_time_sec = time_out_timer.GetElapsedTime();
		  time_out_timer.Stop();
		  time_out_timer.Restart();
		  if( elapsed_time_sec >= 1.0 )
		  {
			elapsed_time_total += elapsed_time_sec;
			elapsed_time_sec = 0.0;
			time_out_timer.Stop();
			time_out_timer.Start();
			elapsed_mins = int(floor(elapsed_time_total)) / 60;
			if( elapsed_mins > old_mins && elapsed_mins > 0 )
			{
			  std::cout << "[RANSAC] elapsed time so far " << elapsed_time_total << " s " << std::endl;
			  old_mins = elapsed_mins;
			}
		  }
		  old_time = sum_time;
		  sum_time = elapsed_time_sec + elapsed_time_total;
		  if( sum_time >= max_time )
		  {
			std::cout << "[RANSAC] Warning: RANSAC reached time limit of " << max_time << " s and was stopped" << std::endl;
			break;
		  }
		  if( old_time == sum_time && old_time > 0.0 )
		  {
			std::cout << "[RANSAC] Warning: Time did not change value (possibly due to nummerical reasons), so we abort after " << old_time << " s" << std::endl;
			break;
		  }
		}
		++taken_samples;
		k_i[current_test] += 1;
		
		if( taken_samples == 0 )
		{
		  std::cerr << "[RANSAC] Error: The number of samples taken exceeds " << --taken_samples << " ( 2^64-1 ) which is the maximal number representable by a uint32_t." << std::endl;
		  std::cerr << "[RANSAC] Error: Therefore, we stop searching for a better model here. (Maybe you want to use SCRAMSAC, if applicable, to get rid of outlier!)" << std::endl;
		  break;
		}
	
		//take a random sample from the set of correspondences
		random_number_gen.generate_pseudorandom_numbers_unique( (uint32_t) 0, (uint32_t) ( nb_correspondences - 1 ), number_of_samples, randomly_choosen_corr_indices );
		
			  
		//add the correspondences
		clear_solver();
		for(std::vector<uint32_t>::iterator itc = randomly_choosen_corr_indices.begin(); itc !=  randomly_choosen_corr_indices.end(); ++itc)
		{
		  add_correspondence(c1,c2,correspondence_indices[*itc]);
		}
		
			  
		//compute hypothesis
		if(!solve_system())
		  continue;
		
		//compute number of inlier to the hypothesis, evaluate the current SPRT test
		inlier_found = 0;
		lambda = 1.0f;
		bad_model = false;
		for(std::vector< uint32_t >::const_iterator it = correspondence_indices.begin() ; it != correspondence_indices.end(); ++it )
		{
		  if(evaluate_correspondece(c1,c2,*it))
		  {
			lambda *= eval_1;
			++inlier_found;
		  }
		  else
		  {
			lambda *= eval_0;
		  }
		  if( lambda > A_i[current_test] )
		  {
			bad_model = true;
			break;
		  }
		}
	
		if( bad_model )
		{
		  // check if we have to design a new test
		  nb_rejected += 1.0f;
		  delta_hat = delta_hat *(nb_rejected-1.0f) / nb_rejected + float(inlier_found)/(float(nb_correspondences) * nb_rejected);
		  
		  if( fabs( delta_hat - delta_i[current_test] ) > 0.05f )
		  {
			++current_test;
			k_i[current_test] = 0;
			epsilon_i[current_test] = epsilon_i[current_test-1];
			delta_i[current_test] = delta_hat;
			A_i[current_test] = sprt_compute_A( epsilon_i[current_test], delta_hat );
		  }
		  
		  continue;
		}
		

		//compare found inliers to the biggest set of correspondences found so far
		if(inlier_found > size_inlier_set )
		{
		  //store hypothesis and update inlier ratio
		  store_hypothesis();
		  
		  
		  // compute inlier
		  inlier.clear();
		  for( std::vector< uint32_t >::const_iterator it = correspondence_indices.begin() ; it != correspondence_indices.end(); ++it )
		  {
			if(evaluate_correspondece(c1,c2,*it))
			  inlier.push_back( *it );  
		  }

		  //do Local Optimization (LO)-steps (if possible!)
		  nb_LO_samples = std::max( number_of_samples, std::min( inlier_found / 2, max_number_of_LO_samples ));
		  for( uint32_t lo_steps = 0; lo_steps < nb_lo_steps; ++lo_steps )
		  {
			random_number_gen.generate_pseudorandom_numbers_unique( (uint32_t) 0, (uint32_t) (inlier.size() - 1), nb_LO_samples, LO_randomly_choosen_corr_indices );
			//generate hypothesis
			//add the correspondences
			clear_solver();
			
			for( counter = 0; counter < nb_LO_samples; ++counter )
			{
			  add_correspondence(c1,c2,inlier[LO_randomly_choosen_corr_indices[counter]]);
			}
			
			//compute hypothesis
			if(!solve_system())
			  continue;
			
			//compute inlier to the hypothesis
			new_found_inlier = 0;
			for(std::vector< uint32_t >::const_iterator it = correspondence_indices.begin(); it != correspondence_indices.end(); ++it )
			{
			  if(evaluate_correspondece(c1,c2,*it))
			++new_found_inlier;
			}
			
			//update found model if new best hypothesis
			if( new_found_inlier > inlier_found )
			{
			  inlier_found = new_found_inlier;
			  store_hypothesis();
			}
		  }
		  
		  old_ratio = inlier_ratio;
		  inlier_ratio = std::max( inlier_ratio, (float) inlier_found / (float) nb_correspondences );
		  max_steps = get_max_ransac_steps( inlier_ratio );
		  size_inlier_set = inlier_found;
		  
	// 	  std::cout << projection_matrix << std::endl << std::endl;
		  
		  // design a new test if needed
		  // the check must be done to avoid designing a test for an inlier ratio below the specified minimal inlier ratio
		  if( old_ratio < inlier_ratio )
		  {
			++current_test;
			k_i[current_test] = 0;
			epsilon_i[current_test] = inlier_ratio;
			delta_i[current_test] = delta_hat;
			A_i[current_test] = sprt_compute_A(inlier_ratio, delta_hat);
		  }
		}
      }
      
      // adjust the number of steps SPRT RANSAC has to take 
      if( old_test != current_test )
      {
		old_test = current_test;
		max_steps = SPRT_get_max_sprt_ransac_steps( inlier_ratio, current_test );
      }
      else
		break;
    }
    
      
    if( !silent )
      std::cout << "[RANSAC] SPRT-LO-RANSAC took " << taken_samples << " samples using " << current_test+1 << " SPRTs, found " << size_inlier_set << " inlier ( " << inlier_ratio << " % ) " << std::endl;
  }
}


//-----------------------------------

void RANSAC::compute_final_correspondences( const std::vector< float > &c1, const std::vector< float > &c2, std::vector< uint32_t > &indices )
{
  if( size_inlier_set < minimal_consensus_set_size )
  {
    inlier_ratio = (float) size_inlier_set / (float) indices.size();
    return;
  }
  
//   std::cout << " nb inlier : " << size_inlier_set << " ratio : " << inlier_ratio << std::endl;
  
  uint32_t nb_correspondences = (uint32_t) indices.size();
  
  inlier_correspondences.clear();
  
 
  size_inlier_set = 0;
  
 
  set_hypothesis();
  
  
  for(std::vector< uint32_t >::const_iterator it = indices.begin() ; it != indices.end(); ++it )
  {
    if(evaluate_correspondece(c1,c2,*it))
    {
      inlier_correspondences.push_back( *it );
      ++size_inlier_set;
    }
  }

  inlier_ratio = float(size_inlier_set)/float(nb_correspondences);
  
  if( !silent )
  {
    std::cout << "[RANSAC] Percentage Inlier found on all (reduced) correspondences : " << (float) inlier_ratio << std::endl;
    std::cout << "[RANSAC] Inlier found on set of (reduced) correspondences : " << size_inlier_set << std::endl;
  }
}


//-----------------------------------

void RANSAC::initialize()
{
  if( computation_type == P6pt )
  {
    index_multiplicator = 3;
  }
  else
    index_multiplicator = 2;
  //clear solvers
  clear_solver();

  //clear resulting correspondences
  correspondence_indices.clear();
  inlier_correspondences.clear();
  inlier.clear();
  outlier.clear();
  inlier.clear();
  outlier.clear();

  //reset counting variables
  size_inlier_set = 0;
  inlier_ratio = 0.0;
  taken_samples = 0;
  elapsed_time = 0.0;
  initialization_time = 0.0;
	  
	  
  //set parameters
  SPRT_m_s = 1.0f;
  switch( computation_type )
  {
    case P6pt:
    {
      minimal_consensus_set_size = number_of_samples = 6;
      break;
    }
  }
  
  if( max_number_of_LO_samples == 0 )
  {
    switch( computation_type )
    {
      case P6pt:
      {
		max_number_of_LO_samples = 12; // NOTE: This does not have to be an optimal setting
		break;
      }
    }
  }
  
  if( max_number_of_LO_samples == 0 )
  {
    switch( computation_type )
    {
      case P6pt:
      {
		t_M = 200.0f;
		break;
      }
    }
  }
  
  internat_error = error;
  
}


//-----------------------------------

void RANSAC::clear_solver()
{
  switch(computation_type)
  {
    case P6pt:
    {
      solverP6pt.clear();
      break;
    }
    
    default:
    {
      break;
    }
  }
}

//-----------------------------------

void RANSAC::add_correspondence( const std::vector< float > &c1, const std::vector< float > &c2, uint32_t index )
{
//   std::cout << " add " << index << std::endl;
  uint32_t real_index = index * index_multiplicator;
  switch(computation_type)
  {
    case P6pt:
    {
      solverP6pt.addCorrespondence(Util::CorrSolver::Vector2D(c1[2*index], c1[2*index+1]), Util::CorrSolver::Vector3D(c2[real_index], c2[real_index+1], c2[real_index+2]) );
      break;
    }
    
    default:
    {
      break;
    }
  }
}

//-----------------------------------

bool RANSAC::evaluate_correspondece( const std::vector< float > &c1, const std::vector< float > &c2, uint32_t index )
{
  uint32_t real_index = index * index_multiplicator;
  switch(computation_type)
  {
    
    case P6pt:
    {
      Util::CorrSolver::Vector3D p3d( c2[real_index], c2[real_index+1], c2[real_index+2] );
	  
      if( ( (p3d - cam_position_tmp) | cam_orientation_tmp ) >= 0.0 )
		return (solverP6pt.evaluateCorrespondence(Util::CorrSolver::Vector2D(c1[2*index], c1[2*index+1]), p3d ) <= internat_error);
	  
      break;
    }
    
    default:
    {
      return false;
      break;
    }
  }
  return false;
}

//-----------------------------------


bool RANSAC::solve_system()
{
  bool solved = false;
  switch(computation_type)
  {
    case P6pt:
    {
      solved = solverP6pt.computeLinearNew();
      
      if( solved )
		solved = solverP6pt.getPositionAndOrientation( cam_position_tmp, cam_orientation_tmp );

      break;
    }
    
    default:
    {
      solved = false;
      break;
    }
  }
  return solved;
}

//-----------------------------------


void RANSAC::store_hypothesis( const bool get_epipoles )
{  
  switch(computation_type)
  {
    case P6pt:
    {
      solverP6pt.getProjectionMatrix( projection_matrix );
      break;
    }
    
    default:
    {
      break;
    }
  }
}

//-----------------------------------

void RANSAC::set_hypothesis()
{
  switch(computation_type)
  {
    case P6pt:
    {
      solverP6pt.setProjectionMatrix( projection_matrix );
      solverP6pt.getPositionAndOrientation( cam_position_tmp, cam_orientation_tmp );
      break;
    }
    
	default:
    {
      break;
    }
  }
}

//-----------------------------------

uint32_t RANSAC::get_max_ransac_steps( float inlier_ratio )
{
  // get the maximal number of steps RANSAC has to take
  // make sure that it does not exceed the maximal value of uint32_t
  if(inlier_ratio == 1.0)
  {
      return 1;
  }
  float real_ratio = inlier_ratio;
  switch(number_of_samples)
  {
    case 3:
    {
      real_ratio = std::max( real_ratio, 0.000887f );
      return uint32_t(ceil(LOG_5_PER/log(1.0-real_ratio*real_ratio*real_ratio)));
      break;
    }
    
    case 4:
    {
      real_ratio = std::max( real_ratio, 0.005139f );
      return uint32_t(ceil(LOG_5_PER/log(1.0-real_ratio*real_ratio*real_ratio*real_ratio)));
      break;
    }
    
    case 6:
    {
      real_ratio = std::max( real_ratio, 0.029780f );
      return uint32_t(ceil(LOG_5_PER/log(1.0-real_ratio*real_ratio*real_ratio*real_ratio*real_ratio*real_ratio)));
      break;
    }   

    case 8:
    {
      real_ratio = std::max( real_ratio, 0.071687f );
      return uint32_t(ceil(LOG_5_PER/log(1.0-real_ratio*real_ratio*real_ratio*real_ratio*real_ratio*real_ratio*real_ratio*real_ratio)));
      break;
    }

    default:
    {
      return 0;
      break;
    }
  }
}

//-----------------------------------

uint32_t RANSAC::get_max_ransac_steps( float inlier_ratio, int nb_samples_to_take )
{
  // get the maximal number of steps RANSAC has to take
  // make sure that it does not exceed the maximal value of uint32_t
  if(inlier_ratio == 1.0){
      return 1;
  }
  float real_ratio = inlier_ratio;
  switch(nb_samples_to_take)
  {
    case 2:
    {
      real_ratio = std::max( real_ratio, 0.000001f );
      return uint32_t(ceil(LOG_5_PER/log(1.0-real_ratio*real_ratio)));
      break;
    }
    
    case 3:
    {
      real_ratio = std::max( real_ratio, 0.000887f );
      return uint32_t(ceil(LOG_5_PER/log(1.0-real_ratio*real_ratio*real_ratio)));
      break;
    }
    
    case 4:
    {
      real_ratio = std::max( real_ratio, 0.005139f );
      return uint32_t(ceil(LOG_5_PER/log(1.0-real_ratio*real_ratio*real_ratio*real_ratio)));
      break;
    }
    
    case 5:
    {
      real_ratio = std::max( real_ratio, 0.014747f );
      return uint32_t(ceil(LOG_5_PER/log(1.0-real_ratio*real_ratio*real_ratio*real_ratio*real_ratio)));
      break;
    }
    
    case 6:
    {
      real_ratio = std::max( real_ratio, 0.029780f );
      return uint32_t(ceil(LOG_5_PER/log(1.0-real_ratio*real_ratio*real_ratio*real_ratio*real_ratio*real_ratio)));
      break;
    } 
    
    case 7:
    {
      real_ratio = std::max( real_ratio, 0.049197f );
      return uint32_t(ceil(LOG_5_PER/log(1.0-real_ratio*real_ratio*real_ratio*real_ratio*real_ratio*real_ratio*real_ratio)));
      break;
    }  

    default:
    {
      return 0;
      break;
    }
  }
}

//-----------------------------------

uint32_t RANSAC::SPRT_get_max_sprt_ransac_steps( float epsilon, int l )
{
  {
    // levenberg-marquardt variables
    float a,b,c,d,e,f,ff,fff,f_1,A,g,x_old,x_new,h_lm,F_x_new,F_x_old,t,upsilon,mu, tmp_flt1, tmp_flt2, rho;
    int k;
    bool found;
   
    // compute all the h_i
    int idx_eps, idx_eps_i, idx_delta_i;
    Timer h_i_timer;
    h_i_timer.Init();
    h_i_timer.Start();
    for( int i=0; i<l; ++i )
    {
      // do a table look-up to get a approximate starting value for h_i
      
      // first find a similar value for epsilon
      tmp_flt1 = fabs( epsilon - SPRT_eps_val[0] );
      idx_eps = 0;
      for( int j=1; j<nb_epsilon_values; ++j )
      {
		tmp_flt2 = fabs( epsilon - SPRT_eps_val[j] );
		if( tmp_flt2 > tmp_flt1 )
		  break;
		else
		  ++idx_eps;
      }
      
      // now for epsilon_i
      tmp_flt1 = fabs( epsilon_i[i] - SPRT_eps_val[0] );
      idx_eps_i = 0;
      for( int j=1; j<nb_epsilon_values; ++j )
      {
		tmp_flt2 = fabs( epsilon_i[i] - SPRT_eps_val[j] );
		if( tmp_flt2 > tmp_flt1 )
		  break;
		else
		  ++idx_eps_i;
      }
      
      // now for delta_i
      tmp_flt1 = fabs( delta_i[i] - SPRT_delta_val[0] );
      idx_delta_i = 0;
      for( int j=1; j<nb_delta_values; ++j )
      {
		tmp_flt2 = fabs( epsilon_i[i] - SPRT_delta_val[j] );
		if( tmp_flt2 > tmp_flt1 )
		  break;
		else
		  ++idx_delta_i;
      }
      
      x_old = SPRT_h_i_val[idx_eps][idx_eps_i][idx_delta_i];
      if( x_old == 0.0f )
		x_old = 100.0f;
      
      x_new = x_old;
      
      // now do levenberg-marquardt for refinement
      a = delta_i[i] / epsilon_i[i];
      b = (1.0f-delta_i[i]) / (1.0f -epsilon_i[i]);
      c = 1.0f - epsilon;
      d = log(a);
      e = log(b);
      
      k = 0;
      ff = pow(a,x_old);
      fff = pow(b,x_old);
      f_1 = epsilon * d * ff + c * e * fff;
      A = f_1 * f_1;
      F_x_new = F_x_old = epsilon*ff*+c*fff-1.0f;
      g = f_1 * F_x_old;
      
      mu = 1e-6f * A;
      upsilon = 2.0f;
      
      found = fabs(g) <= 1e-6f;
      
      while( !found && k<200 )
      {
		++k;
		h_lm = -g/(A+mu);
		found = fabs(h_lm) <= 1e-6f * (fabs(x_old)+1e-6f );
		x_new = x_old + h_lm;
		F_x_new = epsilon * pow(a,x_new) + c * pow(b, x_new) - 1.0f;
		rho = 2.0f * (F_x_old * F_x_old - F_x_new * F_x_new ) / (h_lm*(mu*h_lm-g));
		if( rho > 0.0f )
		{
		  x_old = x_new;
		  f_1 = epsilon * d * pow(a,x_old) * ff + c * e * pow(b,x_old);
		  A = f_1 * f_1;
		  g = f_1 * F_x_new;
		  F_x_old = F_x_new;
		  found = fabs(g) <= 1e-6f;
		  t = 2.0f* rho - 1.0f;
		  mu *= std::max(1.0f/3.0f, 1.0f - t*t*t );
		  upsilon = 2.0f;
		}
		else
		{
		  mu *= upsilon; upsilon *= 2.0f;
		}
      }
      
      h_i[i] = x_old;
    }
    h_i_timer.Stop(); 
    
  }
  
  // now we compute the probability eta(l-1)
  float eta_l_minus_1 = 0.0f;
  float Pg = 1.0f;
  switch(number_of_samples)
  {
    
    case 3:
    {
      Pg = epsilon * epsilon * epsilon;
      break;
    }
    
    case 4:
    {
      Pg = epsilon * epsilon * epsilon * epsilon;
      break;
    }

    case 6:
    {
      Pg = epsilon * epsilon * epsilon * epsilon * epsilon * epsilon;
      break;
    }
    
    case 8:
    {
      Pg = epsilon * epsilon * epsilon * epsilon * epsilon * epsilon * epsilon * epsilon;
      break;
    }
  }
  for( int i=0; i<l; ++i )
    eta_l_minus_1 += log( 1.0f - Pg * (1.0f - pow( A_i[i], -h_i[i] )) ) * k_i[i];
  
  float numerator = LOG_5_PER - eta_l_minus_1;
  
  if( numerator >= 0.0f )
    return 0;
  
  return uint32_t( ceil( numerator / log(1.0f - Pg * (1.0f - 1.0f/ A_i[l] ) ) ) );
}


//-----------------------------------

uint32_t RANSAC::SPRT_get_max_sprt_ransac_steps( float epsilon, int l, int nb_samples_to_take )
{
  {
    // levenberg-marquardt variables
    float a,b,c,d,e,f,ff,fff,f_1,A,g,x_old,x_new,h_lm,F_x_new,F_x_old,t,upsilon,mu, tmp_flt1, tmp_flt2, rho;
    int k;
    bool found;
   
    // compute all the h_i
    int idx_eps, idx_eps_i, idx_delta_i;
    Timer h_i_timer;
    h_i_timer.Init();
    h_i_timer.Start();
    for( int i=0; i<l; ++i )
    {
      // do a table look-up to get a approximate starting value for h_i
      
      // first find a similar value for epsilon
      tmp_flt1 = fabs( epsilon - SPRT_eps_val[0] );
      idx_eps = 0;
      for( int j=1; j<nb_epsilon_values; ++j )
      {
		tmp_flt2 = fabs( epsilon - SPRT_eps_val[j] );
		if( tmp_flt2 > tmp_flt1 )
		  break;
		else
		  ++idx_eps;
      }
      
      // now for epsilon_i
      tmp_flt1 = fabs( epsilon_i[i] - SPRT_eps_val[0] );
      idx_eps_i = 0;
      for( int j=1; j<nb_epsilon_values; ++j )
      {
		tmp_flt2 = fabs( epsilon_i[i] - SPRT_eps_val[j] );
		if( tmp_flt2 > tmp_flt1 )
		  break;
		else
		  ++idx_eps_i;
      }
      
      // now for delta_i
      tmp_flt1 = fabs( delta_i[i] - SPRT_delta_val[0] );
      idx_delta_i = 0;
      for( int j=1; j<nb_delta_values; ++j )
      {
		tmp_flt2 = fabs( epsilon_i[i] - SPRT_delta_val[j] );
		if( tmp_flt2 > tmp_flt1 )
		  break;
		else
		  ++idx_delta_i;
      }
      
      x_old = SPRT_h_i_val[idx_eps][idx_eps_i][idx_delta_i];
      if( x_old == 0.0f )
		x_old = 100.0f;
      
      x_new = x_old;
      
      // now do levenberg-marquardt for refinement
      a = delta_i[i] / epsilon_i[i];
      b = (1.0f-delta_i[i]) / (1.0f -epsilon_i[i]);
      c = 1.0f - epsilon;
      d = log(a);
      e = log(b);
      
      k = 0;
      ff = pow(a,x_old);
      fff = pow(b,x_old);
      f_1 = epsilon * d * ff + c * e * fff;
      A = f_1 * f_1;
      F_x_new = F_x_old = epsilon*ff*+c*fff-1.0f;
      g = f_1 * F_x_old;
      
      mu = 1e-6f * A;
      upsilon = 2.0f;
      
      found = fabs(g) <= 1e-6f;
      
      while( !found && k<200 )
      {
		++k;
		h_lm = -g/(A+mu);
		found = fabs(h_lm) <= 1e-6f * (fabs(x_old)+1e-6f );
		x_new = x_old + h_lm;
		F_x_new = epsilon * pow(a,x_new) + c * pow(b, x_new) - 1.0f;
		rho = 2.0f * (F_x_old * F_x_old - F_x_new * F_x_new ) / (h_lm*(mu*h_lm-g));
		if( rho > 0.0f )
		{
		  x_old = x_new;
		  f_1 = epsilon * d * pow(a,x_old) * ff + c * e * pow(b,x_old);
		  A = f_1 * f_1;
		  g = f_1 * F_x_new;
		  F_x_old = F_x_new;
		  found = fabs(g) <= 1e-6f;
		  t = 2.0f* rho - 1.0f;
		  mu *= std::max(1.0f/3.0f, 1.0f - t*t*t );
		  upsilon = 2.0f;
		}
		else
		{
		  mu *= upsilon; upsilon *= 2.0f;
		}
      }
      
      h_i[i] = x_old;
    }
    h_i_timer.Stop(); 
    
  }
  
  // now we compute the probability eta(l-1)
  float eta_l_minus_1 = 0.0f;
  float Pg = 1.0f;
  switch(nb_samples_to_take)
  {
    case 2:
    {
      Pg = epsilon * epsilon;
      break;
    }
    
    case 3:
    {
      Pg = epsilon * epsilon * epsilon;
      break;
    }
    
    case 4:
    {
      Pg = epsilon * epsilon * epsilon * epsilon;
      break;
    }
    
    case 5:
    {
      Pg = epsilon * epsilon * epsilon * epsilon * epsilon;
      break;
    }
    
    case 6:
    {
      Pg = epsilon * epsilon * epsilon * epsilon * epsilon * epsilon;
      break;
    }
    
    case 7:
    {
      Pg = epsilon * epsilon * epsilon * epsilon * epsilon * epsilon * epsilon;
      break;
    }
    
    case 8:
    {
      Pg = epsilon * epsilon * epsilon * epsilon * epsilon * epsilon * epsilon * epsilon;
      break;
    }
  }
  for( int i=0; i<l; ++i )
    eta_l_minus_1 += log( 1.0f - Pg * (1.0f - pow( A_i[i], -h_i[i] )) ) * k_i[i];
  
  float numerator = LOG_5_PER - eta_l_minus_1;
  
  if( numerator >= 0.0f )
    return 0;
  
  return uint32_t( ceil( numerator / log(1.0f - Pg * (1.0f - 1.0f/ A_i[l] ) ) ) );
}


//-----------------------------------

float RANSAC::sprt_compute_A( float eps, float delta )
{
  float Pg = 1.0f;
  switch(computation_type)
  {
    case P6pt:
    {
      Pg = eps * eps * eps * eps * eps * eps;
      break;
    }
    
    default:
    {
      break;
    }
  }
  float a = 1.0f - delta;
  float b = 1.0f - eps;
  float C = a * log(a/b) + delta * log(delta/eps);
  a = t_M * C / SPRT_m_s + 1.0f;
  float A_0 = a;
  float A_1 = a + log(A_0);
  while( fabs( A_1 - A_0 ) > 1e-6 )
  {
    A_0 = A_1;
    A_1 = a + log(A_0);
  }
  return A_1;
}
