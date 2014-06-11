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


#ifndef PSEUDORANDOMNRGEN_HH
#define PSEUDORANDOMNRGEN_HH


/*
 * Simple uniform Pseudorandomnumber generators, these classes basically encapsulated the SIMD-oriented Fast Mersenne Twister (SFMT) algorithms provided from
 * http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/index.html
 *
 * Author: Torsten Sattler <tsattler@informatik.rwth-aachen.de>
 *
 * Version: 1.2
 * Date:    2011-10-02
 */



#include <stdint.h>
#include "SFMT_src/SFMT.hh"
#include <vector>
#include <list>
#include <iostream>
#include <ctime>


namespace Util {

namespace Math {


//! base class

class PRNG_base
{
  public:
	  
    //! standard constructor
    PRNG_base()
    {
      init_gen_rand((uint64_t) time(0));
    }
    
    //! constructor specifying a seed for the pseudorandom number generator
    PRNG_base(uint64_t seed)
    {
      init_gen_rand(seed);
    }

    //! generate a double pseudorandom number from the interval [0, 1]
    double generate_pseudorandom_double_number()
    {
      return genrand_real1();
    }
    
    //! generate a double pseudorandom number from the interval (0, 1)
    double generate_pseudorandom_double_number_open()
    {
      return genrand_real3();
    }
    
    //! generate a double pseudorandom number from the interval [0, 1)
    double generate_pseudorandom_double_number_rightopen()
    {
      return genrand_real2();
    }
    
    /**
     * generates number many (possibly identical) pseudorandom numbers from the interval [0, 1], 
     * results are stored in the vector. It is assumed that the vector contains at least
     * number many objects.
     **/
    void generate_double_pseudorandom_numbers( uint64_t number, std::vector< double > &numbers )
    {
      for(uint64_t i=0; i<number; i++)
	numbers[i] = genrand_real1();
    }
    
    /**
     * generates number many unique pseudorandom numbers from the interval [0, 1], 
     * results are stored in the vector. It is assumed that the vector contains at least
     * number many objects.
     * WARNING: might become slow for larger values of number
     **/
    void generate_double_pseudorandom_numbers_unique(uint64_t number, std::vector< double > &numbers )
    {
      bool found = true;
      for(uint64_t i=0; i<number; i++)
      {
		double tmp;
		found = true;
		//check for uniqueness. if the number tmp is allready present in numbers we simply generate a new pseudorandom number
		while(found)
		{
		  tmp = genrand_real1();
		  found = false;
		  for(uint64_t j=0; j<i; j++)
		  {
			if(numbers[j] == tmp)
			{
			  found = true;
			  break;
			}
		  }
		}
		numbers[i] = tmp;
      }
    }
    
    /**
     * generates number many (possibly identical) pseudorandom numbers from the interval (0, 1), 
     * results are stored in the vector. It is assumed that the vector contains at least
     * number many objects.
     **/
    void generate_double_pseudorandom_numbers_open(uint64_t number, std::vector< double > &numbers)
    {
      for(uint64_t i=0; i<number; i++)
		numbers[i] = genrand_real3();
    }
    
    /**
     * generates number many unique pseudorandom numbers from the interval (0, 1), 
     * results are stored in the vector. It is assumed that the vector contains at least
     * number many objects.
     * WARNING: might become slow for larger values of number
     **/
    void generate_double_pseudorandom_numbers_open_unique(uint64_t number, std::vector< double > &numbers)
    {
      bool found = true;
      for(uint64_t i=0; i<number; i++)
      {
		double tmp;
		found = true;
		//check for uniqueness. if the number tmp is allready present in numbers we simply generate a new pseudorandom number
		while(found)
		{
		  tmp = genrand_real3();
		  found = false;
		  for(uint64_t j=0; j<i; j++)
		  {
			if(numbers[j] == tmp)
			{
			  found = true;
			  break;
			}
		  }
		}
		numbers[i] = tmp;
      }
    }
    
    /**
     * generates number many (possible) pseudorandom numbers from the interval [0, 1), 
     * results are stored in the vector. It is assumed that the vector contains at least
     * number many objects.
     **/
    void generate_double_pseudorandom_numbers_rightopen(uint64_t number, std::vector< double > &numbers)
    {
      for(uint64_t i=0; i<number; i++)
		numbers[i] = genrand_real2();
    }
    
    /**
     * generates number many unique pseudorandom numbers from the interval [0, 1), 
     * results are stored in the vector. It is assumed that the vector contains at least
     * number many objects.
     * WARNING: might become slow for larger values of number
     **/
    void generate_double_pseudorandom_numbers_rightopen_unique(uint64_t number, std::vector< double > &numbers)
    {
      bool found = true;
      for(uint64_t i=0; i<number; i++)
      {
		double tmp;
		found = true;
		//check for uniqueness. if the number tmp is allready present in numbers we simply generate a new pseudorandom number
		while(found)
		{
		  tmp = genrand_real2();
		  found = false;
		  for(uint64_t j=0; j<i; j++)
		  {
			if(numbers[j] == tmp)
			{
			  found = true;
			  break;
			}
		  }
		}
		numbers[i] = tmp;
      }
    }
};

//! class generating pseudorandom unsigned 64bit integer numbers

class PRNG_u64i : public PRNG_base
{
  public:
    //! constructors
    PRNG_u64i():PRNG_base(){}
    PRNG_u64i(uint64_t seed):PRNG_base(seed){}
    
    //! generate an unsigned 32bit integer pseudorandom number
    uint64_t generate_pseudorandom_number();
    
    //! generate an unsigned 32bit integer pseudorandom number from the interval [start, end]
    uint64_t generate_pseudorandom_number(uint64_t start,uint64_t end);
    
    //! generate an unsigned 32bit integer pseudorandom number from the interval (start, end)
    uint64_t generate_pseudorandom_number_open(uint64_t start,uint64_t end);
    
    //! generate an unsigned 32bit integer pseudorandom number from the interval (start, end]
    uint64_t generate_pseudorandom_number_leftopen(uint64_t start,uint64_t end);
    
    //! generate an unsigned 32bit integer pseudorandom number from the interval [start, end)
    uint64_t generate_pseudorandom_number_rightopen(uint64_t start,uint64_t end);
    
    /**
     * generates number many (possibly identical) pseudorandom numbers from the interval [start, end]
     * assumes that the vector contains at least number many entries
     **/
    void generate_pseudorandom_numbers(uint64_t start, uint64_t end, uint64_t number, std::vector<uint64_t> &numbers);
    
    /**
     * generates number many unique pseudorandom numbers from the interval [start, end]
     * if number > end-start, all elements in the interval are returned
     * assumes that the vector contains at least number many entries
     **/
    void generate_pseudorandom_numbers_unique(uint64_t start, uint64_t end, uint64_t number, std::vector<uint64_t> &numbers);
    
    /**
     * generates number many (possibly identical) pseudorandom numbers from the interval (start, end)
     * assumes that the vector contains at least number many entries
     **/
    void generate_pseudorandom_numbers_open(uint64_t start, uint64_t end, uint64_t number, std::vector<uint64_t> &numbers);
    
    /**
     * generates number many unique pseudorandom numbers from the interval (start, end)
     * if number > end-start, all elements in the interval are returned
     * assumes that the vector contains at least number many entries
     **/
    void generate_pseudorandom_numbers_open_unique(uint64_t start, uint64_t end, uint64_t number, std::vector<uint64_t> &numbers);
    
    /**
     * generates number many (possibly identical) pseudorandom numbers from the interval (start, end]
     * assumes that the vector contains at least number many entries
     **/
    void generate_pseudorandom_numbers_leftopen(uint64_t start, uint64_t end, uint64_t number, std::vector<uint64_t> &numbers);
    
    /**
     * generates number many unique pseudorandom numbers from the interval (start, end]
     * if number > end-start, all elements in the interval are returned
     * assumes that the vector contains at least number many entries
     **/
    void generate_pseudorandom_numbers_leftopen_unique(uint64_t start, uint64_t end, uint64_t number, std::vector<uint64_t> &numbers);
    
    /**
     * generates number many (possibly identical) pseudorandom numbers from the interval [start, end)
     * assumes that the vector contains at least number many entries
     **/
    void generate_pseudorandom_numbers_rightopen(uint64_t start, uint64_t end, uint64_t number, std::vector<uint64_t> &numbers);
    
    /**
     * generates number many unique pseudorandom numbers from the interval [start, end)
     * if number > end-start, all elements in the interval are returned
     * assumes that the vector contains at least number many entries 
     **/
    void generate_pseudorandom_numbers_rightopen_unique(uint64_t start, uint64_t end, uint64_t number, std::vector<uint64_t> &numbers);
};


//! class generating pseudorandom 32bit uint32_teger numbers

class PRNG_u32i : public PRNG_base
{
  public:
    //! constructors
    PRNG_u32i():PRNG_base(){}
    PRNG_u32i(uint64_t seed):PRNG_base(seed){}
    
    //! generate an unsigned 32bit integer pseudorandom number
    uint32_t generate_pseudorandom_number();
    
    //! generate an unsigned 32bit integer pseudorandom number from the interval [start, end]
    uint32_t generate_pseudorandom_number(uint32_t start,uint32_t end);
    
    //! generate an unsigned 32bit integer pseudorandom number from the interval (start, end)
    uint32_t generate_pseudorandom_number_open(uint32_t start,uint32_t end);
    
    //! generate an unsigned 32bit integer pseudorandom number from the interval (start, end]
    uint32_t generate_pseudorandom_number_leftopen(uint32_t start,uint32_t end);
    
    //! generate an unsigned 32bit integer pseudorandom number from the interval [start, end)
    uint32_t generate_pseudorandom_number_rightopen(uint32_t start,uint32_t end);
    
    /**
     * generates number many (possibly identical) pseudorandom numbers from the interval [start, end]
     * assumes that the vector contains at least number many entries
     **/
    void generate_pseudorandom_numbers(uint32_t start, uint32_t end, uint32_t number, std::vector<uint32_t> &numbers);
    
    /**
     * generates number many unique pseudorandom numbers from the interval [start, end]
     * if number > end-start, all elements in the interval are returned
     * assumes that the vector contains at least number many entries
     **/
    void generate_pseudorandom_numbers_unique(uint32_t start, uint32_t end, uint32_t number, std::vector<uint32_t> &numbers);
    
    /**
     * generates number many (possibly identical) pseudorandom numbers from the interval (start, end)
     * assumes that the vector contains at least number many entries
     **/
    void generate_pseudorandom_numbers_open(uint32_t start, uint32_t end, uint32_t number, std::vector<uint32_t> &numbers);
    
    /**
     * generates number many unique pseudorandom numbers from the interval (start, end)
     * if number > end-start, all elements in the interval are returned
     * assumes that the vector contains at least number many entries
     **/
    void generate_pseudorandom_numbers_open_unique(uint32_t start, uint32_t end, uint32_t number, std::vector<uint32_t> &numbers);
    
    /**
     * generates number many (possibly identical) pseudorandom numbers from the interval (start, end]
     * assumes that the vector contains at least number many entries
     **/
    void generate_pseudorandom_numbers_leftopen(uint32_t start, uint32_t end, uint32_t number, std::vector<uint32_t> &numbers);
    
    /**
     * generates number many unique pseudorandom numbers from the interval (start, end]
     * if number > end-start, all elements in the interval are returned
     * assumes that the vector contains at least number many entries
     **/
    void generate_pseudorandom_numbers_leftopen_unique(uint32_t start, uint32_t end, uint32_t number, std::vector<uint32_t> &numbers);
    
    /**
     * generates number many (possibly identical) pseudorandom numbers from the interval [start, end)
     * assumes that the vector contains at least number many entries
     **/
    void generate_pseudorandom_numbers_rightopen(uint32_t start, uint32_t end, uint32_t number, std::vector<uint32_t> &numbers);
    
    /**
     * generates number many unique pseudorandom numbers from the interval [start, end)
     * if number > end-start, all elements in the interval are returned
     * assumes that the vector contains at least number many entries 
     **/
    void generate_pseudorandom_numbers_rightopen_unique(uint32_t start, uint32_t end, uint32_t number, std::vector<uint32_t> &numbers);
};

}

}


#endif
