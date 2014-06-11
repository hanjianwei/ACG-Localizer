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

#include "pseudorandomnrgen.hh"

namespace Util {

namespace Math {

uint64_t PRNG_u64i::generate_pseudorandom_number()
{
  return gen_rand32();
}


uint64_t PRNG_u64i::generate_pseudorandom_number(uint64_t start,uint64_t end)
{
  return ((uint64_t) (generate_pseudorandom_double_number() * (end-start))) + start;
}


uint64_t PRNG_u64i::generate_pseudorandom_number_open(uint64_t start,uint64_t end) 
{
  return ((uint64_t)(generate_pseudorandom_double_number_open() * (end-start))) + start;	
}


uint64_t PRNG_u64i::generate_pseudorandom_number_leftopen(uint64_t start,uint64_t end)
{
  return generate_pseudorandom_number(start+1,end);	
}


uint64_t PRNG_u64i::generate_pseudorandom_number_rightopen(uint64_t start,uint64_t end)
{
  return ((uint64_t) (generate_pseudorandom_double_number_rightopen() * (end-start))) + start;	
}



void PRNG_u64i::generate_pseudorandom_numbers(uint64_t start, uint64_t end, uint64_t number, std::vector<uint64_t> &numbers)
{
  for(uint64_t i=0; i<number; i++)
    numbers[i] = generate_pseudorandom_number(start,end);
}


void PRNG_u64i::generate_pseudorandom_numbers_unique(uint64_t start, uint64_t end, uint64_t number, std::vector<uint64_t> &numbers)
{
  if(number >= (end-start)+1)
  {
	  for(uint64_t i=start; i<=end; i++)
	    numbers[i-start] = i;
  }
  else
  {
    uint64_t interval_length = end - start + 1;
    std::vector<uint64_t> sorted_numbers(number);
    sorted_numbers.clear();
    uint64_t tmp = generate_pseudorandom_number(start, end);
    std::vector<uint64_t>::iterator it;
    for(uint64_t i=0; i<number; i++)
    {
	  tmp = (uint64_t)(generate_pseudorandom_double_number()*(interval_length - number));
      for( it = sorted_numbers.begin(); it != sorted_numbers.end(); it++)
      {
		if(*it > tmp)
		  break;
		else
		  tmp++;
      }
      if(it == sorted_numbers.end())
      {
		sorted_numbers.push_back(start+tmp);
      }
      else
      {
		sorted_numbers.insert(it,start + tmp);
      }
      numbers[i] = start + tmp;
    }
  }
}


void PRNG_u64i::generate_pseudorandom_numbers_open(uint64_t start, uint64_t end, uint64_t number, std::vector<uint64_t> &numbers)
{
  generate_pseudorandom_numbers(start+1,end-1,number,numbers);
}


void PRNG_u64i::generate_pseudorandom_numbers_open_unique(uint64_t start, uint64_t end, uint64_t number, std::vector<uint64_t> &numbers)
{
  generate_pseudorandom_numbers_unique(start+1,end-1,number,numbers);
}


void PRNG_u64i::generate_pseudorandom_numbers_leftopen(uint64_t start, uint64_t end, uint64_t number, std::vector<uint64_t> &numbers)
{
  generate_pseudorandom_numbers(start+1,end,number,numbers);
}


void PRNG_u64i::generate_pseudorandom_numbers_leftopen_unique(uint64_t start, uint64_t end, uint64_t number, std::vector<uint64_t> &numbers)
{
  generate_pseudorandom_numbers_unique(start+1,end,number,numbers);
}


void PRNG_u64i::generate_pseudorandom_numbers_rightopen(uint64_t start, uint64_t end, uint64_t number, std::vector<uint64_t> &numbers)
{
  generate_pseudorandom_numbers(start,end-1,number,numbers);
}


void PRNG_u64i::generate_pseudorandom_numbers_rightopen_unique(uint64_t start, uint64_t end, uint64_t number, std::vector<uint64_t> &numbers)
{
  generate_pseudorandom_numbers_unique(start,end-1,number,numbers);
}

//------------------------------------


uint32_t PRNG_u32i::generate_pseudorandom_number()
{
  return (uint32_t) gen_rand32();
}


uint32_t PRNG_u32i::generate_pseudorandom_number(uint32_t start,uint32_t end)
{
  return ((uint32_t) (generate_pseudorandom_double_number() * (end-start))) + start;
}


uint32_t PRNG_u32i::generate_pseudorandom_number_open(uint32_t start,uint32_t end)
{
  return ((uint64_t)(generate_pseudorandom_double_number_open() * (end-start))) + start;	
}


uint32_t PRNG_u32i::generate_pseudorandom_number_leftopen(uint32_t start,uint32_t end)
{
  return generate_pseudorandom_number(start+1,end);	
}


uint32_t PRNG_u32i::generate_pseudorandom_number_rightopen(uint32_t start,uint32_t end)
{
  return ((uint32_t) (generate_pseudorandom_double_number_rightopen() * (end-start))) + start;	
}



void PRNG_u32i::generate_pseudorandom_numbers(uint32_t start, uint32_t end, uint32_t number, std::vector<uint32_t> &numbers)
{
  for(uint32_t i=0; i<number; i++)
    numbers[i] = generate_pseudorandom_number(start,end);
}


void PRNG_u32i::generate_pseudorandom_numbers_unique(uint32_t start, uint32_t end, uint32_t number, std::vector<uint32_t> &numbers)
{
  if(number >= (end-start)+1)
  {
	  for(uint32_t i=start; i<=end; i++)
	    numbers[i-start] = i;
  }
  else
  {
    uint32_t interval_length = end - start + 1;
    std::vector<uint32_t> sorted_numbers(number);
    sorted_numbers.clear();
    uint32_t tmp = generate_pseudorandom_number(start, end);
    std::vector<uint32_t>::iterator it;
    for(uint32_t i=0; i<number; i++)
    {
	  tmp = (uint32_t)(generate_pseudorandom_double_number()*(interval_length - number));
      for( it = sorted_numbers.begin(); it != sorted_numbers.end(); it++)
      {
		if(*it > tmp)
		  break;
		else
		  tmp++;
      }
      if(it == sorted_numbers.end())
      {
		sorted_numbers.push_back(start+tmp);
      }
      else
      {
		sorted_numbers.insert(it,start + tmp);
      }
      numbers[i] = start + tmp;
    }
  }
}


void PRNG_u32i::generate_pseudorandom_numbers_open(uint32_t start, uint32_t end, uint32_t number, std::vector<uint32_t> &numbers)
{
  generate_pseudorandom_numbers(start+1,end-1,number,numbers);
}


void PRNG_u32i::generate_pseudorandom_numbers_open_unique(uint32_t start, uint32_t end, uint32_t number, std::vector<uint32_t> &numbers)
{
  generate_pseudorandom_numbers_unique(start+1,end-1,number,numbers);
}


void PRNG_u32i::generate_pseudorandom_numbers_leftopen(uint32_t start, uint32_t end, uint32_t number, std::vector<uint32_t> &numbers)
{
  generate_pseudorandom_numbers(start+1,end,number,numbers);
}


void PRNG_u32i::generate_pseudorandom_numbers_leftopen_unique(uint32_t start, uint32_t end, uint32_t number, std::vector<uint32_t> &numbers)
{
  generate_pseudorandom_numbers_unique(start+1,end,number,numbers);
}


void PRNG_u32i::generate_pseudorandom_numbers_rightopen(uint32_t start, uint32_t end, uint32_t number, std::vector<uint32_t> &numbers)
{
  generate_pseudorandom_numbers(start,end-1,number,numbers);
}


void PRNG_u32i::generate_pseudorandom_numbers_rightopen_unique(uint32_t start, uint32_t end, uint32_t number, std::vector<uint32_t> &numbers)
{
  generate_pseudorandom_numbers_unique(start,end-1,number,numbers);
}

}

}
