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

#ifndef TIMER_HH
#define TIMER_HH 

/**
 *    Class to get timing results (basic stopwatch).
 *    Note that this implementation only works for Linux and Mac OS since
 *    Windows has not gettimeofday function.
 *    Reports time in seconds, accurate up to microseconds
 *
 *    Based on the timer implementation of Darko Pavic.
 *  
 *  author : Torsten Sattler (tsattler@cs.rwth-aachen.de)
 *  date : 10-02-2011
**/ 

#include <unistd.h>
#include <sys/time.h>
#include <sys/times.h>
#include <sys/types.h>

#include <string>
#include <sstream>


class Timer
{
  public:
	// constructor
	Timer();
	
	// destructor
	~Timer();
	
	// initialize timer (set elapsed time to 0)
	void Init();
	
	// Start the timer (also set elapsed time to 0)
	void Start();
	
	// Restart the timer (does not reset elapsed time)
	void Restart();
	
	// Stops the timer
	void Stop();
	
	// get the elapsed time in seconds
	double GetElapsedTime();
	
	// get the elapsed time in seconds in a string
	std::string GetElapsedTimeAsString();
	
  private:
	struct timeval time_start;
	struct timeval time_end;
	
	double elapsed_time;
  
};

#endif

