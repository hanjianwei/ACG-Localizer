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

#include "timer.hh"

#include <cmath>

Timer::Timer()
{
  elapsed_time = 0.0;
}

//-----------------------------------

Timer::~Timer()
{}

//-----------------------------------

void Timer::Init()
{
  elapsed_time = 0.0;
}

//-----------------------------------

void Timer::Start()
{
  gettimeofday( &time_start, 0 );
  elapsed_time = 0.0;
}

//-----------------------------------

void Timer::Restart()
{
  gettimeofday( &time_start, 0 );
}

//-----------------------------------

void Timer::Stop()
{
  gettimeofday( & time_end, 0 );
  
  // compute the elapsed time in seconds, i.e. we have to convert from microseconds to seconds
  elapsed_time += ( (double) time_end.tv_sec - (double) time_start.tv_sec + ( (double) time_end.tv_usec - (double) time_start.tv_usec)/1e6 );
}

//-----------------------------------

double Timer::GetElapsedTime()
{
  return elapsed_time;
}

//-----------------------------------

std::string Timer::GetElapsedTimeAsString()
{
  // parse elapsed time into string stream
  std::ostringstream s;
  // get the elapsed minutes and seconds
  double elapsed_minutes = (elapsed_time<60)?0.0:floor(elapsed_time/60);
  double elapsed_seconds = elapsed_time - floor(elapsed_time/60) * 60;
  s << elapsed_minutes << " minutes " << elapsed_seconds << " seconds ";
  return s.str();
}


