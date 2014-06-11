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


/*
 * Simple Matrix definition based on GMM++,
 * supplemented by several computation routines 
 * based on CLAPACK
 *
 * Author: Martin Habbecke <habbecke@cs.rwth-aachen.de>
 *         Alexander Hornung
 *         Torsten Sattler <tsattler@cs.rwth-aachen.de>
 *
 * Version: Revision: 75
 * Date:    2011-10-02 18:00:00
 */



#ifndef _TOOLS_MATH_HH_
#define _TOOLS_MATH_HH_



#ifndef GMM_USES_LAPACK
#define GMM_USES_LAPACK
#endif

#include <algorithm>
#include <gmm/gmm.h>

#include <climits>
#include <cfloat>


//   The Math module contains several type definitions
//   for vectors and matrices which are compatible with GMM++.

//! Collection of general CV tools
namespace Util {

namespace Math {



//! Module-wide pointer to working array
extern double * p_work_double;
extern float * p_work_float;
extern long int * p_work_longint;

//! Module-wide size of working array
extern int m_workSize_double;
extern int m_workSize_float;
extern int m_workSize_longint;




//! Simple vector type of variable length, based on <tt>double</tt>
typedef std::vector< double > Vector;

//! Simple vector type of variable length, based on <tt>double</tt>
typedef std::vector< long int > Vectorli;
typedef std::vector< int > Vectori;


//! Simple vector type of variable length, based on <tt>float</tt>
typedef std::vector< float > VectorSinglePrecision;

//! Matrix type based on <tt>double</tt>, stored column-major
typedef gmm::dense_matrix< double > Matrix;

//! Matrix type based on <tt>float</tt>, stored column-major
typedef gmm::dense_matrix< float > MatrixSinglePrecision;



/////////////////////////////////////////////////////////
// module functions

//! Function to free memory of working array
void FreeWorkingArray( void );




///////////////////////////////////////////////////////////////
// interface to lapack


//! Function that computes SVD A = U * D * V<sup>T</sup> and destroys Matrix A
bool SVD( Matrix & _mat_A, Matrix & _mat_U, 
	  Matrix & _mat_VT, Vector & _vec_D );

//! Function that only computes matrix V<sup>T</sup> of SVD A = U * D * V<sup>T</sup> and destroys Matrix A
bool SVD_V( Matrix & _mat_A, Matrix & _mat_VT );




};
};




#endif 
