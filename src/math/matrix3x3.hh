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


#ifndef _MATRIX_3X3_HH_
#define _MATRIX_3X3_HH_

#include <iostream>
#include <vector>
#include <map>
#include <cmath>

#include "math.hh"
#include "matrixbase.hh"



extern "C" {
    int dgetri_( long int *n, double *a, long int *lda, long int *ipiv,
		 double *work, long int *lwork, long int *info );

    int dgetrf_( long int *m, long int *n, double *a, long int * lda,
		 long int *ipiv, long int *info );
}

namespace Util {

namespace Math {


//! Abstraction of an image homography
template< typename Scalar >
class Matrix3x3T : public MatrixSqrT< Scalar, 3 > {

public:

    //! Base type definition
    typedef MatrixSqrT< Scalar, 3 > Base;

    //! Own type definition
    typedef Matrix3x3T< Scalar > This;

    //! 3-dimensional OpenMesh vector
    typedef OpenMesh::VectorT< Scalar, 3 > Vec3;

    //! Default constructor
    Matrix3x3T()
	: Base()
    {
    }

    //! Copy constructor
    Matrix3x3T( const MatrixT< Scalar, 3, 3 > & _rhs )
	: Base( _rhs )
    {
    }

    //! Assignment operator
    inline const Matrix3x3T & operator=( const MatrixT< Scalar, 3, 3 > & _rhs )
    {
	  if( &_rhs == this )
		  return *this;

	  memcpy( Base::mp_data, _rhs.mp_data, 9 * sizeof( Scalar ) );
	  return *this;
    }




    //! Simple multiplication with a 3-vector from right
    template< typename ScalarA, typename ScalarB >
    inline void multRight( const OpenMesh::VectorT< ScalarA, 3 > & _a,
			   OpenMesh::VectorT< ScalarB, 3 > & _b ) const
    {
	  _b[ 0 ] = (ScalarB) (Base::mp_data[0] * _a[0] 
				  + Base::mp_data[3] * _a[1] 
				  + Base::mp_data[6] * _a[2]);
	  _b[ 1 ] = (ScalarB) (Base::mp_data[1] * _a[0] 
				  + Base::mp_data[4] * _a[1] 
				  + Base::mp_data[7] * _a[2]);
	  _b[ 2 ] = (ScalarB) (Base::mp_data[2] * _a[0] 
				  + Base::mp_data[5] * _a[1] 
				  + Base::mp_data[8] * _a[2]);
    }

    //! Simple multiplication with a 2-vector from right
    template< typename ScalarA, typename ScalarB >
    inline void multRight( const OpenMesh::VectorT< ScalarA, 2 > & _a,
			   OpenMesh::VectorT< ScalarB, 2 > & _b ) const
    {
	  _b[ 0 ] = (ScalarB) (Base::mp_data[0] * _a[0] 
				  + Base::mp_data[3] * _a[1] 
				  + Base::mp_data[6]);
	  _b[ 1 ] = (ScalarB) (Base::mp_data[1] * _a[0] 
				  + Base::mp_data[4] * _a[1] 
				  + Base::mp_data[7]);
    }

    //! Simple multiplication with a 2-vector from right
    template< typename ScalarA, typename ScalarB >
    inline void multRight( const ScalarA _ax, const ScalarA _ay,
			   ScalarB & _bx, ScalarB & _by ) const
    {
	  _bx = (ScalarB) (Base::mp_data[0] * _ax 
			  + Base::mp_data[3] * _ay 
			  + Base::mp_data[6]);
	  _by = (ScalarB) (Base::mp_data[1] * _ax 
			  + Base::mp_data[4] * _ay 
			  + Base::mp_data[7]);
    }




    //! Simple multiplication with a 3-vector from left
    template< typename ScalarA, typename ScalarB >
    inline void multLeft( const OpenMesh::VectorT< ScalarA, 3 > & _a, 
			  OpenMesh::VectorT< ScalarB, 3 > & _b ) const
    {
	  _b[ 0 ] = (ScalarB) (Base::mp_data[0] * _a[0]
				  + Base::mp_data[1] * _a[1] 
				  + Base::mp_data[2] * _a[2]);
	  _b[ 1 ] = (ScalarB) (Base::mp_data[3] * _a[0]
				  + Base::mp_data[4] * _a[1]
				  + Base::mp_data[5] * _a[2]);
	  _b[ 2 ] = (ScalarB) (Base::mp_data[6] * _a[0]
				  + Base::mp_data[7] * _a[1]
				  + Base::mp_data[8] * _a[2]);
    }

    //! Simple multiplication with a 2-vector from right
    template< typename ScalarA, typename ScalarB >
    inline void multLeft( const OpenMesh::VectorT< ScalarA, 2 > & _a, 
			  OpenMesh::VectorT< ScalarB, 2 > & _b ) const
    {
	  _b[ 0 ] = (ScalarB) (Base::mp_data[0] * _a[0]
				  + Base::mp_data[1] * _a[1] 
				  + Base::mp_data[2]);
	  _b[ 1 ] = (ScalarB) (Base::mp_data[3] * _a[0]
				  + Base::mp_data[4] * _a[1]
				  + Base::mp_data[5]);
    }

    //! Simple multiplication with a 2-vector from right
    template< typename ScalarA, typename ScalarB >
    inline void multLeft(  const ScalarA _ax, const ScalarA _ay,
			   ScalarB & _bx, ScalarB & _by ) const
    {
	  _bx = (ScalarB) (Base::mp_data[0] * _ax
			  + Base::mp_data[1] * _ay 
			  + Base::mp_data[2]);
	  _by = (ScalarB) (Base::mp_data[3] * _ax
			  + Base::mp_data[4] * _ay
			  + Base::mp_data[5]);
    }






    //! Multiply with another homography from the right: this := this * B
    template< typename Scalar2 >
    inline void multRight( const MatrixT< Scalar2, 3, 3 > & _rhs )
    {
	  #define T cpy.mp_data
	  #define B _rhs.mp_data

		  Matrix3x3T cpy( *this );

		  Base::mp_data[ 0 ] = (Scalar) (T[0] * B[0] + T[3] * B[1] + T[6] * B[2]);
		  Base::mp_data[ 1 ] = (Scalar) (T[1] * B[0] + T[4] * B[1] + T[7] * B[2]);
		  Base::mp_data[ 2 ] = (Scalar) (T[2] * B[0] + T[5] * B[1] + T[8] * B[2]);

		  Base::mp_data[ 3 ] = (Scalar) (T[0] * B[3] + T[3] * B[4] + T[6] * B[5]);
		  Base::mp_data[ 4 ] = (Scalar) (T[1] * B[3] + T[4] * B[4] + T[7] * B[5]);
		  Base::mp_data[ 5 ] = (Scalar) (T[2] * B[3] + T[5] * B[4] + T[8] * B[5]);

		  Base::mp_data[ 6 ] = (Scalar) (T[0] * B[6] + T[3] * B[7] + T[6] * B[8]);
		  Base::mp_data[ 7 ] = (Scalar) (T[1] * B[6] + T[4] * B[7] + T[7] * B[8]);
		  Base::mp_data[ 8 ] = (Scalar) (T[2] * B[6] + T[5] * B[7] + T[8] * B[8]);
	  #undef T
	  #undef B
    }


    //! Multiply with another homography from the left: this := A * this
    template< typename Scalar2 >
    inline void multLeft( const MatrixT< Scalar2, 3, 3 > & _lhs )
    {
	  #define A _lhs.mp_data
	  #define T cpy.mp_data

		  Matrix3x3T cpy( *this );

		  Base::mp_data[ 0 ] = (Scalar) (A[0] * T[0] + A[3] * T[1] + A[6] * T[2]);
		  Base::mp_data[ 1 ] = (Scalar) (A[1] * T[0] + A[4] * T[1] + A[7] * T[2]);
		  Base::mp_data[ 2 ] = (Scalar) (A[2] * T[0] + A[5] * T[1] + A[8] * T[2]);

		  Base::mp_data[ 3 ] = (Scalar) (A[0] * T[3] + A[3] * T[4] + A[6] * T[5]);
		  Base::mp_data[ 4 ] = (Scalar) (A[1] * T[3] + A[4] * T[4] + A[7] * T[5]);
		  Base::mp_data[ 5 ] = (Scalar) (A[2] * T[3] + A[5] * T[4] + A[8] * T[5]);

		  Base::mp_data[ 6 ] = (Scalar) (A[0] * T[6] + A[3] * T[7] + A[6] * T[8]);
		  Base::mp_data[ 7 ] = (Scalar) (A[1] * T[6] + A[4] * T[7] + A[7] * T[8]);
		  Base::mp_data[ 8 ] = (Scalar) (A[2] * T[6] + A[5] * T[7] + A[8] * T[8]);
	  #undef A
	  #undef T
    }


    //! Initialize skew-symmetric matrix from given 3-vector
    inline void skewSymmetricMatrix( const Vec3 & _vec )
    {
	  this->clear();

	  (*this)( 0, 1 ) = -_vec[ 2 ];
	  (*this)( 0, 2 ) =  _vec[ 1 ];
	  (*this)( 1, 0 ) =  _vec[ 2 ];
	  (*this)( 1, 2 ) = -_vec[ 0 ];
	  (*this)( 2, 0 ) = -_vec[ 1 ];
	  (*this)( 2, 1 ) =  _vec[ 0 ];
    }


    //! Compute rotation parameters from current matrix
    inline void decomposeRotation( Vec3 & _axis, Scalar & _angle ) const
    {
	  Scalar cosPhi, sinPhi;
	  This matEV;
	  Vec3 vecEVR, vecEVI;
	  int idx, maxIdx;
	  Scalar maxEV;

	  // eigen decomposition
	  if( NonsymRightEigenproblem_copy( *this, matEV, vecEVR, vecEVI ) )
	  {	
		  maxIdx = -1;
		  maxEV = -DBL_MAX;

		  for( idx = 0; idx < 3; ++idx )
		  {
  // 		std::cout << "-- EV " << idx << ": " << vecEVR[ idx ]
  // 			  << " + " << vecEVI[ idx ] << "i" << std::endl;

		  if( vecEVR[ idx ] > maxEV )
		  {
			  maxEV = vecEVR[ idx ];
			  maxIdx = idx;
		  }
		  }

  // 	    std::cout << "-- max EV: " << maxIdx << std::endl;

		  _axis[ 0 ] = matEV( 0, maxIdx );
		  _axis[ 1 ] = matEV( 1, maxIdx );
		  _axis[ 2 ] = matEV( 2, maxIdx );
	  }
	  else
	  {
		  std::cout << "********* ROTATION decomposition failed ***********"
				<< std::endl;
		  _axis[ 0 ] = 1.0;
		  _axis[ 1 ] = _axis[ 2 ] = 0.0;
	  }

	  sinPhi = 0.5 * ( _axis[ 0 ] * (Base::mp_data[ 5 ] - Base::mp_data[ 7 ])
			  + _axis[ 1 ] * (Base::mp_data[ 6 ] - Base::mp_data[ 2 ])
			  + _axis[ 2 ] * (Base::mp_data[ 1 ] - Base::mp_data[ 3 ]) );
	  cosPhi = 0.5 * (Base::mp_data[ 0 ] 
			  + Base::mp_data[ 4 ] 
			  + Base::mp_data[ 8 ] - 1.0);

	  _angle = atan2( sinPhi, cosPhi );
    }


    //! Compute rotation parameters from current matrix
    inline void decomposeRotation( Vec3 & _axis ) const
    {
	  Scalar angle;

	  decomposeRotation( _axis, angle );
	  
	  if( angle < 0.0 )
	  {
		  _axis *= -1.0;
		  angle *= -1.0;
	  }

	  if( angle > M_PI )
	  {
		  angle -= 2.0 * M_PI;
	  }
	  
	  _axis *= angle;
    }


    //! Construct 3D rotation matrix from given parameters
    inline void composeRotation( const Vec3 & _axis, const Scalar & _angle )
    {
	  Scalar s = sin( _angle );
	  Scalar c = cos( _angle );
	  Scalar u = 1.0 - c;

	  Base::mp_data[ 0 ] = _axis[ 0 ] * _axis[ 0 ] * u + c;
	  Base::mp_data[ 1 ] = _axis[ 0 ] * _axis[ 1 ] * u + _axis[ 2 ] * s;
	  Base::mp_data[ 2 ] = _axis[ 0 ] * _axis[ 2 ] * u - _axis[ 1 ] * s;

	  Base::mp_data[ 3 ] = _axis[ 0 ] * _axis[ 1 ] * u - _axis[ 2 ] * s;
	  Base::mp_data[ 4 ] = _axis[ 1 ] * _axis[ 1 ] * u + c;
	  Base::mp_data[ 5 ] = _axis[ 1 ] * _axis[ 2 ] * u + _axis[ 0 ] * s;

	  Base::mp_data[ 6 ] = _axis[ 0 ] * _axis[ 2 ] * u + _axis[ 1 ] * s;
	  Base::mp_data[ 7 ] = _axis[ 1 ] * _axis[ 2 ] * u - _axis[ 0 ] * s;
	  Base::mp_data[ 8 ] = _axis[ 2 ] * _axis[ 2 ] * u + c;
    }


    //! Construct 3D rotation matrix from given parameters
    inline void composeRotation( const Vec3 & _axis )
    {
	  Scalar angle, scale;
	  Vec3 axis( _axis );

	  angle = _axis.length();

	  if( angle < 1e-6 )
	  {
  //	    std::cout << "** rotation: short axis" << std::endl;

		  this->skewSymmetricMatrix( _axis );
		  Base::mp_data[ 0 ] = 1.0;
		  Base::mp_data[ 4 ] = 1.0;
		  Base::mp_data[ 8 ] = 1.0;
	  }
	  else
	  {
		  axis *= 1.0 / angle;

		  composeRotation( axis, angle );
	  }
    }





    //! Decompose left 3x3 block into a triangular matrix K and a rotation matrix R
    void decompose( Matrix3x3T< Scalar > & _matK, Matrix3x3T< Scalar > & _matR ) const
//		    Scalar & _scale ) const
    {
	  Matrix3x3T< Scalar > matRot;

	  ////////////////////////////////////
	  // Algorithm taken from 
	  // Hartley, Zisserman: "Multiple View Geometry in Computer Vision",
	  // Appendix 4.1


	  double s, c;

	  ////////////////////////////////////
	  // Cancellation of element 3,2

	  memcpy( _matK.mp_data, this->mp_data, 9 * sizeof( Scalar ) );

	  s = - _matK( 2, 1 ) / sqrt( _matK( 2, 1 ) * _matK( 2, 1 )
					  + _matK( 2, 2 ) * _matK( 2, 2 ) );
	  
	  c = _matK( 2, 2 ) / sqrt( _matK( 2, 1 ) * _matK( 2, 1 )
					+ _matK( 2, 2 ) * _matK( 2, 2 ) );
	  
	  matRot.clear();
	  matRot( 0, 0 ) = 1.0;
	  matRot( 1, 1 ) = matRot( 2, 2 ) = c;
	  matRot( 2, 1 ) = s;
	  matRot( 1, 2 ) = -s;

	  _matK.multRight( matRot );

	  matRot.transpose();
	  _matR = matRot;
	  

	  ////////////////////////////////////
	  // Cancellation of element 3,1
	  
	  s = _matK( 2, 0 ) / sqrt( _matK( 2, 0 ) * _matK( 2, 0 )
					+ _matK( 2, 2 ) * _matK( 2, 2 ) );

	  c = _matK( 2, 2 ) / sqrt( _matK( 2, 0 ) * _matK( 2, 0 )
					+ _matK( 2, 2 ) * _matK( 2, 2 ) );    

	  matRot.clear();
	  matRot( 1, 1 ) = 1.0;
	  matRot( 0, 0 ) = matRot( 2, 2 ) = c;
	  matRot( 0, 2 ) = s;
	  matRot( 2, 0 ) = -s;

	  _matK.multRight( matRot );

	  matRot.transpose();
	  _matR.multLeft( matRot );


	  /////////////////////////////////////
	  // Cancellation of element 2,1

	  s = _matK( 1, 0 ) / sqrt( _matK( 1, 0 ) * _matK( 1, 0 )
					+ _matK( 1, 1 ) * _matK( 1, 1 ) );

	  c = - _matK( 1, 1 ) / sqrt( _matK( 1, 0 ) * _matK( 1, 0 )
					  + _matK( 1, 1 ) * _matK( 1, 1 ) );  

	  matRot.clear();
	  matRot( 2, 2 ) = 1.0;
	  matRot( 0, 0 ) = matRot( 1, 1 ) = c;
	  matRot( 1, 0 ) = s;
	  matRot( 0, 1 ) = -s;

	  _matK.multRight( matRot );

	  matRot.transpose();
	  _matR.multLeft( matRot );


	  // scale calibration matrix such that element 2,2 equals 1
	  double _scale = 1.0 / _matK( 2, 2 );
	  if( !isinf( _scale ) && !isnan( _scale ) )
		  _matK.scale( _scale );
    }



    //! Compute the determinant of the 3x3 matrix
    const Scalar determinant( void ) const;

    //! Compute the inverse
    bool invert( void )
    {
	  const Matrix3x3T< Scalar > matCpy( *this );

	  double det = determinant();

  // 	if( det < 1e-8 && det > -1e-8 )
  // 	    return false;

	  det = 1.0 / det;

  #define A matCpy.mp_data
	  Base::mp_data[ 0 ] = det * (A[8] * A[4] - A[5] * A[7]);
	  Base::mp_data[ 3 ] = det * (-(A[8] * A[3] - A[5] * A[6]));
	  Base::mp_data[ 6 ] = det * (A[7] * A[3] - A[4] * A[6]);

	  Base::mp_data[ 1 ] = det * (-(A[8] * A[1] - A[2] * A[7]));
	  Base::mp_data[ 4 ] = det * (A[8] * A[0] - A[2]* A[6]);
	  Base::mp_data[ 7 ] = det * (-(A[7] * A[0] - A[1] * A[6]));

	  Base::mp_data[ 2 ] = det * (A[5] * A[1] - A[2] * A[4]);
	  Base::mp_data[ 5 ] = det * (-(A[5] * A[0] - A[2] * A[3]));
	  Base::mp_data[ 8 ] = det * (A[4] * A[0] - A[1] * A[3]);
  #undef A

	  return true;
    }
    

};


template< typename Scalar >
const Scalar Matrix3x3T< Scalar >::determinant( void ) const
{
#define A Base::mp_data

    return A[ 0 ] * (A[ 4 ] * A[ 8 ] - A[ 5 ] * A[ 7 ])
 	 + A[ 1 ] * (A[ 6 ] * A[ 5 ] - A[ 8 ] * A[ 3 ])
	 + A[ 2 ] * (A[ 3 ] * A[ 7 ] - A[ 4 ] * A[ 6 ]);

#undef A
}




template< typename Scalar >
const OpenMesh::VectorT< Scalar, 3 >
operator*( const Matrix3x3T< Scalar > & _lhs,
	   const OpenMesh::VectorT< Scalar, 3 > & _vec )
{
    OpenMesh::VectorT< Scalar, 3 > result;

    _lhs.multRight( _vec, result );
    
    return result;
}


template< typename Scalar >
const Matrix3x3T< Scalar >
operator*( const Matrix3x3T< Scalar > & _lhs,
	   const Scalar & _scale )
{
    Matrix3x3T< Scalar > result;

    result.scale( _scale );

    return result;
}




//! Type of a homography
typedef Matrix3x3T< double > Matrix3x3;
typedef Matrix3x3T< float > Matrix3x3f;

//! Type for a map of homographies
typedef std::map< int, Matrix3x3 > Matrix3x3Map;

//! Type for a vector of homographies
typedef std::vector< Matrix3x3 > Matrix3x3Vec;



}

}

#endif // _MATRIX_3X3_HH_
