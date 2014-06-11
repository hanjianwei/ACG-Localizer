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

#ifndef _MATRIX_4X4_HH_
#define _MATRIX_4X4_HH_

#include <iostream>
#include <vector>
#include <map>

#include "math.hh"
#include "matrixbase.hh"



namespace Util {

namespace Math {


template< typename Scalar >
class ProjMatrixT;


//! Abstraction of a 3-space homography, i.e. a 4x4 matrix
template< typename Scalar >
class Matrix4x4T : public MatrixSqrT< Scalar, 4 > {

public:

    //! Base type definition
    typedef MatrixSqrT< Scalar, 4 > Base;

    //! 3-dimensional OpenMesh vector
    typedef OpenMesh::VectorT< Scalar, 3 > Vec3;

    //! 4-dimensional OpenMesh vector
    typedef OpenMesh::VectorT< Scalar, 4 > Vec4;

    Matrix4x4T()
	: Base()
    {
    }

    //! Copy constructor
    Matrix4x4T( const MatrixT< Scalar, 4, 4 > & _rhs )
	: Base( _rhs )
    {
    }

    //! Assignment operator
    inline const Matrix4x4T & operator=( const MatrixT< Scalar, 4, 4 > & _rhs )
    {
	  memcpy( Base::mp_data, _rhs.mp_data, 16 * sizeof( Scalar ) );
	  return *this;
    }


    //! Simple multiplication with a 4-vector from right
    inline void multRight( const Vec4 & _a, Vec4 & _b ) const
    {
	  _b[ 0 ] = Base::mp_data[0] * _a[0] + Base::mp_data[4] * _a[1]
		  + Base::mp_data[8 ] * _a[2] + Base::mp_data[12] * _a[3];

	  _b[ 1 ] = Base::mp_data[1] * _a[0] + Base::mp_data[5] * _a[1] 
		  + Base::mp_data[9 ] * _a[2] + Base::mp_data[13] * _a[3];

	  _b[ 2 ] = Base::mp_data[2] * _a[0] + Base::mp_data[6] * _a[1]
		  + Base::mp_data[10] * _a[2] + Base::mp_data[14] * _a[3];

	  _b[ 3 ] = Base::mp_data[3] * _a[0] + Base::mp_data[7] * _a[1]
		  + Base::mp_data[11] * _a[2] + Base::mp_data[15] * _a[3];
    }


    //! Simple multiplication with a (split) 4-vector from right
    inline void multRight( const Vec3 & _a, const Scalar & _a3,
			   Vec3 & _b, Scalar & _b3  ) const
    {
	  _b[ 0 ] = Base::mp_data[0] * _a[0] + Base::mp_data[4] * _a[1]
		  + Base::mp_data[8 ] * _a[2] + Base::mp_data[12] * _a3;

	  _b[ 1 ] = Base::mp_data[1] * _a[0] + Base::mp_data[5] * _a[1] 
		  + Base::mp_data[9 ] * _a[2] + Base::mp_data[13] * _a3;

	  _b[ 2 ] = Base::mp_data[2] * _a[0] + Base::mp_data[6] * _a[1]
		  + Base::mp_data[10] * _a[2] + Base::mp_data[14] * _a3;

	  _b3     = Base::mp_data[3] * _a[0] + Base::mp_data[7] * _a[1]
		  + Base::mp_data[11] * _a[2] + Base::mp_data[15] * _a3;
    }


    //! Simple multiplication with a 4-vector from right
    inline void multLeft( const Vec4 & _a, Vec4 & _b ) const
    {
	  _b[ 0 ] = Base::mp_data[0] * _a[0] + Base::mp_data[1] * _a[1]
		  + Base::mp_data[2] * _a[2] + Base::mp_data[3] * _a[3];

	  _b[ 1 ] = Base::mp_data[4] * _a[0] + Base::mp_data[5] * _a[1]
		  + Base::mp_data[6] * _a[2] + Base::mp_data[7] * _a[3];

	  _b[ 2 ] = Base::mp_data[8] * _a[0] + Base::mp_data[9] * _a[1]
		  + Base::mp_data[10] * _a[2] + Base::mp_data[11] * _a[3];

	  _b[ 3 ] = Base::mp_data[12] * _a[0] + Base::mp_data[13] * _a[1]
		  + Base::mp_data[14] * _a[2] + Base::mp_data[15] * _a[3];
    }


    //! Simple multiplication with a (split) 4-vector from right
    inline void multLeft( const Vec3 & _a, const Scalar & _a3,
			  Vec3 & _b, Scalar & _b3 ) const
    {
	  _b[ 0 ] = Base::mp_data[0] * _a[0] + Base::mp_data[1] * _a[1]
		  + Base::mp_data[2] * _a[2] + Base::mp_data[3] * _a3;

	  _b[ 1 ] = Base::mp_data[4] * _a[0] + Base::mp_data[5] * _a[1]
		  + Base::mp_data[6] * _a[2] + Base::mp_data[7] * _a3;

	  _b[ 2 ] = Base::mp_data[8] * _a[0] + Base::mp_data[9] * _a[1]
		  + Base::mp_data[10] * _a[2] + Base::mp_data[11] * _a3;

	  _b3     = Base::mp_data[12] * _a[0] + Base::mp_data[13] * _a[1]
		  + Base::mp_data[14] * _a[2] + Base::mp_data[15] * _a3;
    }


    //! Multiply with another matrix from the right: this := this * B
    inline void multRight( const MatrixT< Scalar, 4, 4 > & _rhs )
    {
	  #define T cpy.Base::mp_data
	  #define B _rhs.mp_data

		  Matrix4x4T cpy( *this );

		  Base::mp_data[ 0 ] = T[0] * B[0] + T[4] * B[1] + T[8]  * B[2] + T[12] * B[3];
		  Base::mp_data[ 1 ] = T[1] * B[0] + T[5] * B[1] + T[9]  * B[2] + T[13] * B[3];
		  Base::mp_data[ 2 ] = T[2] * B[0] + T[6] * B[1] + T[10] * B[2] + T[14] * B[3];
		  Base::mp_data[ 3 ] = T[3] * B[0] + T[7] * B[1] + T[11] * B[2] + T[15] * B[3];

		  Base::mp_data[ 4 ] = T[0] * B[4] + T[4] * B[5] + T[8]  * B[6] + T[12] * B[7];
		  Base::mp_data[ 5 ] = T[1] * B[4] + T[5] * B[5] + T[9]  * B[6] + T[13] * B[7];
		  Base::mp_data[ 6 ] = T[2] * B[4] + T[6] * B[5] + T[10] * B[6] + T[14] * B[7];
		  Base::mp_data[ 7 ] = T[3] * B[4] + T[7] * B[5] + T[11] * B[6] + T[15] * B[7];

		  Base::mp_data[ 8 ] = T[0] * B[8] + T[4] * B[9] + T[8]  * B[10] + T[12] * B[11];
		  Base::mp_data[ 9 ] = T[1] * B[8] + T[5] * B[9] + T[9]  * B[10] + T[13] * B[11];
		  Base::mp_data[10 ] = T[2] * B[8] + T[6] * B[9] + T[10] * B[10] + T[14] * B[11];
		  Base::mp_data[11 ] = T[3] * B[8] + T[7] * B[9] + T[11] * B[10] + T[15] * B[11];

		  Base::mp_data[ 12 ] = T[0] * B[12] + T[4] * B[13] + T[8]  * B[14] + T[12] * B[15];
		  Base::mp_data[ 13 ] = T[1] * B[12] + T[5] * B[13] + T[9]  * B[14] + T[13] * B[15];
		  Base::mp_data[ 14 ] = T[2] * B[12] + T[6] * B[13] + T[10] * B[14] + T[14] * B[15];
		  Base::mp_data[ 15 ] = T[3] * B[12] + T[7] * B[13] + T[11] * B[14] + T[15] * B[15];
	  #undef T
	  #undef B
    }


    //! Multiply with another matrix from the left: this := A * this
    inline void multLeft( const MatrixT< Scalar, 4, 4 > & _lhs )
    {
	  #define A _lhs.mp_data
	  #define T cpy.Base::mp_data

		  Matrix4x4T cpy( *this );

		  Base::mp_data[ 0 ] = A[0] * T[0] + A[4] * T[1] + A[8]  * T[2] + A[12] * T[3];
		  Base::mp_data[ 1 ] = A[1] * T[0] + A[5] * T[1] + A[9]  * T[2] + A[13] * T[3];
		  Base::mp_data[ 2 ] = A[2] * T[0] + A[6] * T[1] + A[10] * T[2] + A[14] * T[3];
		  Base::mp_data[ 3 ] = A[3] * T[0] + A[7] * T[1] + A[11] * T[2] + A[15] * T[3];

		  Base::mp_data[ 4 ] = A[0] * T[4] + A[4] * T[5] + A[8]  * T[6] + A[12] * T[7];
		  Base::mp_data[ 5 ] = A[1] * T[4] + A[5] * T[5] + A[9]  * T[6] + A[13] * T[7];
		  Base::mp_data[ 6 ] = A[2] * T[4] + A[6] * T[5] + A[10] * T[6] + A[14] * T[7];
		  Base::mp_data[ 7 ] = A[3] * T[4] + A[7] * T[5] + A[11] * T[6] + A[15] * T[7];

		  Base::mp_data[ 8 ] = A[0] * T[8] + A[4] * T[9] + A[8]  * T[10] + A[12] * T[11];
		  Base::mp_data[ 9 ] = A[1] * T[8] + A[5] * T[9] + A[9]  * T[10] + A[13] * T[11];
		  Base::mp_data[10 ] = A[2] * T[8] + A[6] * T[9] + A[10] * T[10] + A[14] * T[11];
		  Base::mp_data[11 ] = A[3] * T[8] + A[7] * T[9] + A[11] * T[10] + A[15] * T[11];

		  Base::mp_data[ 12 ] = A[0] * T[12] + A[4] * T[13] + A[8]  * T[14] + A[12] * T[15];
		  Base::mp_data[ 13 ] = A[1] * T[12] + A[5] * T[13] + A[9]  * T[14] + A[13] * T[15];
		  Base::mp_data[ 14 ] = A[2] * T[12] + A[6] * T[13] + A[10] * T[14] + A[14] * T[15];
		  Base::mp_data[ 15 ] = A[3] * T[12] + A[7] * T[13] + A[11] * T[14] + A[15] * T[15];
	  #undef A
	  #undef T
    }


    //! Construct a transformation matrix B with P=(M|m) => PB=(Id|0)
    void setupTransformToIdentity( const ProjMatrixT< Scalar > & _matProj );

    //! Construct a transformation matrix B with P=(M|m) => (Id|0) * B = P
    void setupInvTransformToIdentity( const ProjMatrixT< Scalar > & _matProj );


    //! Construct rotation matrix from given parameters
    void rotation( const Vec3 & _axis, const Scalar & _angle );

    //! Construct rotation matrix from given parameters
    void rotation( const Vec3 & _axis );

    //! Construct translation matrix from given parameters
    void translation( const Vec3 & _trans );


};
    
}

}


#include "projmatrix.hh"


namespace Util {

namespace Math {


template< typename Scalar >
void Matrix4x4T< Scalar >::setupTransformToIdentity( const ProjMatrixT< Scalar > & _matProj )
{
    this->clear();
    copy( _matProj.m_inverse, *this, 0, 3, 0, 3 );
    (*this)( 0, 3 ) = _matProj.m_center[ 0 ];
    (*this)( 1, 3 ) = _matProj.m_center[ 1 ];
    (*this)( 2, 3 ) = _matProj.m_center[ 2 ];
    (*this)( 3, 3 ) = 1.0;
}



template< typename Scalar >
void Matrix4x4T< Scalar >::setupInvTransformToIdentity( const ProjMatrixT< Scalar > & _matProj )
{
    this->clear();
    copy( _matProj, *this, 0, 3, 0, 4 );
    (*this)( 3, 3 ) = 1.0;
}


template< typename Scalar >
void Matrix4x4T< Scalar >::rotation( const Vec3 & _axis, const Scalar & _angle )
{
    Scalar s = sin( _angle );
    Scalar c = cos( _angle );
    Scalar u = 1.0 - c;
    
    this->clear();
    
    Base::mp_data[ 0 ] = _axis[ 0 ] * _axis[ 0 ] * u + c;
    Base::mp_data[ 1 ] = _axis[ 0 ] * _axis[ 1 ] * u + _axis[ 2 ] * s;
    Base::mp_data[ 2 ] = _axis[ 0 ] * _axis[ 2 ] * u - _axis[ 1 ] * s;
    
    Base::mp_data[ 4 ] = _axis[ 0 ] * _axis[ 1 ] * u - _axis[ 2 ] * s;
    Base::mp_data[ 5 ] = _axis[ 1 ] * _axis[ 1 ] * u + c;
    Base::mp_data[ 6 ] = _axis[ 1 ] * _axis[ 2 ] * u + _axis[ 0 ] * s;
    
    Base::mp_data[ 8 ] = _axis[ 0 ] * _axis[ 2 ] * u + _axis[ 1 ] * s;
    Base::mp_data[ 9 ] = _axis[ 1 ] * _axis[ 2 ] * u - _axis[ 0 ] * s;
    Base::mp_data[ 10 ] = _axis[ 2 ] * _axis[ 2 ] * u + c;
    
    Base::mp_data[ 15 ] = 1.0;
}



template< typename Scalar >
void Matrix4x4T< Scalar >::rotation( const Vec3 & _axis )
{
    Scalar angle, scale;
    Vec3 axis( _axis );

    angle = _axis.length();
    
    if( angle < 1e-6 )
    {
	  this->clear();

	  (*this)( 0, 1 ) = -_axis[ 2 ];
	  (*this)( 0, 2 ) =  _axis[ 1 ];
	  (*this)( 1, 0 ) =  _axis[ 2 ];
	  (*this)( 1, 2 ) = -_axis[ 0 ];
	  (*this)( 2, 0 ) = -_axis[ 1 ];
	  (*this)( 2, 1 ) =  _axis[ 0 ];

	  Base::mp_data[ 0 ] = 1.0;
	  Base::mp_data[ 5 ] = 1.0;
	  Base::mp_data[ 10 ] = 1.0;

	  Base::mp_data[ 15 ] = 1.0;
    }
    else
    {
	  axis *= 1.0 / angle;
	  
	  rotation( axis, angle );
    }
}



template< typename Scalar >
void Matrix4x4T< Scalar >::translation( const Vec3 & _trans )
{
    this->clear();
    
    Base::mp_data[ 0 ] = 1.0;
    Base::mp_data[ 5 ] = 1.0;
    Base::mp_data[ 10 ] = 1.0;
    Base::mp_data[ 15 ] = 1.0;
    
    Base::mp_data[ 12 ] = _trans[ 0 ];
    Base::mp_data[ 13 ] = _trans[ 1 ];
    Base::mp_data[ 14 ] = _trans[ 2 ];
}



//! Type of a homography
typedef Matrix4x4T< double > Matrix4x4;
typedef Matrix4x4T< float > Matrix4x4f;

//! Type for a map of homographies
typedef std::map< int, Matrix4x4 > Matrix4x4Map;

//! Type for a vector of homographies
typedef std::vector< Matrix4x4 > Matrix4x4Vec;


}

}

#endif // _MATRIX_4X4_HH_
