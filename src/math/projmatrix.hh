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

#ifndef _PROJ_MATRIX_HH_
#define _PROJ_MATRIX_HH_

#include <iostream>
#include <vector>
#include <map>

#include "math.hh"
#include "matrixbase.hh"
#include "matrix3x3.hh"
#include "matrix4x4.hh"



extern "C" {
    int dgetri_( long int *n, double *a, long int *lda, long int *ipiv,
		 double *work, long int *lwork, long int *info );

    int dgetrf_( long int *m, long int *n, double *a, long int * lda,
		 long int *ipiv, long int *info );
}

namespace Util {

namespace Math {




//! Abstraction of a projection matrix
template< typename Scalar >
class ProjMatrixT : public MatrixT< Scalar, 3, 4 > {

public:

    //! Base type definition
    typedef MatrixT< Scalar, 3, 4 > Base;

    //! 2-dimensional vector
    typedef OpenMesh::VectorT< Scalar, 2 > Vec2;

    //! 2-dimensional vector
    typedef OpenMesh::VectorT< Scalar, 3 > Vec3;

    //! 3x3 matrix type
    typedef Matrix3x3T< Scalar > Matrix3x3;

    //! 4x4 matrix type
    typedef Matrix4x4T< Scalar > Matrix4x4;


    //! Default constructor
    ProjMatrixT()
	: Base()
    {
	  #ifdef COUNT_MATRIX_MEMORY
		  g_totalMatrixMemory += sizeof( Scalar ) * 3;
	  #endif
    }

    //! Copy constructor
    ProjMatrixT( const ProjMatrixT & _rhs )
	: Base( _rhs )
    {
	  m_inverse = _rhs.m_inverse;
	  m_center = _rhs.m_center;

	  #ifdef COUNT_MATRIX_MEMORY
		  g_totalMatrixMemory += sizeof( Scalar ) * 3;
	  #endif
    }

    //! Destructor
    ~ProjMatrixT()
    {
	  #ifdef COUNT_MATRIX_MEMORY
		  g_totalMatrixMemory -= sizeof( Scalar ) * 3;
	  #endif
    }
    

    //! Assignment operator
    const ProjMatrixT & operator=( const ProjMatrixT & _rhs )
    {
	  Base::operator=( _rhs );

	  m_inverse = _rhs.m_inverse;
	  m_center = _rhs.m_center;

	  return *this;
    }


    //! Set matrix to identity
    inline void setIdentity( void )
    {
	  memset( Base::mp_data, 0, 12 * sizeof( Scalar ) );
	  Base::mp_data[ 0 ] = Base::mp_data[ 4 ] = Base::mp_data[ 8 ] = (Scalar) 1.0;

	  m_inverse.setIdentity();

	  m_center[ 0 ] = 0.0;
	  m_center[ 1 ] = 0.0;
	  m_center[ 2 ] = 0.0;
    }


    //! Assign given matrix to the left 3x3 block
    inline void setLeftBlock( const Matrix3x3 & _block )
    {
	  memcpy( Base::mp_data, _block.mp_data, 9 * sizeof( Scalar ) );
    }

    //! Access to the left 3x3 block of matrix
    inline void getLeftBlock( const Matrix3x3 & _block ) const
    {
	  memcpy( _block.mp_data, Base::mp_data, 9 * sizeof( Scalar ) );
    }


    //! Assign given vector to the right-most column of the matrix
    inline void setRightColumn( const Vec3 & _col )
    {
	  Base::mp_data[ 9 ] = _col[ 0 ];
	  Base::mp_data[ 10 ] = _col[ 1 ];
	  Base::mp_data[ 11 ] = _col[ 2 ];
    }

    //! Access to the right-most column of the matrix
    inline void getRightColumn( Vec3 & _col ) const
    {
	  _col[ 0 ] = Base::mp_data[ 9 ];
	  _col[ 1 ] = Base::mp_data[ 10 ];
	  _col[ 2 ] = Base::mp_data[ 11 ];
    }


    //! Compute inverse of left 3x3 block
    bool computeInverse( void )
    {
	  // copy left 3x3 block
	  memcpy( m_inverse.mp_data, Base::mp_data, 9 * sizeof( Scalar ) );

	  return m_inverse.invert();
    }



    //! Compute camera center
    void computeCenter( void )
    {
	  OpenMesh::VectorT< Scalar, 3 > hlp;
	  
	  hlp[ 0 ] = -Base::mp_data[ 9 ];
	  hlp[ 1 ] = -Base::mp_data[ 10 ];
	  hlp[ 2 ] = -Base::mp_data[ 11 ];

	  m_inverse.multRight( hlp, m_center );
    }



#define PROJECTION( SX, SY, SZ ) \
	const Scalar num1  = Base::mp_data[0] * (SX) + Base::mp_data[3] * (SY) + Base::mp_data[6] * (SZ) + Base::mp_data[9]; \
	const Scalar num2  = Base::mp_data[1] * (SX) + Base::mp_data[4] * (SY) + Base::mp_data[7] * (SZ) + Base::mp_data[10]; \
	const Scalar denom = Base::mp_data[2] * (SX) + Base::mp_data[5] * (SY) + Base::mp_data[8] * (SZ) + Base::mp_data[11];



    //! Multiply with a 3D double point
    inline void project( const Vec3 & _a, Vec2 & _b ) const
    {
	  PROJECTION( _a[ 0 ], _a[ 1 ], _a[ 2 ] );

	  _b[ 0 ] = num1 / denom;
	  _b[ 1 ] = num2 / denom;
    }


    //! Multiply with a 3D double point
    inline void project( const Vec3 & _a, Scalar & _x, Scalar & _y ) const
    {
	  PROJECTION( _a[ 0 ], _a[ 1 ], _a[ 2 ] );

	  _x = num1 / denom;
	  _y = num2 / denom;
    }


    //! Multiply with a 3D double point, without de-homogenization
    inline void project( const Vec3 & _a, Vec3 & _b ) const
    {
	  _b[0] = Base::mp_data[0] * _a[0] + Base::mp_data[3] * _a[1] 
		  + Base::mp_data[6] * _a[2] + Base::mp_data[9];
	  _b[1] = Base::mp_data[1] * _a[0] + Base::mp_data[4] * _a[1] 
		  + Base::mp_data[7] * _a[2] + Base::mp_data[10];
	  _b[2] = Base::mp_data[2] * _a[0] + Base::mp_data[5] * _a[1] 
		  + Base::mp_data[8] * _a[2] + Base::mp_data[11];
    }


    //! Multiply with a 3D double point, without de-homogenization
    inline void project( const Vec3 & _a, Scalar & _x, Scalar & _y, Scalar & _z ) const
    {
	  _x = Base::mp_data[0] * _a[0] + Base::mp_data[3] * _a[1]
		  + Base::mp_data[6] * _a[2] + Base::mp_data[9];
	  _y = Base::mp_data[1] * _a[0] + Base::mp_data[4] * _a[1]
		  + Base::mp_data[7] * _a[2] + Base::mp_data[10];
	  _z = Base::mp_data[2] * _a[0] + Base::mp_data[5] * _a[1]
		  + Base::mp_data[8] * _a[2] + Base::mp_data[11];
    }


    //! Multiply with a 3D double point
    inline void project( const Vec3 & _a, OpenMesh::Vec2i & _b ) const
    {
	  PROJECTION( _a[ 0 ], _a[ 1 ], _a[ 2 ] );

	  _b[ 0 ] = (int) round( num1 / denom );
	  _b[ 1 ] = (int) round( num2 / denom );
    }


    //! Multiply with a 3D double point
    inline void project( const Vec3 & _a, int & _x, int & _y ) const
    {
	  PROJECTION( _a[ 0 ], _a[ 1 ], _a[ 2 ] );

	  _x = (int) round( num1 / denom );
	  _y = (int) round( num2 / denom );
    }


    //! Multiply with a 3D double point, without de-homogenization
    inline void project( const OpenMesh::Vec4d & _a, Vec2 & _b ) const
    {
	  _b[ 0 ] = Base::mp_data[0] * _a[0] + Base::mp_data[3] * _a[1] + Base::mp_data[6] * _a[2] + Base::mp_data[9]  * _a[3];
	  _b[ 1 ] = Base::mp_data[1] * _a[0] + Base::mp_data[4] * _a[1] + Base::mp_data[7] * _a[2] + Base::mp_data[10] * _a[3];
	  const Scalar denom = Base::mp_data[2] * _a[0] + Base::mp_data[5] * _a[1] + Base::mp_data[8] * _a[2] + Base::mp_data[11] * _a[3];

	  _b[ 0 ] /= denom;
	  _b[ 1 ] /= denom;
    }


    //! Multiply with a 3D double point, without de-homogenization
    inline void project( const OpenMesh::Vec4d & _a, Scalar & _x, Scalar & _y ) const
    {
	  _x = Base::mp_data[0] * _a[0] + Base::mp_data[3] * _a[1] + Base::mp_data[6] * _a[2] + Base::mp_data[9]  * _a[3];
	  _y = Base::mp_data[1] * _a[0] + Base::mp_data[4] * _a[1] + Base::mp_data[7] * _a[2] + Base::mp_data[10]  * _a[3];
	  const Scalar denom = Base::mp_data[2] * _a[0] + Base::mp_data[5] * _a[1] + Base::mp_data[8] * _a[2] + Base::mp_data[11]  * _a[3];

	  _x /= denom;
	  _y /= denom;
    }


    //! Multiply with a 3D double point, without de-homogenization
    inline void project( const OpenMesh::Vec4d & _a, Scalar & _x, Scalar & _y, Scalar & _z ) const
    {
	  _x = Base::mp_data[0] * _a[0] + Base::mp_data[3] * _a[1] + Base::mp_data[6] * _a[2] + Base::mp_data[9]  * _a[3];
	  _y = Base::mp_data[1] * _a[0] + Base::mp_data[4] * _a[1] + Base::mp_data[7] * _a[2] + Base::mp_data[10]  * _a[3];
	  _z = Base::mp_data[2] * _a[0] + Base::mp_data[5] * _a[1] + Base::mp_data[8] * _a[2] + Base::mp_data[11]  * _a[3];
    }

#undef PROJECTION

#define UNPROJECTION( SX, SY, DX, DY, DZ ) \
        (DX) = m_inverse.mp_data[0] * (SX) + m_inverse.mp_data[3] * (SY) + m_inverse.mp_data[6]; \
	(DY) = m_inverse.mp_data[1] * (SX) + m_inverse.mp_data[4] * (SY) + m_inverse.mp_data[7]; \
	(DZ) = m_inverse.mp_data[2] * (SX) + m_inverse.mp_data[5] * (SY) + m_inverse.mp_data[8];


    //! Unproject given image point
    inline void unproject( const Vec2 & _a, Vec3 & _b ) const
    {
	  UNPROJECTION( _a[ 0 ], _a[ 1 ], _b[ 0 ], _b[ 1 ], _b[ 2 ] );
    }

    //! Unproject given image point
    inline void unproject( const OpenMesh::Vec2i & _a, Vec3 & _b ) const
    {
	  UNPROJECTION( _a[ 0 ], _a[ 1 ], _b[ 0 ], _b[ 1 ], _b[ 2 ] );
    }


    //! Unproject given image point
    template< typename ScalarB >
    inline void unproject( const ScalarB & _x, const ScalarB & _y, Vec3 & _b ) const
    {
	  UNPROJECTION( (Scalar)_x, (Scalar)_y, _b[ 0 ], _b[ 1 ], _b[ 2 ] );
    }

#undef UNPROJECTION



    //! Multiply with a 4x4 matrix from the right: this := this * B
    inline void multRight( const Matrix4x4 & _rhs )
    {
	  #define T cpy.mp_data
	  #define B _rhs.mp_data

		  ProjMatrixT cpy( *this );

		  Base::mp_data[ 0 ] = T[0] * B[0] + T[3] * B[1] + T[6] * B[2] + T[9]  * B[3];
		  Base::mp_data[ 1 ] = T[1] * B[0] + T[4] * B[1] + T[7] * B[2] + T[10] * B[3];
		  Base::mp_data[ 2 ] = T[2] * B[0] + T[5] * B[1] + T[8] * B[2] + T[11] * B[3];

		  Base::mp_data[ 3 ] = T[0] * B[4] + T[3] * B[5] + T[6] * B[6] + T[9]  * B[7];
		  Base::mp_data[ 4 ] = T[1] * B[4] + T[4] * B[5] + T[7] * B[6] + T[10] * B[7];
		  Base::mp_data[ 5 ] = T[2] * B[4] + T[5] * B[5] + T[8] * B[6] + T[11] * B[7];

		  Base::mp_data[ 6 ] = T[0] * B[8] + T[3] * B[9] + T[6] * B[10] + T[9]  * B[11];
		  Base::mp_data[ 7 ] = T[1] * B[8] + T[4] * B[9] + T[7] * B[10] + T[10] * B[11];
		  Base::mp_data[ 8 ] = T[2] * B[8] + T[5] * B[9] + T[8] * B[10] + T[11] * B[11];

		  Base::mp_data[ 9 ] = T[0] * B[12] + T[3] * B[13] + T[6] * B[14] + T[9]  * B[15];
		  Base::mp_data[10 ] = T[1] * B[12] + T[4] * B[13] + T[7] * B[14] + T[10] * B[15];
		  Base::mp_data[11 ] = T[2] * B[12] + T[5] * B[13] + T[8] * B[14] + T[11] * B[15];
	  #undef T
	  #undef B
    }


    //! Multiply with a 3x3 matrix from the left: this := A * this
    inline void multLeft( const Matrix3x3 & _lhs )
    {
	  #define A _lhs.mp_data
	  #define T cpy.mp_data

		  ProjMatrixT cpy( *this );

		  Base::mp_data[ 0 ] = A[0] * T[0] + A[3] * T[1] + A[6] * T[2];
		  Base::mp_data[ 1 ] = A[1] * T[0] + A[4] * T[1] + A[7] * T[2];
		  Base::mp_data[ 2 ] = A[2] * T[0] + A[5] * T[1] + A[8] * T[2];

		  Base::mp_data[ 3 ] = A[0] * T[3] + A[3] * T[4] + A[6] * T[5];
		  Base::mp_data[ 4 ] = A[1] * T[3] + A[4] * T[4] + A[7] * T[5];
		  Base::mp_data[ 5 ] = A[2] * T[3] + A[5] * T[4] + A[8] * T[5];

		  Base::mp_data[ 6 ] = A[0] * T[6] + A[3] * T[7] + A[6] * T[8];
		  Base::mp_data[ 7 ] = A[1] * T[6] + A[4] * T[7] + A[7] * T[8];
		  Base::mp_data[ 8 ] = A[2] * T[6] + A[5] * T[7] + A[8] * T[8];

		  Base::mp_data[ 9  ] = A[0] * T[9] + A[3] * T[10] + A[6] * T[11];
		  Base::mp_data[ 10 ] = A[1] * T[9] + A[4] * T[10] + A[7] * T[11];
		  Base::mp_data[ 11 ] = A[2] * T[9] + A[5] * T[10] + A[8] * T[11];
	  #undef A
	  #undef T
    }


    //! Decompose left 3x3 block into a triangular matrix K and a rotation matrix R
    void decompose( Matrix3x3 & _matK, Matrix3x3 & _matR ) const
    {
	  Scalar dummy;
	  decompose( _matK, _matR, dummy );
    }

    //! Decompose left 3x3 block into a triangular matrix K and a rotation matrix R
    void decompose( Matrix3x3 & _matR ) const
    {
	  Matrix3x3 matK;
	  Scalar dummy;
	  decompose( matK, _matR, dummy );
    }


    void compose( const Matrix3x3 & _matK,
		  const Matrix3x3 & _matR )
    {
	  static Vec3 vecM;
	  static Matrix3x3 matM;

	  matM = _matK;
	  matM.multRight( _matR );

	  memcpy( this->mp_data, matM.mp_data, 9 * sizeof( Scalar ) );

	  matM.multRight( m_center, vecM );
	  Base::mp_data[ 9  ] = -vecM[ 0 ];
	  Base::mp_data[ 10 ] = -vecM[ 1 ];
	  Base::mp_data[ 11 ] = -vecM[ 2 ];
    }


    //! Normalize scale of projection matrix
    void normalizeScale( void )
    {
		Matrix3x3 matK, matR;
		Scalar scale;

		decompose( matK, matR, scale );

		this->scale( scale );
    }
    
    //! Get position of camera and its orientation (see Hartley & Zisserman, 2nd edition, pages 158-161), returns false if one of the could not be computed
    bool getPositionAndOrientation( OpenMesh::Vec3d &position, OpenMesh::Vec3d &orientation )
    {
      // get the position of the camera as the vector spanning the right null-space of the 
      // projection matrix, see Hartley & Zisserman, 2nd edition, pages 158-159
      Matrix mat_VT( 4, 4 ), mat_P( 3, 4 );
      for( int i=0; i<4; ++i )
      {
		for( int j=0; j<3; ++j )
		  mat_P(j,i) = Base::mp_data[ 3 * i + j ];
      }
      
      if( ! SVD_V( mat_P, mat_VT ) )
		return false;
      
    
      for( int i=0; i<3; ++i )
		position[i] = mat_VT( 3,i ) / mat_VT( 3,3 );
      
      // get the viewing direction of the camera, see Hartley & Zisserman, 2nd ed., pages 160 - 161
      // compute the determinant of the 3x3 part of the projection matrix
      double det = Base::mp_data[0] * Base::mp_data[4] * Base::mp_data[8] + Base::mp_data[3] * Base::mp_data[7] * Base::mp_data[2] + Base::mp_data[6] * Base::mp_data[1] * Base::mp_data[5] - Base::mp_data[6] * Base::mp_data[4] * Base::mp_data[2] - Base::mp_data[3] * Base::mp_data[1] * Base::mp_data[8] - Base::mp_data[0] * Base::mp_data[7] * Base::mp_data[5];
      
	  // remember that the camera in reconstructions computed by Bundler looks 
	  // down the negative z-axis instead of the positive z-axis.
	  // So we have to multiply the orientation with -1.0
      orientation[0] = - Base::mp_data[2] * det;
      orientation[1] = - Base::mp_data[5] * det;
      orientation[2] = - Base::mp_data[8] * det;
      
      orientation.normalize();
      
      return true;
    }



    //! Storage for the inverse left 3x3 block, column-major
    Matrix3x3  m_inverse;

    //! Storage for camera center
    OpenMesh::VectorT< Scalar, 3 > m_center;



protected:

    //! Decompose left 3x3 block into a triangular matrix K and a rotation matrix R
    void decompose( Matrix3x3 & _matK, Matrix3x3 & _matR,
		    Scalar & _scale ) const
    {
	  Matrix3x3 matRot;

	  ////////////////////////////////////
	  // Algorithm taken from 
	  // Hartley, Zisserman: "Multiple View Geometry in Computer Vision",
	  // Appendix 4.1


	  Scalar s, c;

	  ////////////////////////////////////
	  // Cancellation of element 3,2

	  copy( *this, _matK, 0, 3, 0, 3 );

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
	  _scale = 1.0 / _matK( 2, 2 );
	  if( !isinf( _scale ) && !isnan( _scale ) )
		  _matK.scale( _scale );
    }

};
    



//! Type of a projection matrix
typedef ProjMatrixT< double > ProjMatrix;

//! Type for a map of projection matrices
typedef std::map< int, ProjMatrix > ProjMatrixMap;

//! Type for a vector of projection matrices
typedef std::vector< ProjMatrix > ProjMatrixVec;






//! Abstraction of a decomposed projection matrix
template< typename Scalar >
class ProjMatrixDecompT : public ProjMatrixT< Scalar > {


public:

    //! Base type definition
    typedef ProjMatrixT< Scalar > Base;

    //! 2-dimensional vector
    typedef OpenMesh::VectorT< Scalar, 2 > Vec2;

    //! 2-dimensional vector
    typedef OpenMesh::VectorT< Scalar, 3 > Vec3;

    //! 3x3 matrix type
    typedef Matrix3x3T< Scalar > Matrix3x3;

    //! 4x4 matrix type
    typedef Matrix4x4T< Scalar > Matrix4x4;


    //! Default constructor
    ProjMatrixDecompT()
	: Base()
    {
	  #ifdef COUNT_MATRIX_MEMORY
		  g_totalMatrixMemory += sizeof( Scalar ) * 3;
	  #endif
    }

    //! Copy constructor
    ProjMatrixDecompT( const ProjMatrixDecompT & _rhs )
	: Base( _rhs ),
	  m_matrixK( _rhs.m_matrixK ),
	  m_matrixKInv( _rhs.m_matrixKInv ),
	  m_matrixR( _rhs.m_matrixR ),
	  m_vecR( _rhs.m_vecR )
    {
	  #ifdef COUNT_MATRIX_MEMORY
		  g_totalMatrixMemory += sizeof( Scalar ) * 3;
	  #endif
    }

    //! Initializing constructor
    ProjMatrixDecompT( const Base & _rhs )
	: Base( _rhs )
    {
	  #ifdef COUNT_MATRIX_MEMORY
		  g_totalMatrixMemory += sizeof( Scalar ) * 3;
	  #endif
    }

    //! Destructor
    ~ProjMatrixDecompT()
    {
	  #ifdef COUNT_MATRIX_MEMORY
		  g_totalMatrixMemory -= sizeof( Scalar ) * 3;
	  #endif
    }

    ProjMatrixDecompT< Scalar > & operator=( const ProjMatrixDecompT< Scalar > & _rhs )
    {
	  Base::operator=( _rhs );

	  m_matrixK = _rhs.m_matrixK;
	  m_matrixKInv = _rhs.m_matrixKInv;
	  m_matrixR = _rhs.m_matrixR;
	  m_vecR = _rhs.m_vecR;

	  return *this;
    }

    ProjMatrixDecompT< Scalar > & operator=( const Base & _rhs )
    {
	  Base::operator=( _rhs );
	  return *this;
    }


    //! Decomposition of current matrix
    inline void decompose( void )
    {
	  Base::decompose( m_matrixK, m_matrixR );
	  
	  m_matrixR.decomposeRotation( m_vecR );

	  m_matrixKInv = m_matrixK;
	  m_matrixKInv.invert();
    }

    //! Composition from current intrinsic and extrinsic matrices
    inline void compose( void )
    {
	  m_matrixR.composeRotation( m_vecR );
	  
	  Base::compose( m_matrixK, m_matrixR );
    }
    


    //! Matrix of intrinsic parameters
    Matrix3x3 m_matrixK;

    //! Matrix of intrinsic parameters
    Matrix3x3 m_matrixKInv;

    //! Matrix of intrinsic parameters
    Matrix3x3 m_matrixR;

    //! Rotation matrix
    Vec3 m_vecR;

};



//! Type of a decomposed projection matrix
typedef ProjMatrixDecompT< double > ProjMatrixDecomp;

//! Type for a map of projection matrices
typedef std::map< int, ProjMatrixDecomp > ProjMatrixDecompMap;

//! Type for a vector of projection matrices
typedef std::vector< ProjMatrixDecomp > ProjMatrixDecompVec;


}

}


#endif // _PROJ_MATRIX_HH_
