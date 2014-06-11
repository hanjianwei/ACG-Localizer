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

#ifndef _MATRIX_BASE_HH_
#define _MATRIX_BASE_HH_

#include <iostream>
#include <vector>
#include <map>

#include "math.hh"

#include <OpenMesh/Core/Geometry/VectorT.hh>


#define COUNT_MATRIX_MEMORY



extern "C" {

    int dgesvd_( char *jobu, char *jobvt, long int *m, long int *n,
		 double *a, long int *lda, double *s, double *u, long int *ldu,
		 double *vt, long int *ldvt, double *work, long int *lwork,
		 long int *info );

    int dsyev_( char *jobz, char *uplo, long int *n, double *a,
		long int *lda, double *w, double *work, long int *lwork, 
		long int *info);

    int dgeev_( char *jobvl, char *jobvr, long int *n, double *a,
		long int *lda, double *wr, double *wi, double *vl,
		long int *ldvl, double *vr, long int *ldvr, double *work,
		long int *lwork, long int *info );

    int dpotrf_( char *uplo, long int *n, double *a, long int *lda,
		 long int *info);
}


namespace Util {

namespace Math {


#ifdef COUNT_MATRIX_MEMORY
extern long int g_totalMatrixMemory;
#endif


//! Abstraction of an image homography
template< typename Scalar, int ROWS, int COLS >
class MatrixT {

public:

    //! Row vector type
    typedef OpenMesh::VectorT< Scalar, COLS > VecRow;

    //! Column vector type
    typedef OpenMesh::VectorT< Scalar, ROWS > VecCol;


    //! Default constructor
    MatrixT()
    {
	  mp_data = new Scalar[ ROWS * COLS ];
	  
  #ifdef COUNT_MATRIX_MEMORY
	  g_totalMatrixMemory += sizeof( Scalar ) * ROWS * COLS;
  #endif
    }

    //! Copy constructor
    MatrixT( const MatrixT & _rhs )
    {
	  mp_data = new Scalar[ ROWS * COLS ];
	  memcpy( mp_data, _rhs.mp_data, ROWS * COLS * sizeof( Scalar ) );

  #ifdef COUNT_MATRIX_MEMORY
	  g_totalMatrixMemory += sizeof( Scalar ) * ROWS * COLS;
  #endif
    }

    //! Destructor
    ~MatrixT()
    {
	  delete [] mp_data;
	  mp_data = 0;

  #ifdef COUNT_MATRIX_MEMORY
	  g_totalMatrixMemory -= sizeof( Scalar ) * ROWS * COLS;
  #endif
    }
    

    //! Assignment operator
    inline const MatrixT & operator=( const MatrixT & _rhs )
    {
	  if( &_rhs == this )
		  return *this;

	  memcpy( mp_data, _rhs.mp_data, ROWS * COLS * sizeof( Scalar ) );
	  return *this;
    }


    //! Read access to the elements
    inline Scalar operator()( const int _row, const int _col ) const
    {
	  return mp_data[ ROWS * _col + _row ];
    }

    //! Write access to the elements
    inline Scalar & operator()( const int _row, const int _col )
    {
	  return mp_data[ ROWS * _col + _row ];
    }

    //! Add two matrices
    inline const MatrixT & operator+=( const MatrixT & _rhs )
    {
	  for( unsigned int x = 0; x < COLS * ROWS; ++x )
		  mp_data[ x ] += _rhs.mp_data[ x ];

	  return *this;
    }

    //! Subtract two matrices
    inline const MatrixT & operator-=( const MatrixT & _rhs )
    {
	  for( unsigned int x = 0; x < COLS * ROWS; ++x )
		  mp_data[ x ] -= _rhs.mp_data[ x ];

	  return *this;
    }


    //! Set all elements to zero
    inline void clear( void )
    {
	  memset( mp_data, 0, ROWS * COLS * sizeof( Scalar ) );
    }

    
    //! Scale all elements with given factor
    inline void scale( const Scalar & _scale )
    {
	  const int endI = ROWS * COLS;
	  for( int i = 0; i < endI; ++i )
		  mp_data[ i ] *= _scale;
    }


    //! Add two matrices
    inline void add( const MatrixT & _mat )
    {
	  for( unsigned int x = 0; x < COLS * ROWS; ++x )
		  mp_data[ x ] += _mat.mp_data[ x ];
    }

    //! Add two matrices
    inline void subtract( const MatrixT & _mat )
    {
	  int x;
	  for( x = 0; x < COLS * ROWS; ++x )
		  mp_data[ x ] -= _mat.mp_data[ x ];
    }


    //! Compute Frobenius norm of matrix
    inline double frobeniusNormSqr( void ) const
    {
	  Scalar norm = 0.0;
	  int x;
	  for( x = 0; x < COLS * ROWS; ++x )
		  norm += mp_data[ x ] * mp_data[ x ];

	  return norm;
    }


    //! Copy to our GMM matrix type
    void copyTo( Util::Math::Matrix & _matrix ) const
    {
	  int x, y;
	  for( x = 0; x < COLS; ++x )
	  {
		  for( y = 0; y < ROWS; ++y )
		  {
			_matrix( y, x ) = mp_data[ y + ROWS*x ];
		  }
	  }
    }

    //! Copy from our GMM matrix type
    void copyFrom( const Util::Math::Matrix & _matrix )
    {
	  int x, y;

	  for( x = 0; x < COLS; ++x )
	  {
		  for( y = 0; y < ROWS; ++y )
		  {
			mp_data[ y + ROWS*x ] = _matrix( y, x );
		  }
	  }
    }

    
    //! Construct transposed matrix
    const MatrixT transposed( void ) const
    {
	  MatrixT< Scalar, COLS, ROWS > matT;
	  int x, y;

	  for( y = 0; y < ROWS; ++y )
	  {
		  for( x = 0; x < COLS; ++x )
		  {
			matT( y, x ) = (*this)( x, y );
		  }
	  }
	  
	  return matT;
    }



    //! Access to the columns of matrix
    void getColumn( const int _col, VecCol & _vec ) const
    {
	  for( int y = 0; y < ROWS; ++y )
		  _vec[ y ] = mp_data[ y + ROWS*_col ];
    }

    //! Access to the rows of matrix
    void getRow( const int _row, VecRow & _vec ) const
    {
	  for( int x = 0; x < COLS; ++x )
		  _vec[ x ] = mp_data[ _row + ROWS*x ];
    }


    //! Access to the columns of matrix
    void setColumn( const int _col, const VecCol & _vec )
    {
	  for( int y = 0; y < ROWS; ++y )
		  mp_data[ y + ROWS*_col ] = _vec[ y ];
    }

    //! Access to the rows of matrix
    void setRow( const int _row, const VecRow & _vec )
    {
	  for( int x = 0; x < COLS; ++x )
		  mp_data[ _row + ROWS*x ] = _vec[ x ];
    }


    //! Storage for matrix, column-major
    Scalar * mp_data;
};




template< typename Scalar, int ROWS, int COLS >
std::ostream & operator<<( std::ostream & _os, const MatrixT< Scalar, ROWS, COLS > & _mat )
{
    int r, c;

    for( r = 0; r < ROWS; ++r )
    {
	  _os << "[ ";
	  for( c = 0; c < COLS; ++c )
	  {
		  _os << _mat( r, c );
		  if( c < COLS - 1 )
		  _os << ", ";
	  }

	  _os << " ]" << std::endl;
    }

    return _os;
}



template< typename Scalar, int ROWS, int COLS >
const MatrixT< Scalar, ROWS, COLS >
operator+( const MatrixT< Scalar, ROWS, COLS > & _lhs,
	   const MatrixT< Scalar, ROWS, COLS > & _rhs )
{
    MatrixT< Scalar, ROWS, COLS > result( _lhs );

    for( int x = 0; x < COLS * ROWS; ++x )
	result.mp_data[ x ] += _rhs.mp_data[ x ];

    return result;
}



template< typename Scalar, int ROWS, int COLS >
const MatrixT< Scalar, ROWS, COLS >
operator-( const MatrixT< Scalar, ROWS, COLS > & _lhs,
	   const MatrixT< Scalar, ROWS, COLS > & _rhs )
{
    MatrixT< Scalar, ROWS, COLS > result( _lhs );

    for( int x = 0; x < COLS * ROWS; ++x )
	result.mp_data[ x ] -= _rhs.mp_data[ x ];

    return result;
}





/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////



template< typename Scalar, int DIM >
class MatrixSqrT : public MatrixT< Scalar, DIM, DIM > {

public:

    //! Base type definition
    typedef MatrixT< Scalar, DIM, DIM > Base;

    //! Default constructor
    MatrixSqrT()
	: Base()
    {
    }

    //! Copy constructor
    MatrixSqrT( const MatrixT< Scalar, DIM, DIM > & _rhs )
	: Base( _rhs )
    {
    }


    //! Set matrix to identity
    inline void setIdentity( void )
    {
	  memset( Base::mp_data, 0, DIM * DIM * sizeof( Scalar ) );
	  for( int i = 0; i < DIM; ++i )
		  Base::mp_data[ i * DIM + i ] = (Scalar) 1.0;
    }

    //! Invert matrix
    bool invert( void )
    {
	  static double work[ 100 ];
	  static long int lwork = 100;
	  static long int ipiv[ DIM ];
	  static long int info = 0;
	  static long int n = DIM;

	  dgetrf_( &n, &n, Base::mp_data, &n, ipiv, &info );

	  if( info != 0 )
		  return false;

	  dgetri_( &n, Base::mp_data, &n, ipiv,
		  work, &lwork, &info );

	  return (info == 0);
    }

    //! Transpose matrix
    void transpose( void )
    {
	  int x, y;
	  Scalar hlp;

	  for( y = 0; y < DIM-1; ++y ) // last row can be neglected
	  {
		  for( x = y + 1; x < DIM; ++x )
		  {
		  hlp = Base::mp_data[ DIM * x + y ];
		  Base::mp_data[ DIM * x + y ] = Base::mp_data[ DIM * y + x ];
		  Base::mp_data[ DIM * y + x ] = hlp;
		  }
	  }
    }

    //! Create a new, transposed matrix
    const MatrixSqrT transposed( void ) const
    {
	  MatrixSqrT transp( *this );
	  transp.transpose();
	  return transp;
    }

};





//! Copy a sub-matrix from _src to _dst
template< typename Scalar, int ROWS1, int COLS1, int ROWS2, int COLS2 >
void copy( const MatrixT< Scalar, ROWS1, COLS1 > & _src,
	   MatrixT< Scalar, ROWS2, COLS2 > & _dst,
	   const int _sourceRow, const int _numRows,
	   const int _sourceCol, const int _numCols,
	   const int _destRow = 0, const int _destCol = 0 )
{
    int sx, sy, dx, dy;
    
    const int esx = _sourceCol + _numCols;
    const int esy = _sourceRow + _numRows;
    
    dy = _destRow;
    for( sy = _sourceRow; sy < esy; ++sy, ++dy )
    {
	  dx = _destCol;
	  for( sx = _sourceCol; sx < esx; ++sx, ++dx )
	  {
		  _dst( dy, dx ) = _src( sy, sx );
	  }
    }
}


}

}

#endif // _MATRIX_BASE_HH_
