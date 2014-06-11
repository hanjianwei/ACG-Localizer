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
 * Definition of Solver for the Projection Matrix
 *
 * Author: Martin Habbecke <M.Habbecke@gmx.de>
 * Version: Revision: 1.12
 * Date: 2011-10-03 17:31
 */


#include <vector>
#include <float.h>
#include <math.h>

#include "solverproj.hh"



using namespace Util::Math;

namespace Util {
namespace CorrSolver {



/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////


//! Global pointer to correspondences used by fcn_dist()
static SolverProj::CorrespondenceVec * gp_corresp;




/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////



SolverProj::SolverProj()
{
    m_totalNumIterations = 0;
}


SolverProj::~SolverProj()
{
}



void SolverProj::clear( void )
{
    mv_corresp.clear();
    mv_correspScaled.clear();
}


void SolverProj::addCorrespondence( const Vector2D & _point2D,
				    const Vector3D & _point3D )
{
    // add correspondence
    mv_corresp.push_back( Correspondence( _point2D, _point3D ) );
}


void SolverProj::addCorrespondence( const Correspondence & _corresp )
{
    // add correspondence
    mv_corresp.push_back( _corresp );
}




bool SolverProj::scaleCorrespondences( void )
{
    CorrespondenceVec::iterator cIter, cEnd;
    cEnd = mv_corresp.end();

    Vector2D pt1;
    Vector3D pt2;

    double n;

    Vector2D center2D( 0.0, 0.0 );
    Vector3D center3D( 0.0, 0.0, 0.0 );
    double scale1, scale2;


    //////////////////////////////////////////////
    // determine c.o.g. in images and (projective) world space

    n = 1.0;

    for( cIter = mv_corresp.begin();
	 cIter != cEnd;
	 ++cIter )
    {
	  center2D = cIter->m_point2D * (1.0/n) + center2D * ((n-1.0) / n);
	  center3D = cIter->m_point3D * (1.0/n) + center3D * ((n-1.0) / n);

	  n += 1.0;
    }


    ///////////////////////////////////////////////
    // move features to COG and compute scaling factors

    mv_correspScaled.clear();

    n = 1.0;
    scale1 = 0.0;
    scale2 = 0.0;

    for( cIter = mv_corresp.begin();
	 cIter != cEnd;
	 ++cIter )
    {
	  pt1 = cIter->m_point2D - center2D;
	  pt2 = cIter->m_point3D - center3D;

	  scale1 = scale1 * ((n-1.0) / n) + pt1.norm() * (1.0 / n);
	  scale2 = scale2 * ((n-1.0) / n) + pt2.norm() * (1.0 / n);
	  
	  mv_correspScaled.push_back( Correspondence( pt1, pt2 ) );

	  n += 1.0;
    }
    
    if( fabs(scale1) < 1e-12 || fabs(scale2) < 1e-12 )
      return false;

    scale1 = 1.41421 / scale1; // sqrt(2)
    scale2 = 1.73205 / scale2; // sqrt(3)


    //////////////////////////////////////////////
    // scale coords in both images

    cEnd = mv_correspScaled.end();
    for( cIter = mv_correspScaled.begin();
	 cIter != cEnd;
	 ++cIter )
    {
	  cIter->m_point2D *= scale1;
	  cIter->m_point3D *= scale2;
    }


    ////////////////////////////////////////////
    // create scaling matrices that will be multiplied
    // to result matrix

    m_matScaleImgInv.clear();
    m_matScaleImgInv( 0, 0 ) = 1.0 / scale1;
    m_matScaleImgInv( 1, 1 ) = 1.0 / scale1;
    m_matScaleImgInv( 0, 2 ) = center2D[ 0 ];
    m_matScaleImgInv( 1, 2 ) = center2D[ 1 ];
    m_matScaleImgInv( 2, 2 ) = 1.0;
    
    m_matScaleWorld.clear();
    m_matScaleWorld( 0, 0 ) = scale2;
    m_matScaleWorld( 1, 1 ) = scale2;
    m_matScaleWorld( 2, 2 ) = scale2;
    m_matScaleWorld( 0, 3 ) = -center3D[ 0 ]*scale2;
    m_matScaleWorld( 1, 3 ) = -center3D[ 1 ]*scale2;
    m_matScaleWorld( 2, 3 ) = -center3D[ 2 ]*scale2;
    m_matScaleWorld( 3, 3 ) = 1.0;
    
    return true;
}



void SolverProj::setInitialMatrix( const Matrix & _mat )
{
    m_initialMatrix.copyFrom( _mat );
}


void SolverProj::setInitialMatrix( const ProjMatrix & _mat )
{
    m_initialMatrix = _mat;
}


bool SolverProj::computeLinearNoRescaling( void )
{
    int nrows, ncols;

    nrows = 2 * mv_corresp.size();
    ncols = 12;

    Matrix mat_A( nrows, ncols );
    Matrix mat_VT( ncols, ncols );

    gmm::clear( mat_A );

    ///////////////////////////////////////////////////
    // coordinate scaling

    if( !scaleCorrespondences() )
      return false;

    int corr, endCorr = mv_correspScaled.size();

    ///////////////////////////////////////////////////
    // generate matrix A, see Hartley & Zisserman, 2nd edition, page 178-179

    for( corr = 0; corr < endCorr; ++corr )
    {
	  mat_A( 2*corr, 4 ) = - mv_correspScaled[ corr ].m_point3D[ 0 ];
	  mat_A( 2*corr, 5 ) = - mv_correspScaled[ corr ].m_point3D[ 1 ];
	  mat_A( 2*corr, 6 ) = - mv_correspScaled[ corr ].m_point3D[ 2 ];
	  mat_A( 2*corr, 7 ) = - 1.0;

	  mat_A( 2*corr, 8 ) = mv_correspScaled[ corr ].m_point3D[ 0 ]
		  * mv_correspScaled[ corr ].m_point2D[ 1 ];
	  mat_A( 2*corr, 9 ) = mv_correspScaled[ corr ].m_point3D[ 1 ]
		  * mv_correspScaled[ corr ].m_point2D[ 1 ];
	  mat_A( 2*corr, 10 ) = mv_correspScaled[ corr ].m_point3D[ 2 ]
		  * mv_correspScaled[ corr ].m_point2D[ 1 ];
	  mat_A( 2*corr, 11 ) = 1.0 * mv_correspScaled[ corr ].m_point2D[ 1 ];


	  mat_A( 2*corr + 1, 0  ) = mv_correspScaled[ corr ].m_point3D[ 0 ];
	  mat_A( 2*corr + 1, 1  ) = mv_correspScaled[ corr ].m_point3D[ 1 ];
	  mat_A( 2*corr + 1, 2  ) = mv_correspScaled[ corr ].m_point3D[ 2 ];
	  mat_A( 2*corr + 1, 3  ) = 1.0;
	  
	  mat_A( 2*corr + 1, 8  ) = - mv_correspScaled[ corr ].m_point3D[ 0 ]
		  * mv_correspScaled[ corr ].m_point2D[ 0 ];
	  mat_A( 2*corr + 1, 9  ) = - mv_correspScaled[ corr ].m_point3D[ 1 ]
		  * mv_correspScaled[ corr ].m_point2D[ 0 ];
	  mat_A( 2*corr + 1, 10  ) = - mv_correspScaled[ corr ].m_point3D[ 2 ]
		  * mv_correspScaled[ corr ].m_point2D[ 0 ];
	  mat_A( 2*corr + 1, 11  ) = - 1.0 * mv_correspScaled[ corr ].m_point2D[ 0 ];
    }




    ////////////////////////////////////////////////////
    // solve system A * p = 0

    if( ! SVD_V( mat_A, mat_VT ) )
    {
	return false;
    }

    // original code from Martin, replaced by simply using the eigenvector belonging
    // to the smallest eigenvalue.
    ////////////////////////////////////////////////////
    // find matrix that generates least projection error
    
    ProjMatrix currMatrix;
    Vector2D projPoint;

    double bestError = DBL_MAX, errorSum, error;
    int row, col;

    for( int matrix = 11; matrix >= 5; --matrix )
    {
	  for( row = 0; row < 3; ++row )
	  {
		  for( col = 0; col < 4; ++col )
		  {
			currMatrix( row, col ) = mat_VT( matrix, 4*row + col  );
		  }
	  }

	  errorSum = 0.0;
	  
	  for( corr = 0; corr < endCorr; ++corr )
	  {
		  currMatrix.project( mv_correspScaled[ corr ].m_point3D, projPoint );

		  error = (projPoint - mv_correspScaled[ corr ].m_point2D).length();
		  errorSum += error;
	  }
	  
      if( errorSum < bestError )
	  {
		  bestError = errorSum;
		  m_projectionMatrix = currMatrix;
	  }
    }
    

    return true;  
}



bool SolverProj::computeLinear( void )
{
    if( !computeLinearNoRescaling( ) )
      return false;

    ////////////////////////////////////////////////////////
    // undo coordinate scaling

    m_projectionMatrix.multLeft( m_matScaleImgInv );
    m_projectionMatrix.multRight( m_matScaleWorld );

    return true;
}





bool SolverProj::computeLinearNew( void )
{
    int nrows, ncols;

    nrows = 3 * mv_corresp.size();
    ncols = 12;

    Matrix mat_A( nrows, ncols );
    Matrix mat_VT( ncols, ncols );
    Vector vec_point( 4 );

    gmm::clear( mat_A );


    ///////////////////////////////////////////////////
    // coordinate scaling

    if( !scaleCorrespondences() )
      return false;



    int corr, endCorr = mv_correspScaled.size();

    ///////////////////////////////////////////////////
    // generate matrix A

    for( corr = 0; corr < endCorr; ++corr )
    {
	  vec_point[ 0 ] = mv_correspScaled[ corr ].m_point3D[ 0 ];
	  vec_point[ 1 ] = mv_correspScaled[ corr ].m_point3D[ 1 ];
	  vec_point[ 2 ] = mv_correspScaled[ corr ].m_point3D[ 2 ];
	  vec_point[ 3 ] = 1.0;

	  mat_A( 3*corr, 0 ) = - vec_point[ 0 ];
	  mat_A( 3*corr, 1 ) = - vec_point[ 1 ];
	  mat_A( 3*corr, 2 ) = - vec_point[ 2 ];
	  mat_A( 3*corr, 3 ) = - vec_point[ 3 ];

	  mat_A( 3*corr, 8 ) = vec_point[ 0 ]
		  * mv_correspScaled[ corr ].m_point2D[ 0 ];
	  mat_A( 3*corr, 9 ) = vec_point[ 1 ]
		  * mv_correspScaled[ corr ].m_point2D[ 0 ];
	  mat_A( 3*corr, 10 ) = vec_point[ 2 ]
		  * mv_correspScaled[ corr ].m_point2D[ 0 ];
	  mat_A( 3*corr, 11 ) = vec_point[ 3 ]
		  * mv_correspScaled[ corr ].m_point2D[ 0 ];


	  mat_A( 3*corr + 1, 4  ) = - vec_point[ 0 ];
	  mat_A( 3*corr + 1, 5  ) = - vec_point[ 1 ];
	  mat_A( 3*corr + 1, 6  ) = - vec_point[ 2 ];
	  mat_A( 3*corr + 1, 7  ) = - vec_point[ 3 ];

	  mat_A( 3*corr + 1, 8  ) = vec_point[ 0 ]
		  * mv_correspScaled[ corr ].m_point2D[ 1 ];
	  mat_A( 3*corr + 1, 9  ) = vec_point[ 1 ]
		  * mv_correspScaled[ corr ].m_point2D[ 1 ];
	  mat_A( 3*corr + 1, 10  ) = vec_point[ 2 ]
		  * mv_correspScaled[ corr ].m_point2D[ 1 ];
	  mat_A( 3*corr + 1, 11  ) = vec_point[ 3 ]
		  * mv_correspScaled[ corr ].m_point2D[ 1 ];



	  mat_A( 3*corr + 2, 0  ) = -vec_point[ 0 ]
		  * mv_correspScaled[ corr ].m_point2D[ 1 ];
	  mat_A( 3*corr + 2, 1  ) = -vec_point[ 1 ]
		  * mv_correspScaled[ corr ].m_point2D[ 1 ];
	  mat_A( 3*corr + 2, 2  ) = -vec_point[ 2 ]
		  * mv_correspScaled[ corr ].m_point2D[ 1 ];
	  mat_A( 3*corr + 2, 3  ) = -vec_point[ 3 ]
		  * mv_correspScaled[ corr ].m_point2D[ 1 ];

	  mat_A( 3*corr + 2, 4  ) = vec_point[ 0 ]
		  * mv_correspScaled[ corr ].m_point2D[ 0 ];
	  mat_A( 3*corr + 2, 5  ) = vec_point[ 1 ]
		  * mv_correspScaled[ corr ].m_point2D[ 0 ];
	  mat_A( 3*corr + 2, 6  ) = vec_point[ 2 ]
		  * mv_correspScaled[ corr ].m_point2D[ 0 ];
	  mat_A( 3*corr + 2, 7  ) = vec_point[ 3 ]
		  * mv_correspScaled[ corr ].m_point2D[ 0 ];
    }


    ////////////////////////////////////////////////////
    // solve system A * p = 0

    if( ! SVD_V( mat_A, mat_VT ) )
    {
	  return false;
    }


    // last column of V is solution ( = last row of VT)
    int col, row;
    for( row = 0; row < 3; ++row )
    {
	  for( col = 0; col < 4; ++col )
	  {
		  m_projectionMatrix( row, col ) = mat_VT( 11, 4*row + col  );
	  }
    }


    ////////////////////////////////////////////////////////
    // undo coordinate scaling

    m_projectionMatrix.multLeft( m_matScaleImgInv );
    m_projectionMatrix.multRight( m_matScaleWorld );

    return true;
}



double SolverProj::evaluateCorrespondence( const Vector2D & _point2D,
					   const Vector3D & _point3D )
{
    static double retval, val;
    Vector2D imgPoint;

    m_projectionMatrix.project( _point3D, imgPoint );
    
    val = _point2D[ 0 ] - imgPoint[ 0 ];
    retval = val * val;
    val = _point2D[ 1 ] - imgPoint[ 1 ];
    retval += val * val;

    return retval;
}


bool SolverProj::getPositionAndOrientation( Vector3D &position, Vector3D &orientation )
{
  // get the position of the camera as the vector spanning the right null-space of the 
  // projection matrix, see Hartley & Zisserman, 2nd edition, pages 158-159
  Matrix mat_VT( 4, 4 ), mat_P( 3, 4 );
  for( int i=0; i<3; ++i )
  {
    for( int j=0; j<4; ++j )
      mat_P(i,j) = m_projectionMatrix(i,j);
  }
  
  if( ! SVD_V( mat_P, mat_VT ) )
  {
      return false;
  }
  
  for( int i=0; i<3; ++i )
    position[i] = mat_VT( 3,i ) / mat_VT( 3,3 );
  
  // get the viewing direction of the camera, see Hartley & Zisserman, 2nd ed., pages 160 - 161
  // compute the determinant of the 3x3 part of the projection matrix
  double det = m_projectionMatrix(0,0) * m_projectionMatrix(1,1) * m_projectionMatrix(2,2) + m_projectionMatrix(0,1) * m_projectionMatrix(1,2) * m_projectionMatrix(2,0) + m_projectionMatrix(0,2) * m_projectionMatrix(1,0) * m_projectionMatrix(2,1) - m_projectionMatrix(0,2) * m_projectionMatrix(1,1) * m_projectionMatrix(2,0) - m_projectionMatrix(0,1) * m_projectionMatrix(1,0) * m_projectionMatrix(2,2) - m_projectionMatrix(0,0) * m_projectionMatrix(1,2) * m_projectionMatrix(2,1);
  
  // remember that the camera in reconstructions computed by Bundler looks 
  // down the negative z-axis instead of the positive z-axis.
  // So we have to multiply the orientation with -1.0
  for( int i=0; i<3; ++ i )
    orientation[i] = - m_projectionMatrix(2,i) * det;
  
  return true;
}

}

}
