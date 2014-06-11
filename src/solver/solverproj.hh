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
 * Declaration of solver classes based on point correspondences
 *
 * Authors: Martin Habbecke <M.Habbecke@gmx.de>, Torsten Sattler <tsattler@cs.rwth-aachen.de>
 * Date: 10-03-2011
 */



#ifndef _FT_SOLVER_PROJ_HH_
#define _FT_SOLVER_PROJ_HH_


#include <vector>

#include "../math/math.hh"
#include "../math/matrix3x3.hh"
#include "../math/projmatrix.hh"
#include "solverbase.hh"


using namespace Util::Math;

namespace Util {
namespace CorrSolver {



//! Solver that computes 
class SolverProj {

public:

    //! Default constructor
    SolverProj();

    //! Destructor
    virtual ~SolverProj();


    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////

    //! Correspondence class for 2-image correspondences
    class Correspondence {
	
    public:
	
	  //! Initializing constructor
	  Correspondence( const Vector2D & _p2d, const Vector3D & _p3d )
		  : m_point2D( _p2d ), m_point3D( _p3d )
		  {}

	  inline const Vector2D sortPoint( void ) const
		  { return m_point2D; }

	  //! Feature in first image
	  Vector2D m_point2D;
	  
	  //! Feature in second image
	  Vector3D m_point3D;
    };

    //! Type for a vector of correspondences
    typedef std::vector< Correspondence > CorrespondenceVec;


    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////


    //! Clear all point correspondences
    void clear( void );

    //! Add a point correspondence
    void addCorrespondence( const Vector2D & _point2D,
			    const Vector3D & _point3D );

    //! Add a correspondence
    void addCorrespondence( const Correspondence & _corres );


    //! Compute projection matrix using SVD
    bool computeLinear( void );
    
    //! Compute projection matrix with extended linear matrix
    bool computeLinearNew( void );

    //! Get computed projection matrix
    void getProjectionMatrix( Matrix & _mat )
	{ m_projectionMatrix.copyTo( _mat ); }

    //! Get computed projection matrix
    void getProjectionMatrix( ProjMatrix & _mat )
	{ _mat = m_projectionMatrix; }
	
    //! Set the projection matrix
    void setProjectionMatrix( ProjMatrix & _mat )
    	{ m_projectionMatrix = _mat; }


    //! Set initial projection matrix for iterative solver
    void setInitialMatrix( const Matrix & _mat );

    //! Set initial projection matrix for iterative solver
    void setInitialMatrix( const ProjMatrix & _mat );


    //! Evaluate a point correspondence (compute squared reprojection error)
    virtual double evaluateCorrespondence( const Vector2D & _point2D,
					   const Vector3D & _point3D );
					   
    //! get the position and orientation of the camera as two 3D vectors
    bool getPositionAndOrientation( Vector3D &position, Vector3D &orientation );


    //! For testing only
    int m_totalNumIterations;

protected:
  
    //! Compute projection matrix using SVD without rescaling the resulting projection matrix
    bool computeLinearNoRescaling( void );

    //! Function that scales original features, returns false if one of the scaling factors is too close to 0
    bool scaleCorrespondences( void );


    //! Vector of all point correspondences
    CorrespondenceVec mv_corresp;

    //! Vector of all point correspondences, scaled such that x \in [-1.0, 1.0]
    CorrespondenceVec mv_correspScaled;



    //! Computed projection matrix
    ProjMatrix m_projectionMatrix;

    //! Initial projection matrix
    ProjMatrix m_initialMatrix;

    //! Scaling matrix for image space
    Matrix3x3 m_matScaleImgInv;
    
    //! Scaling matrix for world space
    Matrix4x4 m_matScaleWorld;
};


}

}


#endif // _FT_SOLVER_PROJ_HH_
