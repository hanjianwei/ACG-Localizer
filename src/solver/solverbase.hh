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
 * Author: Martin Habbecke <M.Habbecke@gmx.de>
 * Version: Revision: 1.10
 * Date: 2011-10-03 17:30
 */



#ifndef _FT_SOLVER_BASE_HH_
#define _FT_SOLVER_BASE_HH_


#include <vector>

#include <OpenMesh/Core/Geometry/VectorT.hh>

#include "../math/math.hh"
#include "../math/matrix3x3.hh"


using namespace Util::Math;

namespace Util {
namespace CorrSolver {


//! 2D point type
typedef OpenMesh::Vec2d Vector2D;

//! 3D point type
typedef OpenMesh::Vec3d Vector3D;

//! 4D point type
typedef OpenMesh::Vec4d Vector4D;





///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////



//! Base class for solvers based on 2-image correspondences
class SolverBase {

public:
    //! Default constructor
    SolverBase();

    //! Destructor
    virtual ~SolverBase();


    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////

    //! Correspondence class for 2-image correspondences
    class Correspondence {
	
    public:
	
	  //! Default constructor
	  Correspondence() {}

	  //! Initializing constructor
	  Correspondence( const Vector2D & _f1, const Vector2D & _f2 )
		  : m_feature1( _f1 ), m_feature2( _f2 )
		  {}

	  inline const Vector2D sortPoint( void ) const
		  { return m_feature1; }

	  //! Feature in first image
	  Vector2D m_feature1;
	  
	  //! Feature in second image
	  Vector2D m_feature2;
    };

    //! Type for a vector of correspondences
    typedef std::vector< Correspondence > CorrespondenceVec;

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////


    //! Clear all point correspondences
    void clear( void );

    //! Add a correspondence
    void addCorrespondence( const Vector2D & _point1,
			    const Vector2D & _point2 )
	{ mv_corresp.push_back( Correspondence( _point1, _point2 ) ); }

    //! Add a correspondence
    void addCorrespondence( const Correspondence & _corresp )
	{ mv_corresp.push_back( _corresp ); }

    //! Set vector of correspondences
    void setCorrespondences( const CorrespondenceVec & _corrVec )
	{ mv_corresp = _corrVec; }



protected:

    //! Function that scales original features
    void scaleCorrespondences( void );



    //! Vector of all point correspondences
    CorrespondenceVec mv_corresp;

    //! Vector of all point correspondences, scaled such that x \in [-1.0, 1.0]
    CorrespondenceVec mv_correspScaled;


    //! Matrix to scale coords of first image (orig -> scaled)
    Matrix3x3 m_matScale1;

    //! Transposed matrix to scale coords of second image
    Matrix3x3 m_matScale2T;


    //! C.o.g. of points in first image
    Vector2D m_center1;

    //! C.o.g. of points in second image
    Vector2D m_center2;


    //! Scale factor for points in first image
    double m_scale1;

    //! Scale factor for points in second image
    double m_scale2;


};


}

}


#endif // _FT_SOLVER_BASE_HH_
