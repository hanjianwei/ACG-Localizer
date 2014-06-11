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
 * Definition of Solver base class
 *
 * Author: Martin Habbecke <M.Habbecke@gmx.de>
 * Version: Revision: 1.8
 * Date: 2011-10-03 17:30
 */


#include <float.h>
#include <math.h>

#include "solverbase.hh"


using namespace Util::Math;

namespace Util {
namespace CorrSolver {



SolverBase::SolverBase()
{

}


SolverBase::~SolverBase()
{

}




void SolverBase::clear( void )
{
    mv_corresp.clear();
    mv_correspScaled.clear();
}




void SolverBase::scaleCorrespondences( void )
{
    CorrespondenceVec::iterator cIter, cEnd;
    cEnd = mv_corresp.end();

    Vector2D feat1, feat2;
    double  n;

    m_center1 = Vector2D( 0.0, 0.0 );
    m_center2 = Vector2D( 0.0, 0.0 );


    //////////////////////////////////////////////
    // determine c.o.g. in both images and also x- and y-ranges

    n = 1.0;

    for( cIter = mv_corresp.begin();
	 cIter != cEnd;
	 ++cIter )
    {
	  m_center1 = cIter->m_feature1 * (1.0/n) + m_center1 * ((n-1.0) / n);
	  m_center2 = cIter->m_feature2 * (1.0/n) + m_center2 * ((n-1.0) / n);

	  n += 1.0f;
    }

    ///////////////////////////////////////////////
    // move features to COG and compute scaling factors

    mv_correspScaled.clear();

    n = 1.0;
    m_scale1 = m_scale2 = 0.0;

    for( cIter = mv_corresp.begin();
	 cIter != cEnd;
	 ++cIter )
    {
	  feat1 = cIter->m_feature1 - m_center1;
	  feat2 = cIter->m_feature2 - m_center2;

	  m_scale1 = m_scale1 * ((n-1.0) / n) + feat1.norm() * (1.0 / n);
	  m_scale2 = m_scale2 * ((n-1.0) / n) + feat2.norm() * (1.0 / n);
	  mv_correspScaled.push_back( Correspondence( feat1, feat2 ) );

	  n += 1.0f;
    }
    m_scale1 = 1.41421 / m_scale1;
    m_scale2 = 1.41421 / m_scale2;

    //////////////////////////////////////////////
    // scale coords in both images

    cEnd = mv_correspScaled.end();
    for( cIter = mv_correspScaled.begin();
	 cIter != cEnd;
	 ++cIter )
    {
	  cIter->m_feature1 *= m_scale1;
	  cIter->m_feature2 *= m_scale2;
    }


    ////////////////////////////////////////////
    // create scaling matrices that will be multiplied
    // to result matrix

    m_matScale1.clear();
    m_matScale1( 0, 0 ) = m_scale1;
    m_matScale1( 1, 1 ) = m_scale1;
    m_matScale1( 0, 2 ) = -m_center1[ 0 ] * m_scale1;
    m_matScale1( 1, 2 ) = -m_center1[ 1 ] * m_scale1;
    m_matScale1( 2, 2 ) = 1.0;

    m_matScale2T.clear();
    m_matScale2T( 0, 0 ) = m_scale2;
    m_matScale2T( 1, 1 ) = m_scale2;
    m_matScale2T( 2, 0 ) = -m_center2[ 0 ] * m_scale2;
    m_matScale2T( 2, 1 ) = -m_center2[ 1 ] * m_scale2;
    m_matScale2T( 2, 2 ) = 1.0;
}



}

}
