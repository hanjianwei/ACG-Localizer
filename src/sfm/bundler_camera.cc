/*===========================================================================*\
 *                                                                           *
 *                            ACG Localizer                                  *
 *      Copyright (C) 2011-2012 by Computer Graphics Group, RWTH Aachen      *
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

#include "bundler_camera.hh"

bundler_camera::bundler_camera()
{
  rotation.setIdentity();
  focal_length = 1.0;
  kappa_1 = kappa_2 = 0.0;
  translation[0] = translation[1] = translation[2] = 0.0;
  width = height = 0;
}

bundler_camera::~bundler_camera()
{}

//-----------------------------------------------------------------------------

OpenMesh::Vec2f bundler_camera::project_f( const OpenMesh::Vec3d &p ) const
{
  OpenMesh::Vec3d tmp;
  rotation.multRight( p, tmp );
  tmp += translation;
  tmp[2] *= -1.0;
  for( int i=0; i<2; ++i )
    tmp[i] /= tmp[2];
  
  double length_sqr = pow( tmp[0], 2.0 ) + pow( tmp[1], 2.0 );
  // compute radial distortion term
  double r = 1.0 + kappa_1 * length_sqr + kappa_2 * length_sqr * length_sqr;
  return OpenMesh::Vec2f( focal_length * r * tmp[0], focal_length * r * tmp[1] );
}

//-----------------------------------------------------------------------------

OpenMesh::Vec2f bundler_camera::project_f_undist( const OpenMesh::Vec3d &p ) const
{
  OpenMesh::Vec3d tmp;
  rotation.multRight( p, tmp );
  tmp += translation;
  tmp[2] *= -1.0;
  for( int i=0; i<2; ++i )
    tmp[i] /= tmp[2];
  
  return OpenMesh::Vec2f( focal_length * tmp[0], focal_length * tmp[1] );
}

//-----------------------------------------------------------------------------

OpenMesh::Vec2d bundler_camera::project_d( const OpenMesh::Vec3d &p ) const
{
  OpenMesh::Vec3d tmp;
  rotation.multRight( p, tmp );
  tmp += translation;
  tmp[2] *= -1.0;
  for( int i=0; i<2; ++i )
    tmp[i] /= tmp[2];
  
  double length_sqr = pow( tmp[0], 2.0 ) + pow( tmp[1], 2.0 );
  // compute radial distortion term
  double r = 1.0 + kappa_1 * length_sqr + kappa_2 * length_sqr * length_sqr;
  return OpenMesh::Vec2d( focal_length * r * tmp[0], focal_length * r * tmp[1] );
}


//-----------------------------------------------------------------------------

OpenMesh::Vec2d bundler_camera::project_d_undist( const OpenMesh::Vec3d &p ) const
{
  OpenMesh::Vec3d tmp;
  rotation.multRight( p, tmp );
  tmp += translation;
  tmp[2] *= -1.0;
  for( int i=0; i<2; ++i )
    tmp[i] /= tmp[2];
  
  return OpenMesh::Vec2d( focal_length * tmp[0], focal_length * tmp[1] );
}

//-----------------------------------------------------------------------------



OpenMesh::Vec3f bundler_camera::get_cam_global_vec_f(const OpenMesh::Vec3d &_local) const
{
    OpenMesh::Vec3d g_pos;
    rotation.multLeft( _local, g_pos );
    return OpenMesh::Vec3f( g_pos[0], g_pos[1], g_pos[2] );
}

//-----------------------------------------------------------------------------

OpenMesh::Vec3d bundler_camera::get_cam_global_vec_d(const OpenMesh::Vec3d &_local) const
{
    OpenMesh::Vec3d g_pos;
    rotation.multLeft( _local, g_pos );
    return g_pos;
}

//-----------------------------------------------------------------------------

OpenMesh::Vec3f bundler_camera::get_cam_position_f() const
{
  OpenMesh::Vec3d g_pos;
  rotation.multLeft( translation, g_pos );
  return OpenMesh::Vec3f( -g_pos[0], -g_pos[1], -g_pos[2] );
}

//-----------------------------------------------------------------------------

OpenMesh::Vec3d bundler_camera::get_cam_position_d() const
{
  OpenMesh::Vec3d g_pos;
  rotation.multLeft( translation, g_pos );
  return OpenMesh::Vec3d( -g_pos[0], -g_pos[1], -g_pos[2] );
}

//-----------------------------------------------------------------------------

double bundler_camera::compute_reprojection_error_f( OpenMesh::Vec3f &p, float x, float y ) const
{
  double error = 0.0;
  
  OpenMesh::Vec3d tmp, tmp2( p[0], p[1], p[2] );
  rotation.multRight( tmp2, tmp );
  tmp += translation;
  tmp[2] *= -1.0;
  for( int i=0; i<2; ++i )
    tmp[i] /= tmp[2];
  
  double length_sqr = pow( tmp[0], 2.0 ) + pow( tmp[1], 2.0 );
  // compute radial distortion term
  double r = 1.0 + kappa_1 * length_sqr + kappa_2 * length_sqr * length_sqr;
  
  OpenMesh::Vec2d proj_point( focal_length * r * tmp[0], focal_length * r * tmp[1] );
  
  error += ( x - proj_point[0] ) * ( x - proj_point[0] );
  error += ( y - proj_point[1] ) * ( y - proj_point[1] );
  
//   std::cout << proj_point << " vs ( " << x << " , " << y << " ) " << std::endl;
  
  return sqrt(error);
}

//-----------------------------------------------------------------------------

OpenMesh::Vec3d bundler_camera::to_cam_coords_d( const OpenMesh::Vec3d &p ) const
{
  OpenMesh::Vec3d tmp;
  rotation.multRight( p, tmp );
  tmp += translation;
  
  return tmp;
}

