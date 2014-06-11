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
 *	       Torsten Sattler
 *
 * Version: Revision: 75
 * Date:    2011-10-02 18:00:00
 */



#include "math.hh"



namespace Util {
namespace Math {


extern "C" {

    typedef void (*U_fp)(...);

    int dgesvd_( char *jobu, char *jobvt, long int *m, long int *n,
		 double *a, long int *lda, double *s, double *u, long int *ldu,
		 double *vt, long int *ldvt, double *work, long int *lwork,
		 long int *info );

    int sgesvd_( char *jobu, char *jobvt, long int *m, long int *n,
		 float *a, long int *lda, float *s, float *u, long int *ldu,
		 float *vt, long int *ldvt, float *work, long int *lwork,
		 long int *info );

    int dgesdd_( char *jobz, long int *m, long int *n, double *a,
		 long int *lda, double *s, double *u, long int *ldu, 
		 double *vt, long int *ldvt, double *work, long int *lwork, 
		 long int *iwork, long int *info);


    int dggev_( char *jobvl, char *jobvr, long int *n, double *a, 
		long int *lda, double *b, long int *ldb, double *alphar,
		double *alphai, double *beta, double *vl, long int *ldvl,
		double *vr, long int *ldvr, double *work, long int *lwork,
		long int *info );

    int dsygv_( long int *itype, char *jobz, char *uplo, long int *n, 
		double *a, long int *lda, double *b, long int *ldb, 
		double *w, double *work, long int *lwork, long int *info);

    int dgeev_( char *jobvl, char *jobvr, long int *n, double *a,
		long int *lda, double *wr, double *wi, double *vl,
		long int *ldvl, double *vr, long int *ldvr, double *work,
		long int *lwork, long int *info );
    
    int sgeev_( char *jobvl, char *jobvr, long int *n, float *a,
		long int *lda, float *wr, float *wi, float *vl,
  		long int *ldvl, float *vr, long int *ldvr, float *work,
 		long int *lwork, long int *info );

    int dgesv_( long int *n, long int *nrhs, double *a, long int *lda, 
		long int *ipiv, double *b, long int *ldb, long int *info );

    int sgesv_( long int *n, long int *nrhs, float *a, long int *lda, 
		long int *ipiv, float *b, long int *ldb, long int *info );

    int dsyev_( char *jobz, char *uplo, long int *n, double *a,
		long int *lda, double *w, double *work, long int *lwork, 
		long int *info);

    int dgeqrf_( long int *m, long int *n, double * a, long int * lda,
		 double * tau, double * work, long int *lwork, long int *info );

    int dposv_( char *uplo, long int *n, long int *nrhs, double *a,
		long int *lda, double *b, long int *ldb, long int *info );

    int dsysv_( char *uplo, long int *n, long int *nrhs, double *a,
		long int *lda, long int *ipiv, double *b, long int *ldb, 
		double *work, long int *lwork, long int *info );

    int dgels_( char *trans, long int *m, long int *n, long int *nrhs,
		double *a, long int* lda, 
		double *b, long int* ldb, double *work, long int *lwork, 
		long int *info);

    int sgels_( char *trans, long int *m, long int *n, long int *nrhs,
		float *a, long int* lda, 
		float *b, long int* ldb, float *work, long int *lwork, 
		long int *info);

    int dgetri_( long int *n, double *a, long int *lda, long int *ipiv,
		 double *work, long int *lwork, long int *info );

    int dgetrf_( long int *m, long int *n, double *a, long int * lda,
		 long int *ipiv, long int *info );

    int dpotrf_( char *uplo, long int *n, double *a, long int *lda, long int *info );
    
    int dgeqp3_( long int *m, long int *n, double *a, long int *lda, int *jvpt, double *tau, double *work, long int *lwork, long int*info );

};



//! Module-wide pointer to working array
double * p_work_double = 0;
float * p_work_float = 0;
long int * p_work_longint = 0;

//! Module-wide size of working array
int m_workSize_double = 0;
int m_workSize_float = 0;
int m_workSize_longint = 0;




void FreeWorkingArray( void )
{
    if( p_work_double )
	delete [] p_work_double;
    if( p_work_float )
	delete [] p_work_float;
    if( p_work_longint )
	delete [] p_work_longint;

    p_work_double = 0;
    m_workSize_double = 0;
    p_work_float = 0;
    m_workSize_float = 0;
    p_work_longint = 0;
    m_workSize_longint = 0;
}






///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////





bool SVD( Matrix & _mat_A, Matrix & _mat_U, 
	  Matrix & _mat_VT, Vector & _vec_D )
{
    long int info = 0, lda, ldu, ldvt, lwork, rows, cols;
    char jobu, jobvt;
    double size;


    jobu = 'A';
    jobvt = 'A';
    rows = gmm::mat_nrows( _mat_A );
    cols = gmm::mat_ncols( _mat_A );
    lda = gmm::mat_nrows( _mat_A );
    ldu = rows;
    ldvt = cols;

    // compute optimal size for work-array first
    lwork = -1;
    dgesvd_( &jobu, &jobvt, &rows, &cols, &_mat_A( 0, 0 ), 
	     &lda, &_vec_D[ 0 ], &_mat_U( 0, 0 ), &ldu, &_mat_VT( 0, 0 ), 
	     &ldvt, &size, &lwork, &info );
    
    if( info != 0 )
	return false;

    // re-allocate working array if necessary
    if( (long int) size > m_workSize_double )
    {
	if( p_work_double )
	    delete [] p_work_double;
	
	m_workSize_double = (long int) size;
	p_work_double = new double[ m_workSize_double ];
    }
    lwork = m_workSize_double;

    // compute SVD
    dgesvd_( &jobu, &jobvt, &rows, &cols, &_mat_A( 0, 0 ), 
	     &lda, &_vec_D[ 0 ], &_mat_U( 0, 0 ), &ldu, &_mat_VT( 0, 0 ), 
	     &ldvt, p_work_double, &lwork, &info );
    
    if( info != 0 )
	return false;
    
    return true;
}



bool SVD_V( Matrix & _mat_A, Matrix & _mat_VT )
{
    Vector vec_D( std::min( gmm::mat_ncols( _mat_A ), gmm::mat_nrows( _mat_A ) ) );

    long int info = 0, lda, ldu, ldvt, lwork, rows, cols;
    char jobu, jobvt;
    double size;


    jobu = 'N';
    jobvt = 'A';
    rows = gmm::mat_nrows( _mat_A );
    cols = gmm::mat_ncols( _mat_A );
    lda = gmm::mat_nrows( _mat_A );
    ldu = 1;
    ldvt = cols;

    // compute optimal size for work-array first
    lwork = -1;
    dgesvd_( &jobu, &jobvt, &rows, &cols, &_mat_A( 0, 0 ), 
	     &lda, &vec_D[ 0 ], 0, &ldu, &_mat_VT( 0, 0 ), 
	     &ldvt, &size, &lwork, &info );
    
    if( info != 0 )
    {
	std::cout << "info (1): " << info << std::endl;
	return false;
    }

    // re-allocate working array if necessary
    if( (long int) size > m_workSize_double )
    {
	  if( p_work_double )
		  delete [] p_work_double;
	  m_workSize_double = (long int) size;
	  p_work_double = new double[ m_workSize_double ];
    }
    lwork = m_workSize_double;

	// compute SVD
    dgesvd_( &jobu, &jobvt, &rows, &cols, &_mat_A( 0, 0 ), 
	     &lda, &vec_D[ 0 ], 0, &ldu, &_mat_VT( 0, 0 ), 
	     &ldvt, p_work_double, &lwork, &info );

		 
    if( info != 0 )
    {
	  std::cout << "info (2): " << info << std::endl;
	  return false;
    }
    
    return true;
}




int dgetrf_( long int *m, long int *n, double *a, long int * lda,
		 long int *ipiv, long int *info );






} 
}
