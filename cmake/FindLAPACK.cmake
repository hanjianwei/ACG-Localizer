# - Try to find LAPACK
# Once done this will define
#  
# LAPACK_FOUND           - system has LAPACK
# LAPACK_INCLUDE_DIR - theLAPACK include directory
# LAPACK_LIBRARY         - Link these to use LAPACK
# LAPACK_LIBRARY_DIR  - Library DIR of LAPACK
#   

IF (LAPACK_LIBRARIES)
  # Already in cache, be silent
  SET(LAPACK_FIND_QUIETLY TRUE)
ENDIF (LAPACK_LIBRARIES)


if( WIN32 )
  
  FIND_LIBRARY(LAPACK_LIB_opt
    NAMES lapack 
    PATHS "${CMAKE_SOURCE_DIR}/../win32libs")
  
  FIND_LIBRARY(BLAS_LIB_opt 
    NAMES blas 
    PATHS "${CMAKE_SOURCE_DIR}/../win32libs")
  
  
  FIND_LIBRARY(F2C_LIB_opt
    NAMES libf2c 
    PATHS "${CMAKE_SOURCE_DIR}/../win32libs")
  
  GET_FILENAME_COMPONENT( LAPACK_LIB_opt ${LAPACK_LIB_opt} NAME )
  GET_FILENAME_COMPONENT( BLAS_LIB_opt ${BLAS_LIB_opt} NAME )
  GET_FILENAME_COMPONENT( F2C_LIB_opt ${F2C_LIB_opt} NAME )
  
  
  
  FIND_LIBRARY(LAPACK_LIB_debug
    NAMES lapackd 
    PATHS "${CMAKE_SOURCE_DIR}/../win32libs")
  
  FIND_LIBRARY(BLAS_LIB_debug 
    NAMES blasd 
    PATHS "${CMAKE_SOURCE_DIR}/../win32libs")
  
  
  FIND_LIBRARY(F2C_LIB_debug
    NAMES libf2cd 
    PATHS "${CMAKE_SOURCE_DIR}/../win32libs")
  
  GET_FILENAME_COMPONENT( LAPACK_LIB_debug ${LAPACK_LIB_debug} NAME )
  GET_FILENAME_COMPONENT( BLAS_LIB_debug ${BLAS_LIB_debug} NAME )
  GET_FILENAME_COMPONENT( F2C_LIB_debug ${F2C_LIB_debug} NAME )
  
  
  if( LAPACK_LIB_opt AND BLAS_LIB_opt AND F2C_LIB_opt
      AND LAPACK_LIB_debug AND BLAS_LIB_debug AND F2C_LIB_debug )

    SET( LAPACK_FOUND TRUE )
    SET( LAPACK_LIBRARIES
      optimized ${LAPACK_LIB_opt} debug ${LAPACK_LIB_debug}
      optimized ${BLAS_LIB_opt} debug ${BLAS_LIB_debug}
      optimized ${F2C_LIB_opt} debug ${F2C_LIB_debug} )
    
  else()
    
    SET(LAPACK_FOUND FALSE )
    MESSAGE(STATUS "Lapack not found")
    
  endif()

else()


  FIND_LIBRARY(LAPACK_LIB
    NAMES lapack 
    PATHS /usr/lib /usr/local/lib)
  
  FIND_LIBRARY(BLAS_LIB 
    NAMES blas 
    PATHS /usr/lib /usr/local/lib)
  
  FIND_LIBRARY(F2C_LIB
    NAMES f2c 
    PATHS /usr/lib /usr/local/lib)

 GET_FILENAME_COMPONENT( LAPACK_LIBRARY_PATH ${LAPACK_LIB} PATH )

  if( LAPACK_LIB AND BLAS_LIB AND F2C_LIB )

    SET( LAPACK_FOUND TRUE )
      
    SET( LAPACK_LIBRARIES
      ${LAPACK_LIB}
      ${BLAS_LIB}
      ${F2C_LIB_opt} )
        
  else()
    
    SET(LAPACK_FOUND FALSE )
    MESSAGE(STATUS "Lapack not found")
    
  endif()

endif()




