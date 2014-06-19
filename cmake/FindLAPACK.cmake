# - Find LAPACK
# Find the native LAPACK includes and library
# This module defines
#
#  LAPACK_LIBRARIES, the libraries needed to use LAPACK.
#  LAPACK_FOUND, If false, do not try to use LAPACK.
# also defined, but not for general use are
#  LAPACK_LIBRARY, where to find the LAPACK library.

include(FindPackageHandleStandardArgs)

find_library(LAPACK_LIB lapack)
find_library(BLAS_LIB blas)
find_library(F2C_LIB f2c)

find_package_handle_standard_args(LAPACK DEFAULT_MSG
  LAPACK_LIB BLAS_LIB F2C_LIB)

if(LAPACK_FOUND)
  set(LAPACK_LIBRARIES ${LAPACK_LIB} ${BLAS_LIB} ${F2C_LIB})
  set(LAPACK_LIBRARY ${LAPACK_LIBRARIES})
endif()
