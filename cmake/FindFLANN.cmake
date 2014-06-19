# - Find FLANN
# Find the native FLANN includes and library
# This module defines
#
#  FLANN_INCLUDE_DIR, where to find ANN/ANN.h, etc.
#  FLANN_LIBRARIES, the libraries needed to use ANN.
#  FLANN_FOUND, If false, do not try to use ANN.
# also defined, but not for general use are
#  FLANN_LIBRARY, where to find the ANN library.

include(FindPackageHandleStandardArgs)

find_path(FLANN_INCLUDE_DIR flann/flann.h)
find_library(FLANN_LIBRARY flann_cpp
  PATHS /usr/local/lib64)

find_package_handle_standard_args(FLANN DEFAULT_MSG
  FLANN_LIBRARY FLANN_INCLUDE_DIR)

if(FLANN_FOUND)
  set(FLANN_INCLUDE_DIRS ${FLANN_INCLUDE_DIR})
  set(FLANN_LIBRARIES ${FLANN_LIBRARY})
endif()
