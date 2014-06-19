# - Find GMM
# Find the native GMM includes and library
# This module defines
#
#  GMM_INCLUDE_DIR, where to find GMM/GMM.h, etc.
#  GMM_FOUND, If false, do not try to use GMM.

include(FindPackageHandleStandardArgs)

find_path(GMM_INCLUDE_DIR gmm/gmm.h)

find_package_handle_standard_args(GMM DEFAULT_MSG
  GMM_INCLUDE_DIR)

if(GMM_FOUND)
  set(GMM_INCLUDE_DIRS ${GMM_INCLUDE_DIR})
endif()
