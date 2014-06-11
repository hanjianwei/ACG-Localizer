# - Try to find OPENMESH_VEC (the simple version including only the vector class)
# Once done this will define
#  
# OPENMESH_VEC_FOUND           - system has OPENMESH_vec
# OPENMESH_VEC_INCLUDE_DIR - theOPENMESH include directory

# IF (OPENMESH_INCLUDE_DIR)
 # Already in cache, be silent
# SET(OPENMESH_FIND_QUIETLY TRUE)
# ENDIF (OPENMESH_INCLUDE_DIR)


FIND_PATH(OPENMESH_INCLUDE_DIR OpenMesh/Core/Geometry/VectorT.hh 
	PATHS "${CMAKE_SOURCE_DIR}/src/OpenMesh_vec" )
	
if( OPENMESH_INCLUDE_DIR )
    set( OPENMESH_FOUND TRUE )
else( OPENMESH_INCLUDE_DIR )
    set( OPENMESH_FOUND FALSE )
endif( OPENMESH_INCLUDE_DIR )