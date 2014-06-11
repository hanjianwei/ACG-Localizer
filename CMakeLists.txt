cmake_minimum_required (VERSION 2.6)

project (acg_localizer)

# add our macro directory to cmake search path
if (EXISTS "${CMAKE_SOURCE_DIR}/cmake")
  set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
endif ()

if ( NOT ENABLE_TESTS )
    set (ENABLE_TESTS OFF CACHE BOOL "Turn on testing framework.")
endif()

if ( NOT STL_VECTOR_CHECKS )
    set ( STL_VECTOR_CHECKS false CACHE BOOL "Include full stl vector checks in debug mode (This option is only used in debug Mode!)" )
endif ( NOT STL_VECTOR_CHECKS )

if ( STL_VECTOR_CHECKS )
    add_definitions (
        -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC
    )
endif ()


# use sse4_2
if (CMAKE_C_COMPILER MATCHES ".*gcc.*$")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.2")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse4.2")
endif()


if (CMAKE_C_COMPILER MATCHES ".*icc.*$")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ip -axSSE4.2,SSSE3,SSE2 -fp-model fast=2")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -ip -axSSE4.2,SSSE3,SSE2 -fp-model fast=2")
endif()



set( CMAKE_DEBUG_POSTFIX "d" )


add_subdirectory (src)

# testing
if( ENABLE_TESTS )
    add_subdirectory (test)
endif()

