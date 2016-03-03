cmake_minimum_required (VERSION 2.6)
project (libcubs C)

include_directories( "${CMAKE_SOURCE_DIR}/include" )
set( CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -O3" )

add_library( libcubs
	"${CMAKE_SOURCE_DIR}/src/spline.c"
	"${CMAKE_SOURCE_DIR}/include/cubs/spline.h"
)

add_executable( spline_test
	"${CMAKE_SOURCE_DIR}/src/splinetest.c"
)

target_link_libraries ( spline_test libcubs m )
