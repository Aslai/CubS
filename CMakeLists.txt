cmake_minimum_required (VERSION 2.6)
project (libcubs C)

include_directories( "${CMAKE_SOURCE_DIR}/include" )

add_library( libcubs 
	"${CMAKE_SOURCE_DIR}/src/spline.c" 
)
