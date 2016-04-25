cmake_minimum_required( VERSION 2.8 FATAL_ERROR )
project( registration )

#########################################################
# FIND CGAL
#########################################################
find_package(CGAL REQUIRED)
include(${CGAL_USE_FILE})
#include_directories(${CGAL_INCLUDE_DIRS})
#link_directories(${CGAL_LIBRARY_DIRS})
#add_definitions(${CGAL_DEFINITIONS})
message( ${CGAL_LIBRARY} )
if(NOT CGAL_FOUND)
    message(ERROR " CGAL not found!")
endif(NOT CGAL_FOUND)

#########################################################
# FIND EIGEN3
#########################################################
find_path( EIGEN Dense )
include_directories( ${EIGEN} )

#########################################################
# FIND LBFGS
#########################################################
find_path(LBFGS_INCLUDE_DIR lbfgs.h )
include_directories( ${LBFGS_INCLUDE_DIR} )
find_library(LBFGS_LIBRARY NAMES liblbfgs.so )

file(GLOB SOURCES "*.cpp")
message( ${SOURCES} )
add_executable( demon ${SOURCES} )
target_link_libraries( demon ${LBFGS_LIBRARY} )
target_link_libraries( demon ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES})
target_link_libraries( demon ${BOOST_LIBRARY} )
