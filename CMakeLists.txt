cmake_minimum_required( VERSION 2.8 FATAL_ERROR )
project( registration )

set(warnings "-Wall -Werror -Wextra")

#########################################################
# FIND GLUT
#########################################################
find_package(GLUT REQUIRED)
#include_directories(${GLUT_INCLUDE_DIRS})
#link_directories(${GLUT_LIBRARY_DIRS})
#add_definitions(${GLUT_DEFINITIONS})
if(NOT GLUT_FOUND)
    message(ERROR " GLUT not found!")
endif(NOT GLUT_FOUND)

#########################################################
# FIND OPENGL
#########################################################
find_package(OpenGL REQUIRED)
#include_directories(${OpenGL_INCLUDE_DIRS})
#link_directories(${OpenGL_LIBRARY_DIRS})
#add_definitions(${OpenGL_DEFINITIONS})
if(NOT OPENGL_FOUND)
    message(ERROR " OPENGL not found!")
endif(NOT OPENGL_FOUND)

#########################################################
# FIND CGAL
#########################################################
find_package(CGAL REQUIRED)
include(${CGAL_USE_FILE})
#include_directories(${CGAL_INCLUDE_DIRS})
#link_directories(${CGAL_LIBRARY_DIRS})
#add_definitions(${CGAL_DEFINITIONS})
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
#message( ${LBFGS_INCLUDE_DIR} )
include_directories( ${LBFGS_INCLUDE_DIR} )
find_library(LBFGS_LIBRARY NAMES liblbfgs.so )
#message( ${LBFGS_LIBRARY} )

add_executable( demon Demon.cpp )
target_link_libraries( demon ${LBFGS_LIBRARY} )
#target_link_libraries( sample ${LBFGS_LIBRARY} )
