# Created by the script cgal_create_cmake_script_with_options
# This is the CMake script for compiling a set of CGAL applications.

project( Demons )


cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
    if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
        cmake_policy(VERSION 2.8.4)
    else()
        cmake_policy(VERSION 2.6)
    endif()
endif()

set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )

if ( COMMAND cmake_policy )

    cmake_policy( SET CMP0003 NEW )  

endif()

# CGAL and its components
find_package( CGAL )

if ( CGAL_FOUND )

    include( ${CGAL_USE_FILE} )
    include( ${EIGEN3_USE_FILE} )
#     # use either Eigen or BLAS/LAPACK
#     #find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)
#     if (NOT EIGEN3_FOUND)
#         message(STATUS "Eigen not found")
#         #find_package(LAPACK)
#         #if(LAPACK_FOUND)
#         #message("Lapack Found")
#         #  include( ${LAPACK_USE_FILE} )
#         #endif(LAPACK_FOUND)
#     else()
#         message("Eigen found")
#         include( ${EIGEN3_USE_FILE} )
#     endif()
# 
#     if(EIGEN3_FOUND OR LAPACK_FOUND)
#         # Link with Boost.ProgramOptions (optional)
#         find_package(Boost QUIET COMPONENTS program_options)
#         if(Boost_PROGRAM_OPTIONS_FOUND)
#             if( CGAL_AUTO_LINK_ENABLED )
#                 message( STATUS "Boost.ProgramOptions library: found" )
#             else()
#                 message( STATUS "Boost.ProgramOptions library: ${Boost_PROGRAM_OPTIONS_LIBRARY}" )
#             endif()
#             add_definitions( "-DCGAL_USE_BOOST_PROGRAM_OPTIONS" )
#             if ( NOT CGAL_AUTO_LINK_ENABLED )
#                 message("FOOBAR")
#                 link_libraries( ${Boost_PROGRAM_OPTIONS_LIBRARY} )
#             endif()
#         endif()
# 
#     else(EIGEN3_FOUND OR LAPACK_FOUND)
# 
#         message(STATUS "NOTICE: This program requires either Eigen 3.1 (or greater) or LAPACK, and will not be compiled.")
# 
#     endif(EIGEN3_FOUND OR LAPACK_FOUND)

else()

    message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")

endif()

find_package(Boost QUIET COMPONENTS program_options)
if(Boost_PROGRAM_OPTIONS_FOUND)
    if( CGAL_AUTO_LINK_ENABLED )
        message( STATUS "Boost.ProgramOptions library: found" )
    else()
        message( STATUS "Boost.ProgramOptions library: ${Boost_PROGRAM_OPTIONS_LIBRARY}" )
    endif()
    add_definitions( "-DCGAL_USE_BOOST_PROGRAM_OPTIONS" )
    if ( NOT CGAL_AUTO_LINK_ENABLED )
        message("FOOBAR")
        link_libraries( ${Boost_PROGRAM_OPTIONS_LIBRARY} )
    endif()
endif()
# include helper file
include( ${CGAL_USE_FILE} )


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

    message(STATUS "This project requires the Boost library, and will not be compiled.")

    return()  

endif()

# include for local directory

# include for local package

#########################################################
# FIND LBFGS
#########################################################
find_path(LBFGS_INCLUDE_DIR lbfgs.h )
include_directories( ${LBFGS_INCLUDE_DIR} )
find_library(LBFGS_LIBRARY NAMES liblbfgs.so )

# Creating entries for target: Demons
# ############################

add_executable( Demons  ComputationalTools.cpp ComputeGeometricEntity.cpp Demons.cpp Optimization-lbfgs.cpp PolyhedralSurf.cpp SignatureManager.cpp )

add_to_cached_list( CGAL_EXECUTABLE_TARGETS Demons )

# Link the executable to CGAL and third-party libraries
target_link_libraries(Demons   ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
target_link_libraries( Demons ${LBFGS_LIBRARY} )
