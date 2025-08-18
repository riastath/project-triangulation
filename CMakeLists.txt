# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.23)

project( opt_triangulation )

# CGAL and its components
find_package( CGAL QUIET COMPONENTS Qt5 )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

# include for local directory

# include for local package


# Creating entries for target: project
# ############################

add_executable(opt_triangulation algorithms/ex.cpp algorithms/pslg.cpp algorithms/opt_algorithms.cpp)
add_to_cached_list( CGAL_EXECUTABLE_TARGETS project )
target_include_directories(opt_triangulation PRIVATE algorithms)

# Link the executable to CGAL and third-party libraries
target_link_libraries(opt_triangulation PRIVATE CGAL::CGAL)

#ignore BOOST deprecated headers
add_definitions("-DBOOST_ALLOW_DEPRECATED_HEADERS")
add_definitions("-DBOOST_BIND_GLOBAL_PLACEHOLDERS")

if (CGAL_Qt5_FOUND)
  add_definitions(-DCGAL_USE_BASIC_VIEWER)
  target_link_libraries(opt_triangulation PUBLIC CGAL::CGAL CGAL::CGAL_Qt5)
endif()
