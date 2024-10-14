# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.23)

project( project )

find_package(CGAL COMPONENTS Qt5)

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

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

# add_executable(test_draw tests/test_draw.cpp qrc_File.cpp qrc_CGAL.cpp qrc_Input.cpp qrc_Triangulation_2.cpp )

# add_to_cached_list( CGAL_EXECUTABLE_TARGETS test_draw )

add_executable(project  ex.cpp )
# add_executable(test_draw tests/test_draw.cpp)
add_to_cached_list( CGAL_EXECUTABLE_TARGETS project )

# Link the executable to CGAL and third-party libraries
target_link_libraries(project PRIVATE CGAL::CGAL)


# if (CGAL_Qt5_FOUND)
#   add_definitions(-DCGAL_USE_BASIC_VIEWER)
#   target_link_libraries(test_draw PUBLIC CGAL::CGAL CGAL::CGAL_Qt5)
# endif()
