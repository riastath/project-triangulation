# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.23)

project( project )


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
# added these two below
add_subdirectory(includes/utils)
add_subdirectory(includes/graphics)
list(APPEND EXTRA_LIBS utils)

# include for local package


# Creating entries for target: project
# ############################

add_executable( project  ex.cpp )

add_to_cached_list( CGAL_EXECUTABLE_TARGETS project )

# Link the executable to CGAL and third-party libraries
#target_link_libraries(project PRIVATE CGAL::CGAL )
#target_link_libraries(project PRIVATE CGAL::CGAL ${EXTRA_LIBS})
target_link_libraries(project Qt5::Widgets Qt5::Gui Qt5::Core CGAL::CGAL Boost::boost)

# need to find out how to make qt work (if it's needed, we can simulate otherwise i guess)
# changed subdirectories

