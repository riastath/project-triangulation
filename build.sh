#!/bin/bash

echo "Removing..."
rm -rf CMakeCache.txt CMakeFiles

cgal_create_CMakeLists -s project
cmake -DCGAL_DIR=/usr/lib/CGAL

echo "Building..."
make

echo "Complete."

# be sure to run chmod +x build.sh before executing

