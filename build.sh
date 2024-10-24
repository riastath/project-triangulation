#!/bin/bash
cmake -DCGAL_DIR=/usr/lib/CGAL

echo "Building..."
make

echo "Complete."

# be sure to run chmod +x build.sh before executing

