#!/usr/bin/env bash
echo "Building Shape-Guided-Mixed-Metro-Map-Layout"
cd MetroShapes
mkdir build
cd build
cmake ../ 
make 
