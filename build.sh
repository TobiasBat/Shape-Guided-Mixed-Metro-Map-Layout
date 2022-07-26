#!/usr/bin/env bash
echo "Building Shape-Guided-Mixed-Metro-Map-Layout"
cd Deformationprocess
mkdir build
cd build
cmake ../ 
make 
