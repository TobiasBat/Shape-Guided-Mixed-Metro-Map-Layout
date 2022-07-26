#!/usr/bin/env bash
echo "Executing Shape-Guided-Mixed-Metro-Map-Layout"

metro=$1
guide=$2
auto=$3

if [ $# -lt 3 ]
  then
    echo "No enough arguments supplied"
    echo "Taking Defaults: "
    metro="taipei-new.txt"
    guide="heart-guide.txt"
    auto=1
fi

cd MetroShapes/build/bin
echo "Metro Network File: $metro"
echo "Guide Shape File: $guide" 
echo "Execute Automatic Mode: $auto"

./MetroShapes $metro $guide


