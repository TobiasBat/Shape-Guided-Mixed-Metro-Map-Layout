#!/usr/bin/env bash
echo "Executing Shape-Guided-Mixed-Metro-Map-Layout"

testcasename=$1
metro=$2
guide=$3
auto=$4

if [ $# -lt 3 ]
  then
    echo "No enough arguments supplied"
    echo "Taking Defaults: "
    metro="taipei-new.txt"
    guide="heart-guide.txt"
    auto=1
    testcasename="taipei"
fi

cd output-optimisation
mkdir $testcasename

cd ../MetroShapes/build/bin
echo "Metro Network File: $metro"
echo "Guide Shape File: $guide" 
echo "Execute Automatic Mode: $auto"

./MetroShapes $metro $guide $testcasename $auto

echo "Completed the Optimization Stage"

echo "TODO Soeren add command"