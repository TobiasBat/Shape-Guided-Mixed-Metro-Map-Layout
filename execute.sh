#!/usr/bin/env bash
echo "Executing Shape-Guided-Mixed-Metro-Map-Layout"

testcasename=$1
metro=$2
guide=$3
auto=$4

instance=$5
grid_factor=$6
candidate_factor=$7

if [ $# -lt 3 ]
  then
    echo "No enough arguments supplied"
    echo "Taking Defaults: "
    metro="taipei-new.txt"
    guide="heart-guide.txt"
    auto=1
    testcasename="taipei"
    instance=T
    grid_factor=0.3
    candidate_factor=1
fi

cd output-optimisation
mkdir $testcasename

cd ../MetroShapes/build/bin
echo "Metro Network File: $metro"
echo "Guide Shape File: $guide" 
echo "Execute Automatic Mode: $auto"

./MetroShapes $metro $guide $testcasename $auto

echo "Completed the Optimization Stage"
echo "Starting Ocolinearization"

cd ../../../ac-metro-schematization-framework
source .venv38_pg22/bin/activate
cd ac_metro/tests
python3 octi_metro_heuristic.py $instance $grid_factor $candidate_factor 20 1 0 1 0

echo "Octolinearization finished"