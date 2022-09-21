#!/bin/zsh
echo "Executing Shape-Guided-Mixed-Metro-Map-Layout"

testcasename=$1
metro=$2
guide=$3
auto=$4

instance=$5
grid_factor=$6
candidate_factor=$7

maxIter=$8


if [ $# -lt 3 ]
  then
    echo "No enough arguments supplied"
    echo "Taking Defaults: "
    metro="taipei-new.txt"
    guide="heart-guide.txt"
    auto=1
    testcasename="taipei-steps"
    instance=T
    grid_factor=0.3
    candidate_factor=1
    maxIter=192
fi

cd output-optimisation
mkdir $testcasename
cd $testcasename
mkdir png
mkdir svg
cd ../ 


cd ../MetroShapes/build/bin
echo "Metro Network File: $metro"
echo "Guide Shape File: $guide" 
echo "Execute Automatic Mode: $auto"
echo "Max number of Iteraction: $maxIter" 

for i in {1..$maxIter..1};
do
  ./MetroShapes $metro $guide $testcasename $auto "$i"
  echo "##########################"
  echo "Done with iteration: $i"
  echo "##########################"
done
