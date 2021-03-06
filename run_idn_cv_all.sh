#!/bin/bash

declare -a cv_type=("random" "spatial" "spatialkeeppr")
declare -a model_type=("points" "polygons" "joint" "prgp")

for i in "${cv_type[@]}" 
do
   for j in "${model_type[@]}" 
   do
      echo "Run $i $j"
      ./run_idn_cv_single.sh $i $j &
   done
done

