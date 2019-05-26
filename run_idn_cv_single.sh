#!/bin/bash

if [ $# -ne 2 ]; then
  echo $0: usage: $0 cv_type model_type
  exit 1
fi

if [[ "$1" != @(random|spatial|spatialkeeppr) ]]; then
  echo "cv_type must be random or spatial"
  exit 1
fi

if [[ "$2" != @(points|polygons|joint|prgp) ]]; then
  echo "model_type must be points, polygons, joint or prgp"
  exit 1
fi

N=10

if [ "$1" == "random" ]; then
   END=10
elif [ "$1" == "spatial" ]; then
   END=7
elif [ "$1" == "spatialkeeppr" ]; then
   END=7
fi


for fold in $(seq 1 $END); do 
   Rscript run_single_cv_fold.R "$fold" "$1" "$2" & 
done

echo 'finished'
