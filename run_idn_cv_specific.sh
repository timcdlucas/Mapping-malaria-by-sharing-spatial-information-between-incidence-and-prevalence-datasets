#!/bin/bash

###
# If you need to rerun a specific fold that failed
###

if [ $# -ne 3 ]; then
  echo $0: usage: $0 cv_type model_type fold
  exit 1
fi

if [[ "$1" != @(random|spatial|spatialkeeppr) ]]; then
  echo "cv_type must be random or spatial"
  exit 1
fi

if [[ "$2" != @(points|polygons|joint) ]]; then
  echo "model_type must be points, polygons or joint"
  exit 1
fi

Rscript run_single_cv_fold.R "$3" "$1" "$2" & 


echo 'finished'
