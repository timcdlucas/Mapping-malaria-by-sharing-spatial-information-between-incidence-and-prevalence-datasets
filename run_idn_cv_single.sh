#!/bin/bash

if [ $# -ne 2 ]; then
  echo $0: usage: $0 cv_type model_type
  exit 1
fi

N=10

for boot in {1..10}; do 
Rscript run_single_cv_fold.R "$boot" "$1" "$2" & 
  i=$(($boot % $N))
if [ $i == 0 -o $boot == 100 ]; then
wait
fi
done

echo 'finished'
