#!/bin/bash

N=10

for boot in {1..10}; do 
Rscript run_single_cv_fold.R "$boot" "random" "points" & 
  i=$(($boot % $N))
if [ $i == 0 -o $boot == 100 ]; then
wait
fi
done

echo 'finished'