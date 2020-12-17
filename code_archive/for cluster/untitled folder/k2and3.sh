#!/bin/bash

for run in {1..10}
do
    bsub -n 1 -W 04:00 -o output.txt -R "rusage[mem=4000]" "R --vanilla --slave --args '${run}' < k3.R > result.out" 
done

for run in {1..10}
do
    bsub -n 1 -W 04:00 -o output.txt -R "rusage[mem=4000]" "R --vanilla --slave --args '${run}' < k2.R > result.out" 
done



