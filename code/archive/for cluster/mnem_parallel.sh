#!/bin/bash

for type in random networks cluster cluster2 cluster3
do
    for completeness in true false
    do
        bsub -n 32 -W 23:55 -o output.txt -R "rusage[mem=375]" "R --vanilla --slave --args '${type}' '100' '${completeness}' '1' 'kmeans' 'silhouette' 'euclidean' '32' < mnem_parallel.R > result.out" 
    done
done

for type in cluster cluster2 cluster3
do
    for completeness in true false
    do
        bsub -n 1 -W 23:55 -o output.txt -R "rusage[mem=8000]" "R --vanilla --slave --args '${type}' '1' '${completeness}' '1' 'hc' 'silhouette' 'euclidean' '1' < mnem_parallel.R > result.out" 
    done
done

for completeness in true false
do
    bsub -n 32 -W 23:55 -o output.txt -R "rusage[mem=375]" "R --vanilla --slave --args '${completeness}' < mnem_fixk.R > result.out" 
done