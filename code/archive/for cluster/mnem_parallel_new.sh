#!/bin/bash

for type in random networks cluster cluster2 cluster3
do
    for completeness in true false
    do
        if [ $type == "random" ] || [ $type == "networks" ]
        then 
            for run in {1..90}
            do
                bsub -n 1 -W 23:55 -o output.txt -R "rusage[mem=4000]" "R --vanilla --slave --args '${type}' '1' '${completeness}' '${run}' 'kmeans' 'silhouette' 'euclidean' '1' < mnem_parallel.R > result.out" 
            done
        else 
            bsub -n 32 -W 23:55 -o output.txt -R "rusage[mem=4000]" "R --vanilla --slave --args '${type}' '100' '${completeness}' '1' 'kmeans' 'silhouette' 'euclidean' '32' < mnem_parallel.R > result.out" 
        fi
    done
done


for completeness in true false
do
    for run in {1..90}
        do
            bsub -n 1 -W 04:00 -o output.txt -R "rusage[mem=4000]" "R --vanilla --slave --args '${completeness}' '${run}' < mnem_fixk.R > result.out" 
        done
done