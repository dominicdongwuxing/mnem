#!/bin/bash
for type in random networks cluster cluster2 cluster3
do
    for ksel1 in kmeans hc
    do 
        for ksel2 in silhouette BIC AIC
        do
            if ! ( [ $ksel1 == 'hc' ] && [ $ksel2 != 'silhouette' ] )
            then
                for ksel3 in cor euclidean
                do
                    if [ $type == "cluster" ];
                    then
                        run=1
                        bsub -n 1 -W 23:55 -R "rusage[mem=8000]" "R --vanilla --slave --args '${type}' '10' 'True' '${run}' '${ksel1}' '${ksel2}' '${ksel3}' < mnem.R > result.out"  
                    elif [ $type == "cluster2" ] || [ $type == "cluster3" ]
                    then
                        for run in {1..10}
                        do 
                            bsub -n 1 -W 23:55 -R "rusage[mem=8000]" "R --vanilla --slave --args '${type}' '1' 'True' '${run}' '${ksel1}' '${ksel2}' '${ksel3}' < mnem.R > result.out" 
                        done 
                    else
                        for completeness in true false
                        do
                            for run in {1..10}
                            do
                                bsub -n 1 -W 23:55 -R "rusage[mem=8000]" "R --vanilla --slave --args '${type}' '1' '${completeness}' '${run}' '${ksel1}' '${ksel2}' '${ksel3}' < mnem.R > result.out"  
                            done
                        done
                    fi
                done
            fi
        done
    done
done