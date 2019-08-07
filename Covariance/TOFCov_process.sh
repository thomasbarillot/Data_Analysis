#!/bin/bash

list="1139"

for i in $list
do 
    a="20161110"
    b="1139"

    for j in $(seq -f "%02g" 0 0)
    do 
        mpirun -n 8 python TOFCov.py --noe 1000 --nos 8000 --fpath /Volumes/d/labdata/"$a"_Pentanol/"$a"_"$b"_Scan/SSdelay00/ --spath /Volumes/for_mac/007Labdata/PentanolSSRun2/Covariance/"$a"_"$b"_CovMap
        python TOFCov_Merge.py --nos 8000 --fpath /Volumes/for_mac/007Labdata/PentanolSSRun2/Covariance/"$a"_"$b"_CovMap
        rm -f /Volumes/for_mac/007Labdata/PentanolSSRun2/Covariance/"$a"_"$b"_CovMap/*.dat
    done
done
