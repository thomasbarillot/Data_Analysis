#!/bin/bash

list="0553"

for i in $list
do 
    a="20160528"
    b="$i"

    for j in $(seq -f "%02g" 0 53)
    do 
        mpirun -n 16 python iTOFCovariance_block.py --noe 1000 --nos 8000 --binning $1 --x0 $2 --y0 $3 --fpath /Volumes/Seagate\ Expansion\ Drive/labdata/2016/20160526_Isopropanol/"$a"_"$b"_Scan/SSdelay$j/ --spath /Volumes/for_mac/007labdata/IsopropanolSSRun4/Covariance/"$a"_"$b"_CovMap/SSdelay$j/
    done
done
