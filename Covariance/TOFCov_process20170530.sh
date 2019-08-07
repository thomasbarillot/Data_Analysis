#!/bin/bash

#list="1040 1045 1050 1055 1100 1105 1110 1115 1120 1125 1129 1134 1139 1144 1149 1154 1159 1204 1209 1214"
list="1510"

for i in $list
do 
    a="20160420"
    b="$i"

    for j in $(seq -f "%02g" 0 39)
    do 
        python TOFCov_block.py --noe 1000 --nos 5000 --binning $1 --x0 $2 --y0 $3 --fpath /Volumes/Seagate\ Expansion\ Drive\ 1/labdata/2016/"$a"_Isopropanol/"$a"_"$b"_Scan/SSdelay$j/ --spath /Volumes/for_mac/007labdata/IsopropanolSSRun4/Covariance/"$a"_"$b"_CovMap/SSdelay$j/
        python TOFCov_merge_block.py --nos 5000 --binning $1 --x0 $2 --y0 $3 --fpath /Volumes/for_mac/007labdata/IsopropanolSSRun4/Covariance/"$a"_"$b"_CovMap/SSdelay$j/
	rm -f /Volumes/for_mac/007labdata/IsopropanolSSRun4/Covariance/"$a"_"$b"_CovMap/SSdelay$j/*.dat
    done
done
