#!/bin/bash

a="20160420"
b="1510"

echo "Start processingf "$a"_"$b"_Scan delay $i SS2PF" > /Volumes/qolslc/atto/labdata/2016/"$a"_Isopropanol/SS2PF_process_"$b".log

for i in $(seq -f "%02g" 0 39)
do

    mpirun -n 8 --report-bindings --map-by slot python PeakFindingScan.py --noe 1000 --nos 5000 --fpath /Volumes/Seagate\ Expansion\ Drive/labdata/2016/"$a"_Isopropanol/"$a"_"$b"_Scan/SSdelay$i/ --spath /Volumes/qolslc/atto/labdata/2016/"$a"_Isopropanol/"$a"_"$b"_Scan/SSdelay$i/

    echo "$a_$b_Scan delay $i Peak Finding calculation done" >> /Volumes/qolslc/atto/labdata/2016/"$a"_Isopropanol/SS2PF_process_"$b".log

done

echo "$a_$b_Scan Peak Finding calculation done" >> /Volumes/qolslc/atto/labdata/2016/"$a"_Isopropanol/SS2PF_process_"$b".log
