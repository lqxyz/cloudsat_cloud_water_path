#!/bin/bash

tmp='./tmp'
[[ ! -d $tmp ]] && mkdir $tmp
for yr in {2015..2016}
do
    echo $yr
    nohup nice -5 python -u cloudsat_spatial_pattern_mon_flip_lon_t42.py $yr &> $tmp/tmp.mon_t42."$yr".txt &
done

