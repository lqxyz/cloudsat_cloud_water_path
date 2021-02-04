#!/bin/bash

tmp='./tmp'
[[ ! -d $tmp ]] && mkdir $tmp
for yr in {2012..2014}
do
    echo $yr
    nohup nice -5 python -u cloudsat_spatial_pattern_year_flip_lon_t42.py $yr &> $tmp/tmp.yr_t42."$yr".txt &
done

