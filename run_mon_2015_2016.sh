#!/bin/bash

tmp='./tmp'
[[ ! -d $tmp ]] && mkdir $tmp
for yr in {2015..2016}
do
    echo $yr
    nohup nice -5 python -u get_accumulated_gridded_data_mon.py $yr &> $tmp/tmp.mon_t42."$yr".txt &
done
