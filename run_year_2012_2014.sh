#!/bin/bash

tmp='./tmp'
[[ ! -d $tmp ]] && mkdir $tmp
for yr in {2012..2014}
do
    echo $yr
    nohup nice -5 python -u get_accumulated_gridded_data_year.py $yr &> $tmp/tmp.yr_t42."$yr".txt &
done
