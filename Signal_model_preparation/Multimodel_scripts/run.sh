#!/bin/bash
for year in 2016 2016APV 2017 2018 2022 2022EE 2023BPix 2023
do
    
    for cat in ggf1 ggf2 ggf3 ggf4
    do
	for prod in ggf vbf
	do
	    ./signal_fit.py -c $cat -y $year -p $prod -plot ../ggf_plots>../ggf_logs/$cat'_'$year'_'$prod'.log'  2>&1
	    ./config_parser.py -c $cat -p $prod -y $year -log ../ggf_logs/$cat'_'$year'_'$prod'.log' -con ../../Config/config_DSCB_flat.json
	    echo Done $cat $year $prod
	done
    done
done

