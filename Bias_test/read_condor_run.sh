#!/bin/bash

funcs=(bern3 pow1 exp2 lau2 lau3 agg)
#(bern3 exp2 lau2) ggf4
#(bern3 bern4 pow1 exp2 lau2) vbf2
#(bern3 pow1 exp2 lau2 lau3 agg) ggf3

cat=ggf3
for str in ${funcs[@]}
do
    for i in 0 1
    do
	root -l -q 'read_condor_output.c("'$str'", "'$cat'", '$i')'
    done
done
