#!/bin/bash

# cat: category i.e. ggf1
# funcs: functions to test bias, i.e. (bern4 pow1)
# sig: signal injection i.e. (0 1 5)

funcs=$2
cat=$1

for str in ${funcs[@]}
do
    for sig in $3
    do
	direct=$str'_'$cat'_'$sig'sig'
	if [ -d "$direct" ]; then
	    rm $direct/*
	else
	    mkdir $direct
	fi
	sed -i 's/arguments = \$(ProcId) .*$/arguments = \$(ProcId) '$str' '$cat' '$sig'/g' condor_submit.sub
	condor_submit condor_submit.sub
	echo Done $cat $str $sig signal
    done
done
