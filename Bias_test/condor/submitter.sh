#!/bin/bash

funcs=(bern2 pow1 lau2 modg)
cat=vbf4

for str in ${funcs[@]}
do
    for sig in 1 5 10
    do
	sed -i 's/arguments = \$(ProcId) .*$/arguments = \$(ProcId) '$str' '$cat' '$sig'/g' condor_submit.sub
	condor_submit condor_submit.sub
    done
done
