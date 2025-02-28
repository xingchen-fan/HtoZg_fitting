#!/bin/bash
cd /eos/user/f/fanx/CMSSW_14_1_0_pre4/src
#source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsenv
cd ~/CMSSW_12_6_0_patch1/src/HtoZg_fitting/Bias_test/condor/
./Bias_test_combineToy_condor.py -n 5 -i $1 -f $2 -c $3> $2_$3/output$1.txt 2>&1

