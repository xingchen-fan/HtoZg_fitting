#!/bin/bash
# Script that sets up enviornment on UCSB servers

cd /net/cms37/data2/oshiro/analysis/CMSSW_14_1_0_pre4/src
. /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd -
