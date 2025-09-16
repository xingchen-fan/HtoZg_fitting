#!/bin/bash
cd /eos/user/f/fanx/CMSSW_14_1_0_pre4/src
#source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsenv
cd ~/HtoZg_fitting/Bias_test/condor/

./Bias_test_combineToy_condor.py -n 5 -i $1 -f $2 -c $3 -s $4 -conB ../../Config/config_xgboost_final_final.json -conS ../../Config/config_DSCB_double_flat.json > $2_$3_$4sig/output$1.txt 2>&1
#./Bias_test_combineToy_condor_Anders_untag.py -n 5 -i $1 -f /eos/home-f/fanx/anders_troubleshooting/higgsCombine.0sig.pdfidx1.untagged.ntoys1000.GenerateOnly.mH125.1932494908.root -c $3 -s $4 -conB ../../Config/config_xgboost_final_final.json -conS ../../Config/config_DSCB_double_flat.json > $2_$3_$4sig/output$1.txt 2>&1
