#!/bin/bash
cd ~/EOS_space/CMSSW_14_1_0_pre4/src/
cmsenv
cd ~/HtoZg_fitting/Make_combine_workspaces/
IT=$1
combineTool.py -M Impacts -d hzg_datacard_v1p4p0_cleaned.root --setParameters pdfindex_ggf1=0,pdfindex_ggf2=5,pdfindex_ggf3=5,pdfindex_ggf4=1,pdfindex_vbf1=1,pdfindex_vbf2=1,pdfindex_vbf3=0,pdfindex_vbf4=0,pdfindex_untag=0,pdfindex_vh3l=2,pdfindex_vhptmiss=0,pdfindex_tthlep=1,pdfindex_tthhad=0,MH=125.38 -m 125.38 --freezeParameters MH -n .all.impacts.unblind --setParameterRanges r=-2,4 --doFits --robustFit 1 --cminDefaultMinimizerTolerance 0.05 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_maskChannels=2 --named $(sed $(($IT+1))'!d' condor/sys_list.txt)
