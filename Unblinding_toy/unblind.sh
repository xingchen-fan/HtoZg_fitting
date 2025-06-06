#!/bin/bash
#combine -M MultiDimFit -d datacard/ggf1_datacard.txt --algo grid --setParameterRanges r=-2,4 --cminDefaultMinimizerStrategy 0 --saveNLL -n .ggf1  -m 125  --saveWorkspace --points 30
#combine -M MultiDimFit -d datacard/ggf1_datacard_el.txt --algo grid --setParameterRanges r=1,8 --cminDefaultMinimizerStrategy 0 --saveNLL -n .ggf1.el  -m 125  --saveWorkspace --points 30
#combine -M MultiDimFit -d datacard/ggf1_datacard_mu.txt --algo grid --setParameterRanges r=-3,3 --cminDefaultMinimizerStrategy 0 --saveNLL -n .ggf1.mu  -m 125  --saveWorkspace --points 30

#combine -M MultiDimFit -d datacard/ggf2_datacard.txt --algo grid --setParameterRanges r=-1,4 --cminDefaultMinimizerStrategy 0 --saveNLL -n .ggf2  -m 125  --saveWorkspace --points 30
#combine -M MultiDimFit -d datacard/ggf2_datacard_el.txt --algo grid --setParameterRanges r=-2,3.5 --cminDefaultMinimizerStrategy 0 --saveNLL -n .ggf2.el  -m 125  --saveWorkspace --points 50  --cminDefaultMinimizerTolerance 0.05 --freezeParameters pdfindex_ggf2_el --setParameters pdfindex_ggf2_el=0
#combine -M MultiDimFit -d datacard/ggf2_datacard_mu.txt --algo grid --setParameterRanges r=1,4.5 --cminDefaultMinimizerStrategy 0 --saveNLL -n .ggf2.mu  -m 125  --saveWorkspace --points 30
#--freezeParameters pdfindex_ggf2_el --setParameters pdfindex_ggf2_el=0

#combine -M MultiDimFit -d datacard/ggf3_datacard.txt --algo grid --setParameterRanges r=-1,4 --cminDefaultMinimizerStrategy 0 --saveNLL -n .ggf3  -m 125  --saveWorkspace --points 30 --cminDefaultMinimizerTolerance 0.2
#combine -M MultiDimFit -d datacard/ggf3_datacard_el.txt --algo grid --setParameterRanges r=0,6 --cminDefaultMinimizerStrategy 0 --saveNLL -n .ggf3.el  -m 125  --saveWorkspace --points 30 --freezeParameters pdfindex_ggf4_el --setParameters pdfindex_ggf4_el=0 --cminDefaultMinimizerTolerance 0.2
#combine -M MultiDimFit -d datacard/ggf3_datacard_mu.txt --algo grid --setParameterRanges r=-0.5,5.5 --cminDefaultMinimizerStrategy 0 --saveNLL -n .ggf3.mu  -m 125  --saveWorkspace --points 30

#combine -M MultiDimFit -d datacard/ggf4_datacard.txt --algo grid --setParameterRanges r=-3,6.5 --cminDefaultMinimizerStrategy 0 --saveNLL -n .ggf4  -m 125  --saveWorkspace --points 30 --cminDefaultMinimizerTolerance 0.05
#combine -M MultiDimFit -d datacard/ggf4_datacard_el.txt --algo grid --setParameterRanges r=-6,7 --cminDefaultMinimizerStrategy 0 --saveNLL -n .ggf4.el  -m 125  --saveWorkspace --points 30 --cminDefaultMinimizerTolerance 0.05
#combine -M MultiDimFit -d datacard/ggf4_datacard_mu.txt --algo grid --setParameterRanges r=-5,11 --cminDefaultMinimizerStrategy 0 --saveNLL -n .ggf4.mu  -m 125  --saveWorkspace --points 50 --cminDefaultMinimizerTolerance 0.05

#combine -M MultiDimFit -d datacard/ggf3_datacard.txt --cminDefaultMinimizerTolerance 0.02

#combine -M MultiDimFit -d datacard/vbf1_datacard.txt --algo grid --setParameterRanges r=0,5 --cminDefaultMinimizerStrategy 0 --saveNLL -n .vbf1  -m 125  --saveWorkspace --points 30
#combine -M MultiDimFit -d datacard/vbf2_datacard.txt --algo grid --setParameterRanges r=-1,3 --cminDefaultMinimizerStrategy 0 --saveNLL -n .vbf2  -m 125  --saveWorkspace --points 30
#combine -M MultiDimFit -d datacard/vbf3_datacard.txt --algo grid --setParameterRanges r=-3,2 --cminDefaultMinimizerStrategy 0 --saveNLL -n .vbf3  -m 125  --saveWorkspace --points 30
#combine -M MultiDimFit -d datacard/vbf4_datacard.txt --algo grid --setParameterRanges r=1,10 --cminDefaultMinimizerStrategy 0 --saveNLL -n .vbf4  -m 125  --saveWorkspace --points 30

#combine -M MultiDimFit -d datacard/vbf3_datacard_el.txt --algo grid --setParameterRanges r=-4,1 --cminDefaultMinimizerStrategy 0 --saveNLL -n .vbf3.el  -m 125  --saveWorkspace --points 30
#combine -M MultiDimFit -d datacard/vbf3_datacard_mu.txt --algo grid --setParameterRanges pow1_p_vbf3=-15,0 --rMin -1 --rMax 4 --cminDefaultMinimizerStrategy 0 --saveNLL -n .vbf3.mu  -m 125  --saveWorkspace --points 50 --freezeParameters pdfindex_vbf3_mu --setParameters pdfindex_vbf3_mu=1,pow1_p_vbf3=-5.13,pow1_sigma_vbf3=2.2338,pow1_t_vbf3=100,r=1
#--cminDefaultMinimizerTolerance 10 --cminPreFit 1 --cminApproxPreFitTolerance 10 

#combine -M MultiDimFit -d datacard/vbf4_datacard_el.txt --algo grid --setParameterRanges r=-3,11 --cminDefaultMinimizerStrategy 0 --saveNLL -n .vbf4.el  -m 125  --saveWorkspace --points 30
#combine -M MultiDimFit -d datacard/vbf4_datacard_mu.txt --algo grid --setParameterRanges r=1,14 --cminDefaultMinimizerStrategy 0 --saveNLL -n .vbf4.mu  -m 125  --saveWorkspace --points 30

#combine -M MultiDimFit -d datacard/vbf12_datacard_el.txt --algo grid --setParameterRanges r=0,4 --cminDefaultMinimizerStrategy 0 --saveNLL -n .vbf12.el  -m 125  --saveWorkspace --points 30
#combine -M MultiDimFit -d datacard/vbf12_datacard_mu.txt --algo grid --setParameterRanges r=-1,2 --cminDefaultMinimizerStrategy 0 --saveNLL -n .vbf12.mu  -m 125  --saveWorkspace --points 30
#combine -M MultiDimFit -d datacard/vbf3_datacard_mu.txt -m 125 --freezeParameters pdfindex_vbf3_mu --setParameters pdfindex_vbf3_mu=1,pow1_p_vbf3=-5.13,pow1_sigma_vbf3=2.2338,pow1_t_vbf3=100

combine -M MultiDimFit -d datacard/vbf_ggf_datacard.txt --algo grid --setParameterRanges r=-2,5 --cminDefaultMinimizerStrategy 0 --saveNLL -n .full  -m 125  --saveWorkspace --points 30
