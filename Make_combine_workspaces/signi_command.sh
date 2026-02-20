#!/bin/bash
#combine -M Significance ggf1_datacard_bias.txt -t -1 --expectSignal 1 --setParameters pdfindex_ggf1=3,MH=125 -n .ggf1 -m 125 #--toysFrequentist
#combine -M Significance ggf1_datacard_pow1.txt --setParameters pdfindex_ggf1=1,MH=125 -n .ggf1 -m 125 #--toysFrequentist

#combine -M Significance ggf2_datacard_bias.txt -t -1 --expectSignal=1 --setParameters pdfindex_ggf2=4,MH=125 -n .ggf2 -m 125
#combine -M Significance ggf3_datacard_bias.txt -t -1 --expectSignal=1 --setParameters pdfindex_ggf3=4,MH=125 -n .ggf3 -m 125
#combine -M Significance ggf4_datacard_bias.txt -t -1 --expectSignal=1 --setParameters pdfindex_ggf4=1,MH=125 -n .ggf4 -m 125
#combine -M Significance ggf1234_datacard_bias.txt -t -1 --expectSignal=1 --setParameters pdfindex_ggf4=1,pdfindex_ggf3=4,pdfindex_ggf2=4,pdfindex_ggf1=3,MH=125 -n .ggf1234 -m 125 
combine -M Significance hzg_datacard_v1p2_cleaned.root -m 125.38 --freezeParameters MH,pdfindex_ggf3 --cminDefaultMinimizerTolerance 0.05 -t -1 --expectSignal=1 --setParameters pdfindex_ggf1=3,pdfindex_ggf2=0,pdfindex_ggf3=4,pdfindex_ggf4=4,pdfindex_vbf1=1,pdfindex_vbf2=5,pdfindex_vbf3=0,pdfindex_vbf4=2,pdfindex_vh3l=0,pdfindex_tthlep=0,pdfindex_tthhad=0,pdfindex_untag=0,pdfindex_vhptmiss=0,MH=125.38,mask_cat_vh3l=0,mask_cat_vhptmiss=0,mask_cat_tthlep=0,mask_cat_tthhad=0,mask_cat_untag=0,mask_cat_vbf1=1,mask_cat_vbf2=1,mask_cat_vbf3=1,mask_cat_vbf4=1,mask_cat_ggf4=1,mask_cat_ggf3=1,mask_cat_ggf1=1,mask_cat_ggf2=1 --setParameterRanges r=-1,3



