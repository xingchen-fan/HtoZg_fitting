#!/bin/bash
combine -M Significance datacard/vbf1_datacard.txt --setParameters MH=125 -m 125 #--toysFrequentist
#combine -M Significance vbf1_datacard_pow1.txt --setParameters pdfindex_vbf1=3,MH=125 -n .run2mu.vbf1 -m 125 #--toysFrequentist

combine -M Significance datacard/vbf2_datacard.txt --setParameters MH=125 -m 125
combine -M Significance datacard/vbf3_datacard.txt --setParameters MH=125 -m 125
combine -M Significance datacard/vbf4_datacard.txt --setParameters MH=125 -m 125
#combine -M Significance vbf1234_datacard_bias.txt -t -1 --expectSignal=1 --setParameters pdfindex_vbf4=1,pdfindex_vbf3=0,pdfindex_vbf2=2,pdfindex_vbf1=3,MH=125 -m 125 




