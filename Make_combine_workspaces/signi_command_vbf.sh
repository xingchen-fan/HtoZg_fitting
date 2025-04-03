#!/bin/bash
combine -M Significance vbf1_datacard_bias.txt -t -1 --expectSignal 1 --setParameters pdfindex_vbf1=0,MH=125 -n .vbf1 -m 125 #--toysFrequentist
#combine -M Significance ggf1_datacard_pow1.txt --setParameters pdfindex_ggf1=3,MH=125 -n .ggf1 -m 125 #--toysFrequentist

combine -M Significance vbf2_datacard_bias.txt -t -1 --expectSignal=1 --setParameters pdfindex_vbf2=0,MH=125 -n .vbf2 -m 125
combine -M Significance vbf3_datacard_bias.txt -t -1 --expectSignal=1 --setParameters pdfindex_vbf3=3,MH=125 -n .vbf3 -m 125
combine -M Significance vbf4_datacard_bias.txt -t -1 --expectSignal=1 --setParameters pdfindex_vbf4=0,MH=125 -n .vbf4 -m 125
combine -M Significance vbf1234_datacard_bias.txt -t -1 --expectSignal=1 --setParameters pdfindex_vbf4=0,pdfindex_vbf3=3,pdfindex_vbf2=0,pdfindex_vbf1=0,MH=125 -n .vbf1234 -m 125 




