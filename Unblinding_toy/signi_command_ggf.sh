#!/bin/bash
combine -M Significance datacard/ggf1_datacard.txt --setParameters MH=125 -m 125 #--toysFrequentist
#combine -M Significance ggf1_datacard_pow1.txt --setParameters pdfindex_ggf1=3,MH=125 -n .run2mu.ggf1 -m 125 #--toysFrequentist

combine -M Significance datacard/ggf2_datacard.txt --setParameters MH=125 -m 125
combine -M Significance datacard/ggf3_datacard.txt --setParameters MH=125 -m 125
combine -M Significance datacard/ggf4_datacard.txt --setParameters MH=125 -m 125
#combine -M Significance ggf1234_datacard_bias.txt -t -1 --expectSignal=1 --setParameters pdfindex_ggf4=1,pdfindex_ggf3=0,pdfindex_ggf2=2,pdfindex_ggf1=3,MH=125 -m 125 




