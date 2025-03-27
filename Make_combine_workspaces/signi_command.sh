#!/bin/bash
#combine -M Significance ggf1_datacard_bias.txt -t -1 --expectSignal 1 --setParameters pdfindex_ggf1=3,MH=125 -n .ggf1 -m 125 #--toysFrequentist
#combine -M Significance ggf1_datacard_pow1.txt --setParameters pdfindex_ggf1=1,MH=125 -n .ggf1 -m 125 #--toysFrequentist

#combine -M Significance ggf2_datacard_bias.txt -t -1 --expectSignal=1 --setParameters pdfindex_ggf2=4,MH=125 -n .ggf2 -m 125
#combine -M Significance ggf3_datacard_bias.txt -t -1 --expectSignal=1 --setParameters pdfindex_ggf3=4,MH=125 -n .ggf3 -m 125
combine -M Significance ggf4_datacard_bias.txt -t -1 --expectSignal=1 --setParameters pdfindex_ggf4=1,MH=125 -n .ggf4 -m 125
#combine -M Significance ggf1234_datacard_bias.txt -t -1 --expectSignal=1 --setParameters pdfindex_ggf4=1,pdfindex_ggf3=4,pdfindex_ggf2=4,pdfindex_ggf1=3,MH=125 -n .ggf1234 -m 125 




