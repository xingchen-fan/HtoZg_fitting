#!/bin/bash

combine ggf1_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf1=0 --toysNoSystematics -t 1000 --expectSignal 5 --saveToys -m 125 --freezeParameters pdfindex_ggf1 -n .5sig.bern4.ggf1
combine ggf1_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf1=1 --toysNoSystematics -t 1000 --expectSignal 5 --saveToys -m 125 --freezeParameters pdfindex_ggf1 -n .5sig.pow2.ggf1
combine ggf1_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf1=2 --toysNoSystematics -t 1000 --expectSignal 5 --saveToys -m 125 --freezeParameters pdfindex_ggf1 -n .5sig.exp2.ggf1
combine ggf1_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf1=3 --toysNoSystematics -t 1000 --expectSignal 5 --saveToys -m 125 --freezeParameters pdfindex_ggf1 -n .5sig.lau2.ggf1
combine ggf1_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf1=4 --toysNoSystematics -t 1000 --expectSignal 5 --saveToys -m 125 --freezeParameters pdfindex_ggf1 -n .5sig.modg.ggf1

combine ggf2_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf2=0 --toysNoSystematics -t 1000 --expectSignal 5 --saveToys -m 125 --freezeParameters pdfindex_ggf2 -n .5sig.bern5.ggf2
combine ggf2_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf2=1 --toysNoSystematics -t 1000 --expectSignal 5 --saveToys -m 125 --freezeParameters pdfindex_ggf2 -n .5sig.pow2.ggf2
combine ggf2_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf2=2 --toysNoSystematics -t 1000 --expectSignal 5 --saveToys -m 125 --freezeParameters pdfindex_ggf2 -n .5sig.exp2.ggf2
combine ggf2_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf2=3 --toysNoSystematics -t 1000 --expectSignal 5 --saveToys -m 125 --freezeParameters pdfindex_ggf2 -n .5sig.lau4.ggf2
#combine ggf2_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf2=4 --toysNoSystematics -t 1000 --expectSignal 5 --saveToys -m 125 --freezeParameters pdfindex_ggf2 -n .5sig.modg.ggf2


combine ggf3_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf3=0 --toysNoSystematics -t 1000 --expectSignal 5 --saveToys -m 125 --freezeParameters pdfindex_ggf3 -n .5sig.bern4.ggf3
combine ggf3_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf3=1 --toysNoSystematics -t 1000 --expectSignal 5 --saveToys -m 125 --freezeParameters pdfindex_ggf3 -n .5sig.pow2.ggf3
combine ggf3_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf3=2 --toysNoSystematics -t 1000 --expectSignal 5 --saveToys -m 125 --freezeParameters pdfindex_ggf3 -n .5sig.exp2.ggf3
combine ggf3_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf3=3 --toysNoSystematics -t 1000 --expectSignal 5 --saveToys -m 125 --freezeParameters pdfindex_ggf3 -n .5sig.lau3.ggf3
#combine ggf3_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf3=4 --toysNoSystematics -t 1000 --expectSignal 5 --saveToys -m 125 --freezeParameters pdfindex_ggf3 -n .5sig.modg.ggf3

combine ggf4_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf4=0 --toysNoSystematics -t 1000 --expectSignal 5 --saveToys -m 125 --freezeParameters pdfindex_ggf4 -n .5sig.pow1.ggf4
combine ggf4_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf4=1 --toysNoSystematics -t 1000 --expectSignal 5 --saveToys -m 125 --freezeParameters pdfindex_ggf4 -n .5sig.exp2.ggf4
combine ggf4_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf4=2 --toysNoSystematics -t 1000 --expectSignal 5 --saveToys -m 125 --freezeParameters pdfindex_ggf4 -n .5sig.lau2.ggf4



combine ggf1_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf1=0 --toysNoSystematics -t 1000 --expectSignal 1 --saveToys -m 125 --freezeParameters pdfindex_ggf1 -n .1sig.bern4.ggf1
combine ggf1_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf1=1 --toysNoSystematics -t 1000 --expectSignal 1 --saveToys -m 125 --freezeParameters pdfindex_ggf1 -n .1sig.pow2.ggf1
combine ggf1_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf1=2 --toysNoSystematics -t 1000 --expectSignal 1 --saveToys -m 125 --freezeParameters pdfindex_ggf1 -n .1sig.exp2.ggf1
combine ggf1_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf1=3 --toysNoSystematics -t 1000 --expectSignal 1 --saveToys -m 125 --freezeParameters pdfindex_ggf1 -n .1sig.lau2.ggf1
combine ggf1_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf1=4 --toysNoSystematics -t 1000 --expectSignal 1 --saveToys -m 125 --freezeParameters pdfindex_ggf1 -n .1sig.modg.ggf1

combine ggf2_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf2=0 --toysNoSystematics -t 1000 --expectSignal 1 --saveToys -m 125 --freezeParameters pdfindex_ggf2 -n .1sig.bern5.ggf2
combine ggf2_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf2=1 --toysNoSystematics -t 1000 --expectSignal 1 --saveToys -m 125 --freezeParameters pdfindex_ggf2 -n .1sig.pow2.ggf2
combine ggf2_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf2=2 --toysNoSystematics -t 1000 --expectSignal 1 --saveToys -m 125 --freezeParameters pdfindex_ggf2 -n .1sig.exp2.ggf2
combine ggf2_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf2=3 --toysNoSystematics -t 1000 --expectSignal 1 --saveToys -m 125 --freezeParameters pdfindex_ggf2 -n .1sig.lau4.ggf2
#combine ggf2_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf2=4 --toysNoSystematics -t 1000 --expectSignal 1 --saveToys -m 125 --freezeParameters pdfindex_ggf2 -n .1sig.modg.ggf2


combine ggf3_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf3=0 --toysNoSystematics -t 1000 --expectSignal 1 --saveToys -m 125 --freezeParameters pdfindex_ggf3 -n .1sig.bern4.ggf3
combine ggf3_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf3=1 --toysNoSystematics -t 1000 --expectSignal 1 --saveToys -m 125 --freezeParameters pdfindex_ggf3 -n .1sig.pow2.ggf3
combine ggf3_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf3=2 --toysNoSystematics -t 1000 --expectSignal 1 --saveToys -m 125 --freezeParameters pdfindex_ggf3 -n .1sig.exp2.ggf3
combine ggf3_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf3=3 --toysNoSystematics -t 1000 --expectSignal 1 --saveToys -m 125 --freezeParameters pdfindex_ggf3 -n .1sig.lau3.ggf3
#combine ggf3_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf3=4 --toysNoSystematics -t 1000 --expectSignal 1 --saveToys -m 125 --freezeParameters pdfindex_ggf3 -n .1sig.modg.ggf3

combine ggf4_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf4=0 --toysNoSystematics -t 1000 --expectSignal 1 --saveToys -m 125 --freezeParameters pdfindex_ggf4 -n .1sig.pow1.ggf4
combine ggf4_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf4=1 --toysNoSystematics -t 1000 --expectSignal 1 --saveToys -m 125 --freezeParameters pdfindex_ggf4 -n .1sig.exp2.ggf4
combine ggf4_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf4=2 --toysNoSystematics -t 1000 --expectSignal 1 --saveToys -m 125 --freezeParameters pdfindex_ggf4 -n .1sig.lau2.ggf4
#combine ggf4_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf4=3 --toysNoSystematics -t 1000 --expectSignal 1 --saveToys -m 125 --freezeParameters pdfindex_ggf4 -n .1sig.modg.ggf4

combine ggf1_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf1=0 --toysNoSystematics -t 1000 --expectSignal 10 --saveToys -m 125 --freezeParameters pdfindex_ggf1 -n .10sig.bern4.ggf1
combine ggf1_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf1=1 --toysNoSystematics -t 1000 --expectSignal 10 --saveToys -m 125 --freezeParameters pdfindex_ggf1 -n .10sig.pow2.ggf1
combine ggf1_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf1=2 --toysNoSystematics -t 1000 --expectSignal 10 --saveToys -m 125 --freezeParameters pdfindex_ggf1 -n .10sig.exp2.ggf1
combine ggf1_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf1=3 --toysNoSystematics -t 1000 --expectSignal 10 --saveToys -m 125 --freezeParameters pdfindex_ggf1 -n .10sig.lau2.ggf1
combine ggf1_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf1=4 --toysNoSystematics -t 1000 --expectSignal 10 --saveToys -m 125 --freezeParameters pdfindex_ggf1 -n .10sig.modg.ggf1

combine ggf2_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf2=0 --toysNoSystematics -t 1000 --expectSignal 10 --saveToys -m 125 --freezeParameters pdfindex_ggf2 -n .10sig.bern5.ggf2
combine ggf2_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf2=1 --toysNoSystematics -t 1000 --expectSignal 10 --saveToys -m 125 --freezeParameters pdfindex_ggf2 -n .10sig.pow2.ggf2
combine ggf2_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf2=2 --toysNoSystematics -t 1000 --expectSignal 10 --saveToys -m 125 --freezeParameters pdfindex_ggf2 -n .10sig.exp2.ggf2
combine ggf2_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf2=3 --toysNoSystematics -t 1000 --expectSignal 10 --saveToys -m 125 --freezeParameters pdfindex_ggf2 -n .10sig.lau4.ggf2
#combine ggf2_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf2=4 --toysNoSystematics -t 1000 --expectSignal 10 --saveToys -m 125 --freezeParameters pdfindex_ggf2 -n .10sig.modg.ggf2


combine ggf3_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf3=0 --toysNoSystematics -t 1000 --expectSignal 10 --saveToys -m 125 --freezeParameters pdfindex_ggf3 -n .10sig.bern4.ggf3
combine ggf3_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf3=1 --toysNoSystematics -t 1000 --expectSignal 10 --saveToys -m 125 --freezeParameters pdfindex_ggf3 -n .10sig.pow2.ggf3
combine ggf3_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf3=2 --toysNoSystematics -t 1000 --expectSignal 10 --saveToys -m 125 --freezeParameters pdfindex_ggf3 -n .10sig.exp2.ggf3
combine ggf3_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf3=3 --toysNoSystematics -t 1000 --expectSignal 10 --saveToys -m 125 --freezeParameters pdfindex_ggf3 -n .10sig.lau3.ggf3
#combine ggf3_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf3=4 --toysNoSystematics -t 1000 --expectSignal 10 --saveToys -m 125 --freezeParameters pdfindex_ggf3 -n .10sig.modg.ggf3

combine ggf4_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf4=0 --toysNoSystematics -t 1000 --expectSignal 10 --saveToys -m 125 --freezeParameters pdfindex_ggf4 -n .10sig.pow1.ggf4
combine ggf4_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf4=1 --toysNoSystematics -t 1000 --expectSignal 10 --saveToys -m 125 --freezeParameters pdfindex_ggf4 -n .10sig.exp2.ggf4
combine ggf4_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf4=2 --toysNoSystematics -t 1000 --expectSignal 10 --saveToys -m 125 --freezeParameters pdfindex_ggf4 -n .10sig.lau2.ggf4
#combine ggf4_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf4=3 --toysNoSystematics -t 1000 --expectSignal 10 --saveToys -m 125 --freezeParameters pdfindex_ggf4 -n .10sig.modg.ggf4

