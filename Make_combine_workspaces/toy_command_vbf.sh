#!/bin/bash

combine vbf1_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf1=0 --toysNoSystematics -t 1000 --expectSignal 5 --saveToys -m 125 --freezeParameters pdfindex_vbf1 -n .5sig.pow1.vbf1
combine vbf1_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf1=1 --toysNoSystematics -t 1000 --expectSignal 5 --saveToys -m 125 --freezeParameters pdfindex_vbf1 -n .5sig.pow2.vbf1
combine vbf1_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf1=2 --toysNoSystematics -t 1000 --expectSignal 5 --saveToys -m 125 --freezeParameters pdfindex_vbf1 -n .5sig.exp2.vbf1
combine vbf1_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf1=3 --toysNoSystematics -t 1000 --expectSignal 5 --saveToys -m 125 --freezeParameters pdfindex_vbf1 -n .5sig.lau3.vbf1

combine vbf2_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf2=0 --toysNoSystematics -t 1000 --expectSignal 5 --saveToys -m 125 --freezeParameters pdfindex_vbf2 -n .5sig.bern4.vbf2
combine vbf2_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf2=1 --toysNoSystematics -t 1000 --expectSignal 5 --saveToys -m 125 --freezeParameters pdfindex_vbf2 -n .5sig.pow2.vbf2
combine vbf2_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf2=2 --toysNoSystematics -t 1000 --expectSignal 5 --saveToys -m 125 --freezeParameters pdfindex_vbf2 -n .5sig.exp2.vbf2
combine vbf2_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf2=3 --toysNoSystematics -t 1000 --expectSignal 5 --saveToys -m 125 --freezeParameters pdfindex_vbf2 -n .5sig.lau3.vbf2
combine vbf2_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf2=4 --toysNoSystematics -t 1000 --expectSignal 5 --saveToys -m 125 --freezeParameters pdfindex_vbf2 -n .5sig.modg.vbf2


combine vbf3_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf3=0 --toysNoSystematics -t 1000 --expectSignal 5 --saveToys -m 125 --freezeParameters pdfindex_vbf3 -n .5sig.bern3.vbf3
combine vbf3_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf3=1 --toysNoSystematics -t 1000 --expectSignal 5 --saveToys -m 125 --freezeParameters pdfindex_vbf3 -n .5sig.pow1.vbf3
combine vbf3_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf3=2 --toysNoSystematics -t 1000 --expectSignal 5 --saveToys -m 125 --freezeParameters pdfindex_vbf3 -n .5sig.exp2.vbf3
combine vbf3_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf3=3 --toysNoSystematics -t 1000 --expectSignal 5 --saveToys -m 125 --freezeParameters pdfindex_vbf3 -n .5sig.lau2.vbf3
combine vbf3_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf3=4 --toysNoSystematics -t 1000 --expectSignal 5 --saveToys -m 125 --freezeParameters pdfindex_vbf3 -n .5sig.modg.vbf3

combine vbf4_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf4=0 --toysNoSystematics -t 1000 --expectSignal 5 --saveToys -m 125 --freezeParameters pdfindex_vbf4 -n .5sig.bern2.vbf4
combine vbf4_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf4=1 --toysNoSystematics -t 1000 --expectSignal 5 --saveToys -m 125 --freezeParameters pdfindex_vbf4 -n .5sig.pow1.vbf4
combine vbf4_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf4=2 --toysNoSystematics -t 1000 --expectSignal 5 --saveToys -m 125 --freezeParameters pdfindex_vbf4 -n .5sig.lau2.vbf4
combine vbf4_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf4=3 --toysNoSystematics -t 1000 --expectSignal 5 --saveToys -m 125 --freezeParameters pdfindex_vbf4 -n .5sig.modg.vbf4



combine vbf1_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf1=0 --toysNoSystematics -t 1000 --expectSignal 1 --saveToys -m 125 --freezeParameters pdfindex_vbf1 -n .1sig.pow1.vbf1
combine vbf1_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf1=1 --toysNoSystematics -t 1000 --expectSignal 1 --saveToys -m 125 --freezeParameters pdfindex_vbf1 -n .1sig.pow2.vbf1
combine vbf1_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf1=2 --toysNoSystematics -t 1000 --expectSignal 1 --saveToys -m 125 --freezeParameters pdfindex_vbf1 -n .1sig.exp2.vbf1
combine vbf1_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf1=3 --toysNoSystematics -t 1000 --expectSignal 1 --saveToys -m 125 --freezeParameters pdfindex_vbf1 -n .1sig.lau3.vbf1

combine vbf2_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf2=0 --toysNoSystematics -t 1000 --expectSignal 1 --saveToys -m 125 --freezeParameters pdfindex_vbf2 -n .1sig.bern4.vbf2
combine vbf2_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf2=1 --toysNoSystematics -t 1000 --expectSignal 1 --saveToys -m 125 --freezeParameters pdfindex_vbf2 -n .1sig.pow2.vbf2
combine vbf2_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf2=2 --toysNoSystematics -t 1000 --expectSignal 1 --saveToys -m 125 --freezeParameters pdfindex_vbf2 -n .1sig.exp2.vbf2
combine vbf2_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf2=3 --toysNoSystematics -t 1000 --expectSignal 1 --saveToys -m 125 --freezeParameters pdfindex_vbf2 -n .1sig.lau3.vbf2
combine vbf2_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf2=4 --toysNoSystematics -t 1000 --expectSignal 1 --saveToys -m 125 --freezeParameters pdfindex_vbf2 -n .1sig.modg.vbf2


combine vbf3_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf3=0 --toysNoSystematics -t 1000 --expectSignal 1 --saveToys -m 125 --freezeParameters pdfindex_vbf3 -n .1sig.bern3.vbf3
combine vbf3_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf3=1 --toysNoSystematics -t 1000 --expectSignal 1 --saveToys -m 125 --freezeParameters pdfindex_vbf3 -n .1sig.pow1.vbf3
combine vbf3_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf3=2 --toysNoSystematics -t 1000 --expectSignal 1 --saveToys -m 125 --freezeParameters pdfindex_vbf3 -n .1sig.exp2.vbf3
combine vbf3_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf3=3 --toysNoSystematics -t 1000 --expectSignal 1 --saveToys -m 125 --freezeParameters pdfindex_vbf3 -n .1sig.lau2.vbf3
combine vbf3_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf3=4 --toysNoSystematics -t 1000 --expectSignal 1 --saveToys -m 125 --freezeParameters pdfindex_vbf3 -n .1sig.modg.vbf3

combine vbf4_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf4=0 --toysNoSystematics -t 1000 --expectSignal 1 --saveToys -m 125 --freezeParameters pdfindex_vbf4 -n .1sig.bern2.vbf4
combine vbf4_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf4=1 --toysNoSystematics -t 1000 --expectSignal 1 --saveToys -m 125 --freezeParameters pdfindex_vbf4 -n .1sig.pow1.vbf4
combine vbf4_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf4=2 --toysNoSystematics -t 1000 --expectSignal 1 --saveToys -m 125 --freezeParameters pdfindex_vbf4 -n .1sig.lau2.vbf4
combine vbf4_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf4=3 --toysNoSystematics -t 1000 --expectSignal 1 --saveToys -m 125 --freezeParameters pdfindex_vbf4 -n .1sig.modg.vbf4

combine vbf1_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf1=0 --toysNoSystematics -t 1000 --expectSignal 10 --saveToys -m 125 --freezeParameters pdfindex_vbf1 -n .10sig.pow1.vbf1
combine vbf1_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf1=1 --toysNoSystematics -t 1000 --expectSignal 10 --saveToys -m 125 --freezeParameters pdfindex_vbf1 -n .10sig.pow2.vbf1
combine vbf1_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf1=2 --toysNoSystematics -t 1000 --expectSignal 10 --saveToys -m 125 --freezeParameters pdfindex_vbf1 -n .10sig.exp2.vbf1
combine vbf1_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf1=3 --toysNoSystematics -t 1000 --expectSignal 10 --saveToys -m 125 --freezeParameters pdfindex_vbf1 -n .10sig.lau3.vbf1

combine vbf2_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf2=0 --toysNoSystematics -t 1000 --expectSignal 10 --saveToys -m 125 --freezeParameters pdfindex_vbf2 -n .10sig.bern4.vbf2
combine vbf2_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf2=1 --toysNoSystematics -t 1000 --expectSignal 10 --saveToys -m 125 --freezeParameters pdfindex_vbf2 -n .10sig.pow2.vbf2
combine vbf2_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf2=2 --toysNoSystematics -t 1000 --expectSignal 10 --saveToys -m 125 --freezeParameters pdfindex_vbf2 -n .10sig.exp2.vbf2
combine vbf2_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf2=3 --toysNoSystematics -t 1000 --expectSignal 10 --saveToys -m 125 --freezeParameters pdfindex_vbf2 -n .10sig.lau3.vbf2
combine vbf2_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf2=4 --toysNoSystematics -t 1000 --expectSignal 10 --saveToys -m 125 --freezeParameters pdfindex_vbf2 -n .10sig.modg.vbf2


combine vbf3_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf3=0 --toysNoSystematics -t 1000 --expectSignal 10 --saveToys -m 125 --freezeParameters pdfindex_vbf3 -n .10sig.bern3.vbf3
combine vbf3_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf3=1 --toysNoSystematics -t 1000 --expectSignal 10 --saveToys -m 125 --freezeParameters pdfindex_vbf3 -n .10sig.pow1.vbf3
combine vbf3_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf3=2 --toysNoSystematics -t 1000 --expectSignal 10 --saveToys -m 125 --freezeParameters pdfindex_vbf3 -n .10sig.exp2.vbf3
combine vbf3_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf3=3 --toysNoSystematics -t 1000 --expectSignal 10 --saveToys -m 125 --freezeParameters pdfindex_vbf3 -n .10sig.lau2.vbf3
combine vbf3_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf3=4 --toysNoSystematics -t 1000 --expectSignal 10 --saveToys -m 125 --freezeParameters pdfindex_vbf3 -n .10sig.modg.vbf3

combine vbf4_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf4=0 --toysNoSystematics -t 1000 --expectSignal 10 --saveToys -m 125 --freezeParameters pdfindex_vbf4 -n .10sig.bern2.vbf4
combine vbf4_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf4=1 --toysNoSystematics -t 1000 --expectSignal 10 --saveToys -m 125 --freezeParameters pdfindex_vbf4 -n .10sig.pow1.vbf4
combine vbf4_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf4=2 --toysNoSystematics -t 1000 --expectSignal 10 --saveToys -m 125 --freezeParameters pdfindex_vbf4 -n .10sig.lau2.vbf4
combine vbf4_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_vbf4=3 --toysNoSystematics -t 1000 --expectSignal 10 --saveToys -m 125 --freezeParameters pdfindex_vbf4 -n .10sig.modg.vbf4
