#!/bin/bash
combine ggf1_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf1=0 --toysNoSystematics -t 1000 --expectSignal 0 --saveToys -m 125 --freezeParameters pdfindex_ggf1 -n .bern4.ggf1
combine ggf1_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf1=1 --toysNoSystematics -t 1000 --expectSignal 0 --saveToys -m 125 --freezeParameters pdfindex_ggf1 -n .pow1.ggf1
combine ggf1_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf1=2 --toysNoSystematics -t 1000 --expectSignal 0 --saveToys -m 125 --freezeParameters pdfindex_ggf1 -n .exp2.ggf1
combine ggf1_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf1=3 --toysNoSystematics -t 1000 --expectSignal 0 --saveToys -m 125 --freezeParameters pdfindex_ggf1 -n .lau2.ggf1
combine ggf1_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf1=4 --toysNoSystematics -t 1000 --expectSignal 0 --saveToys -m 125 --freezeParameters pdfindex_ggf1 -n .modg.ggf1

combine ggf2_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf2=0 --toysNoSystematics -t 1000 --expectSignal 0 --saveToys -m 125 --freezeParameters pdfindex_ggf2 -n .bern4.ggf2
combine ggf2_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf2=1 --toysNoSystematics -t 1000 --expectSignal 0 --saveToys -m 125 --freezeParameters pdfindex_ggf2 -n .bern5.ggf2
combine ggf2_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf2=2 --toysNoSystematics -t 1000 --expectSignal 0 --saveToys -m 125 --freezeParameters pdfindex_ggf2 -n .pow2.ggf2
combine ggf2_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf2=3 --toysNoSystematics -t 1000 --expectSignal 0 --saveToys -m 125 --freezeParameters pdfindex_ggf2 -n .exp2.ggf2
combine ggf2_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf2=4 --toysNoSystematics -t 1000 --expectSignal 0 --saveToys -m 125 --freezeParameters pdfindex_ggf2 -n .lau2.ggf2
combine ggf2_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf2=5 --toysNoSystematics -t 1000 --expectSignal 0 --saveToys -m 125 --freezeParameters pdfindex_ggf2 -n .modg.ggf2

combine ggf3_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf3=0 --toysNoSystematics -t 1000 --expectSignal 0 --saveToys -m 125 --freezeParameters pdfindex_ggf3 -n .bern4.ggf3
combine ggf3_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf3=1 --toysNoSystematics -t 1000 --expectSignal 0 --saveToys -m 125 --freezeParameters pdfindex_ggf3 -n .pow2.ggf3
combine ggf3_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf3=2 --toysNoSystematics -t 1000 --expectSignal 0 --saveToys -m 125 --freezeParameters pdfindex_ggf3 -n .exp2.ggf3
combine ggf3_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf3=3 --toysNoSystematics -t 1000 --expectSignal 0 --saveToys -m 125 --freezeParameters pdfindex_ggf3 -n .lau2.ggf3
combine ggf3_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf3=4 --toysNoSystematics -t 1000 --expectSignal 0 --saveToys -m 125 --freezeParameters pdfindex_ggf3 -n .lau3.ggf3
combine ggf3_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf3=5 --toysNoSystematics -t 1000 --expectSignal 0 --saveToys -m 125 --freezeParameters pdfindex_ggf3 -n .modg.ggf3

combine ggf4_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf4=0 --toysNoSystematics -t 1000 --expectSignal 0 --saveToys -m 125 --freezeParameters pdfindex_ggf4 -n .bern5.ggf4
combine ggf4_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf4=1 --toysNoSystematics -t 1000 --expectSignal 0 --saveToys -m 125 --freezeParameters pdfindex_ggf4 -n .pow2.ggf4
combine ggf4_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf4=2 --toysNoSystematics -t 1000 --expectSignal 0 --saveToys -m 125 --freezeParameters pdfindex_ggf4 -n .exp2.ggf4
combine ggf4_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf4=3 --toysNoSystematics -t 1000 --expectSignal 0 --saveToys -m 125 --freezeParameters pdfindex_ggf4 -n .lau3.ggf4



