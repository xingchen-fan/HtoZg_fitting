#!/bin/bash
categories=("ttH_lep" "ttH_had" "WH_3l" "ZH_MET" "untagged")
num_pdf_in_envelope=(5 5 5 5 5)
signal_inj=0
num_toys=1000
a=0

for cat in ${categories[@]}
do
  npdf=${num_pdf_in_envelope[${a}]}
  for ipdf in $(seq 0 1 $npdf)
  do
    echo "The cat is $cat and pdf idx is $ipdf" 
    #combine ./datacard/${cat}_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_${cat}=${ipdf} -t ${num_toys} --expectSignal ${signal_inj} --saveToys -m 125 --freezeParameters pdfindex_${cat} -n .0sig.pdfidx${ipdf}.${cat}
  done 

  a=$((a+1))
done

comment='
combine ggf1_datacard_bias.txt -M GenerateOnly --setParameters pdfindex_ggf1=0 --toysNoSystematics -t 1000 --expectSignal 5 --saveToys -m 125 --freezeParameters pdfindex_ggf1 -n .5sig.bern4.ggf1
'
