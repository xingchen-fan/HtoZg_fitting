#!/bin/bash

categories=("ttH_lep" "WH_3l" "ZH_MET" "ttH_had" "VBF")
sig_strength=0

while [ $sig_strength -le 10 ]
do
  for cat in ${categories[@]}
  do
    nohup root -l -b -q "APFitting/signal_extraction_with_generated_data.cxx(\"${cat}\",${sig_strength})" > txt_output/output_sig_inj_${cat}_${sig_strength}.txt &
  done

  wait $!

  sig_strength=$(( $sig_strength + 1 ))
done
