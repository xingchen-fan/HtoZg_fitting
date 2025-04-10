#!/bin/bash
#function run_combine{
#  combine -M Significance ./datacards/$1_datacard_bias.txt -t -1 --expectSignal 1 --setParameters pdfindex_$1=0,MH=125 -n .vbf1 -m 125 #--toysFrequentist
#  combine -M Significance ./datacards/$1_datacard_bias.txt -t -1 --expectSignal 1 --setParameters pdfindex_$1=1,MH=125 -n .vbf1 -m 125 #--toysFrequentist
#  combine -M Significance ./datacards/$1_datacard_bias.txt -t -1 --expectSignal 1 --setParameters pdfindex_$1=2,MH=125 -n .vbf1 -m 125 #--toysFrequentist
#  combine -M Significance ./datacards/$1_datacard_bias.txt -t -1 --expectSignal 1 --setParameters pdfindex_$1=3,MH=125 -n .vbf1 -m 125 #--toysFrequentist
#  combine -M Significance ./datacards/$1_datacard_bias.txt -t -1 --expectSignal 1 --setParameters pdfindex_$1=4,MH=125 -n .vbf1 -m 125 #--toysFrequentist
#}


categories=("ttH_lep" "ttH_had" "WH_3l" "ZH_MET" "untagged")
#num_pdf_in_envelope=(5 5 6 5 6)

for cat in ${categories[@]}
do
  #pdfindex_${cat}=${ipdf}
  echo "${cat} SIGNIFICANCE -----------------------------------------------------------------------------------------------------------------------------------"
  
  echo "FUNC 0 ------------------------------------------------------------------------------------"
  combine -M Significance ./datacards/${cat}_datacard_bias.txt -t -1 --expectSignal 1 --setParameters pdfindex_${cat}=0,MH=125 -n .${cat}0 -m 125 #--toysFrequentist
  
  echo "FUNC 1 ------------------------------------------------------------------------------------"
  combine -M Significance ./datacards/${cat}_datacard_bias.txt -t -1 --expectSignal 1 --setParameters pdfindex_${cat}=1,MH=125 -n .${cat}1 -m 125 #--toysFrequentist
  
  echo "FUNC 2 ------------------------------------------------------------------------------------"
  combine -M Significance ./datacards/${cat}_datacard_bias.txt -t -1 --expectSignal 1 --setParameters pdfindex_${cat}=2,MH=125 -n .${cat}2 -m 125 #--toysFrequentist
 
  echo "FUNC 3 ------------------------------------------------------------------------------------"
  combine -M Significance ./datacards/${cat}_datacard_bias.txt -t -1 --expectSignal 1 --setParameters pdfindex_${cat}=3,MH=125 -n .${cat}3 -m 125 #--toysFrequentist
  
  echo "FUNC 4 ------------------------------------------------------------------------------------"
  combine -M Significance ./datacards/${cat}_datacard_bias.txt -t -1 --expectSignal 1 --setParameters pdfindex_${cat}=4,MH=125 -n .${cat}4 -m 125 #--toysFrequentist
  
  #nohup $(run_combine ${cat}) > out/output_significance_${cat}.txt &
  echo "${cat} SIGNIFICANCE -----------------------------------------------------------------------------------------------------------------------------------"

done

echo "FUNC 5 ------------------------------------------------------------------------------------"
combine -M Significance ./datacards/untagged_datacard_bias.txt -t -1 --expectSignal 1 --setParameters pdfindex_untagged=5,MH=125 -n .untagged5 -m 125 #--toysFrequentist

echo "FUNC 6 ------------------------------------------------------------------------------------"
combine -M Significance ./datacards/untagged_datacard_bias.txt -t -1 --expectSignal 1 --setParameters pdfindex_untagged=6,MH=125 -n .untagged6 -m 125 #--toysFrequentist


echo "WH3l FUNC 5 ------------------------------------------------------------------------------------"
combine -M Significance ./datacards/WH_3l_datacard_bias.txt -t -1 --expectSignal 1 --setParameters pdfindex_WH_3l=5,MH=125 -n .WH_3l5 -m 125 #--toysFrequentist
 

#combine -M Significance vbf1234_datacard_bias.txt -t -1 --expectSignal=1 --setParameters pdfindex_vbf4=0,pdfindex_vbf3=3,pdfindex_vbf2=0,pdfindex_vbf1=0,MH=125 -n .vbf1234 -m 125 



#










#combine -M Significance vbf1_datacard_bias.txt -t -1 --expectSignal 1 --setParameters pdfindex_vbf1=0,MH=125 -n .vbf1 -m 125 #--toysFrequentist
#combine -M Significance ggf1_datacard_pow1.txt --setParameters pdfindex_ggf1=3,MH=125 -n .ggf1 -m 125 #--toysFrequentist

#combine -M Significance vbf2_datacard_bias.txt -t -1 --expectSignal=1 --setParameters pdfindex_vbf2=0,MH=125 -n .vbf2 -m 125
#combine -M Significance vbf3_datacard_bias.txt -t -1 --expectSignal=1 --setParameters pdfindex_vbf3=3,MH=125 -n .vbf3 -m 125
#combine -M Significance vbf4_datacard_bias.txt -t -1 --expectSignal=1 --setParameters pdfindex_vbf4=0,MH=125 -n .vbf4 -m 125
#combine -M Significance vbf1234_datacard_bias.txt -t -1 --expectSignal=1 --setParameters pdfindex_vbf4=0,pdfindex_vbf3=3,pdfindex_vbf2=0,pdfindex_vbf1=0,MH=125 -n .vbf1234 -m 125 




