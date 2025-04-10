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
#specific_pdf=(3 4 3 3 3)

#For TempPruned datacards
specific_pdf=(2 1 2 2 3)


a=0
for cat in ${categories[@]}
do
  #pdfindex_${cat}=${ipdf}
  echo "${cat} SIGNIFICANCE -----------------------------------------------------------------------------------------------------------------------------------"
  
  #combine -M Significance ./datacards/${cat}_datacard_bias.txt -t -1 --expectSignal 1 --setParameters pdfindex_${cat}=${specific_pdf[$a]},MH=125 -n .${cat}${a} -m 125 #--toysFrequentist
  
  #For tempPruned datacards
  combine -M Significance ./datacards/tempPruned/${cat}_datacard_bias.txt -t -1 --expectSignal 1 --setParameters pdfindex_${cat}=${specific_pdf[$a]},MH=125 -n .${cat}${a} -m 125 #--toysFrequentist
  #combine -M Significance ./datacards/tempPruned/${cat}_datacard_bias.txt -t -1 --expectSignal=2.4 --setParameters pdfindex_${cat}=${specific_pdf[$a]},MH=125 -n .${cat}${a} -m 125 #--toysFrequentist
 

  a=$((a + 1))
done


#combine -M Significance ./datacards/tempPruned_VHttH_datacard.txt -t -1 --expectSignal=1 --setParameters pdfindex_ttH_lep=3,pdfindex_ttH_had=4,pdfindex_WH_3l=3,pdfindex_ZH_MET=3,pdfindex_untagged=3,MH=125 -n .VHttHTot -m 125 

#For tempPruned datacard
combine -M Significance ./datacards/tempPruned/VHttH_datacard.txt -t -1 --expectSignal=1 --setParameters pdfindex_ttH_lep=2,pdfindex_ttH_had=1,pdfindex_WH_3l=2,pdfindex_ZH_MET=2,pdfindex_untagged=3,MH=125 -n .VHttHTot -m 125 
#combine -M Significance ./datacards/tempPruned/VHttH_datacard.txt -t -1 --expectSignal=2.4 --setParameters pdfindex_ttH_lep=2,pdfindex_ttH_had=1,pdfindex_WH_3l=2,pdfindex_ZH_MET=2,pdfindex_untagged=3,MH=125 -n .VHttHTot -m 125 



#combine -M Significance vbf1_datacard_bias.txt -t -1 --expectSignal 1 --setParameters pdfindex_vbf1=0,MH=125 -n .vbf1 -m 125 #--toysFrequentist
#combine -M Significance ggf1_datacard_pow1.txt --setParameters pdfindex_ggf1=3,MH=125 -n .ggf1 -m 125 #--toysFrequentist

#combine -M Significance vbf2_datacard_bias.txt -t -1 --expectSignal=1 --setParameters pdfindex_vbf2=0,MH=125 -n .vbf2 -m 125
#combine -M Significance vbf3_datacard_bias.txt -t -1 --expectSignal=1 --setParameters pdfindex_vbf3=3,MH=125 -n .vbf3 -m 125
#combine -M Significance vbf4_datacard_bias.txt -t -1 --expectSignal=1 --setParameters pdfindex_vbf4=0,MH=125 -n .vbf4 -m 125
#combine -M Significance vbf1234_datacard_bias.txt -t -1 --expectSignal=1 --setParameters pdfindex_vbf4=0,pdfindex_vbf3=3,pdfindex_vbf2=0,pdfindex_vbf1=0,MH=125 -n .vbf1234 -m 125 




