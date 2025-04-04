categories=("ttH_lep" "ttH_had" "WH_3l" "ZH_MET" "untagged")
polylists=("BERN" "EXP" "POW" "LAU")

for cat in ${categories[@]}
do
  for plist in ${polylists[@]}
  do
    nohup python3 ./VHttH_F_test.py -c $cat -con vhtth_config.json -l $plist > out/output_ftest_$cat_$plist.txt &
  done
done
