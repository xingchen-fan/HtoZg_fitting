categories=("ttH_lep" "ttH_had" "WH_3l" "ZH_MET" "untagged")

for cat in ${categories[@]}
do
    nohup python3 ./VHttH_plot_envelope.py -c $cat -con vhtth_config.json > out/output_plot_$cat.txt &
done
