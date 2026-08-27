#!/bin/bash 
./multi_plot.py -c ggf1 -con ../Config/config_xgboost_0207.json -b 0 -t FT -yr 0 4
./multi_plot.py -c ggf2 -con ../Config/config_xgboost_0207.json -b 5 -t FT -yr 0 2
./multi_plot.py -c ggf3 -con ../Config/config_xgboost_0207.json -b 5 -t FT -yr 0.8 1.2
./multi_plot.py -c ggf4 -con ../Config/config_xgboost_0207.json -b 1 -t FT -yr 0.9 1.1

./multi_plot.py -c vbf1 -con ../Config/config_xgboost_0207.json -b 1 -t FT -yr 0 4
./multi_plot.py -c vbf2 -con ../Config/config_xgboost_0207.json -b 1 -t FT -yr 0 2
./multi_plot.py -c vbf3 -con ../Config/config_xgboost_0207.json -b 0 -t FT -yr 0.5 1.5 
./multi_plot.py -c vbf4 -con ../Config/config_xgboost_0207.json -b 0 -t FT -yr 0.8 1.2

./multi_plot.py -c vh3l -con ../Config/config_xgboost_vhtth.json -b 2 -t pruneUnstable -yr 0 2
./multi_plot.py -c untag -con ../Config/config_xgboost_vhtth.json -b 1 -t pruneUnstable -yr 0.5 1.5
./multi_plot.py -c tthlep -con ../Config/config_xgboost_vhtth.json -b 1 -t pruneUnstable -yr 0 4
./multi_plot.py -c tthhad -con ../Config/config_xgboost_vhtth.json -b 0 -t pruneUnstable -yr 0 2
./multi_plot.py -c vhptmiss -con ../Config/config_xgboost_vhtth.json -b 0 -t pruneUnstable -yr 0 2



