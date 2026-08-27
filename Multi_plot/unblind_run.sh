#!/bin/bash
./unblind_stage123.py -c ggf1 -con ../Config/config_xgboost_0207.json -f bern3 -i ../Make_combine_workspaces/higgsCombine.all.unblind1.MultiDimFit.mH125.38.root -s 3 -eb True -ed ~/EOS_space/eb_toys_5k
./unblind_stage123.py -c ggf2 -con ../Config/config_xgboost_0207.json -f agg -i ../Make_combine_workspaces/higgsCombine.all.unblind1.MultiDimFit.mH125.38.root -s 3 -eb True -ed ~/EOS_space/eb_toys_5k
./unblind_stage123.py -c ggf3 -con ../Config/config_xgboost_0207.json -f lau3 -i ../Make_combine_workspaces/higgsCombine.all.unblind1.MultiDimFit.mH125.38.root -s 3 -eb True -ed ~/EOS_space/eb_toys_5k
./unblind_stage123.py -c ggf4 -con ../Config/config_xgboost_0207.json -f exp2 -i ../Make_combine_workspaces/higgsCombine.all.unblind1.MultiDimFit.mH125.38.root -s 3 -eb True -ed ~/EOS_space/eb_toys_5k

./unblind_stage123.py -c vbf1 -con ../Config/config_xgboost_0207.json -f exp1 -i ../Make_combine_workspaces/higgsCombine.all.unblind1.MultiDimFit.mH125.38.root -s 3 -eb True -ed ~/EOS_space/eb_toys_5k
./unblind_stage123.py -c vbf2 -con ../Config/config_xgboost_0207.json -f bern4 -i ../Make_combine_workspaces/higgsCombine.all.unblind1.MultiDimFit.mH125.38.root -s 3 -eb True -ed ~/EOS_space/eb_toys_5k
./unblind_stage123.py -c vbf3 -con ../Config/config_xgboost_0207.json -f modg -i ../Make_combine_workspaces/higgsCombine.all.unblind1.MultiDimFit.mH125.38.root -s 3 -eb True -ed ~/EOS_space/eb_toys_5k
./unblind_stage123.py -c vbf4 -con ../Config/config_xgboost_0207.json -f bern2 -i ../Make_combine_workspaces/higgsCombine.all.unblind1.MultiDimFit.mH125.38.root -s 3 -eb True -ed ~/EOS_space/eb_toys_5k

./unblind_stage123.py -c vh3l -con ../Config/config_xgboost_vhtth.json -f pow1 -i ../Make_combine_workspaces/higgsCombine.all.unblind1.MultiDimFit.mH125.38.root -s 3 -eb True -ed ~/EOS_space/eb_toys_5k
./unblind_stage123.py -c untag -con ../Config/config_xgboost_vhtth.json -f pow1 -i ../Make_combine_workspaces/higgsCombine.all.unblind1.MultiDimFit.mH125.38.root -s 3 -eb True -ed ~/EOS_space/eb_toys_5k
./unblind_stage123.py -c tthlep -con ../Config/config_xgboost_vhtth.json -f pow1 -i ../Make_combine_workspaces/higgsCombine.all.unblind1.MultiDimFit.mH125.38.root -s 3 -eb True -ed ~/EOS_space/eb_toys_5k
./unblind_stage123.py -c tthhad -con ../Config/config_xgboost_vhtth.json -f bern2 -i ../Make_combine_workspaces/higgsCombine.all.unblind1.MultiDimFit.mH125.38.root -s 3 -eb True -ed ~/EOS_space/eb_toys_5k
./unblind_stage123.py -c vhptmiss -con ../Config/config_xgboost_vhtth.json -f bern2 -i ../Make_combine_workspaces/higgsCombine.all.unblind1.MultiDimFit.mH125.38.root -s 3 -eb True -ed ~/EOS_space/eb_toys_5k



