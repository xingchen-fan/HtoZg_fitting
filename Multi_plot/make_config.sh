#!/bin/bash
./make_config_for_toy.py -c ggf1 -i ../Unblinding_toy/higgsCombine.ggf1.MultiDimFit.mH125.root -conB ../Config/toy_config_xgboost_nodrop.json -conS ../Config/config_DSCB.json
./make_config_for_toy.py -c ggf2 -i ../Unblinding_toy/higgsCombine.ggf2.MultiDimFit.mH125.root -conB ../Config/toy_config_xgboost_nodrop.json -conS ../Config/config_DSCB.json
./make_config_for_toy.py -c ggf3 -i ../Unblinding_toy/higgsCombine.ggf3.MultiDimFit.mH125.root -conB ../Config/toy_config_xgboost_nodrop.json -conS ../Config/config_DSCB.json
./make_config_for_toy.py -c ggf4 -i ../Unblinding_toy/higgsCombine.ggf4.MultiDimFit.mH125.root -conB ../Config/toy_config_xgboost_nodrop.json -conS ../Config/config_DSCB.json
./make_config_for_toy.py -c vbf1 -i ../Unblinding_toy/higgsCombine.vbf1.MultiDimFit.mH125.root -conB ../Config/toy_config_xgboost_nodrop.json -conS ../Config/config_DSCB.json
./make_config_for_toy.py -c vbf2 -i ../Unblinding_toy/higgsCombine.vbf2.MultiDimFit.mH125.root -conB ../Config/toy_config_xgboost_nodrop.json -conS ../Config/config_DSCB.json
./make_config_for_toy.py -c vbf3 -i ../Unblinding_toy/higgsCombine.vbf3.MultiDimFit.mH125.root -conB ../Config/toy_config_xgboost_nodrop.json -conS ../Config/config_DSCB.json
./make_config_for_toy.py -c vbf4 -i ../Unblinding_toy/higgsCombine.vbf4.MultiDimFit.mH125.root -conB ../Config/toy_config_xgboost_nodrop.json -conS ../Config/config_DSCB.json


