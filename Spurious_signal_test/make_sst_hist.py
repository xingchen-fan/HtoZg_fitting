#!/usr/bin/env python3
import os
import sys
import ROOT
import argparse
import json
sys.path.append(os.path.abspath("../Utilities/"))
sys.path.append(os.path.abspath("../CMS_plotter/"))
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--cat', help = 'Category')
args = parser.parse_args()
CAT = args.cat

year = ['2016', '2016APV', '2017', '2018', '2022', '2022EE', '2023', '2023BPix']
prod = ['DY0', 'DYfilter', 'DYmix', 'SM1', 'SM2']
chain = ROOT.TChain('outtree')
if CAT == 'ggf':
    bdt1 = 0.94
    bdt2 = 0.83
    bdt3 = 0.57
    dir_ = '/eos/project/h/htozg-dy-privatemc/rzou/bdt/BDT_output_redwood/Output_ggF_rui_redwood_v1_ext_val/'
elif CAT == 'vbf':
    bdt1 = 0.91
    bdt2 = 0.81
    bdt3 = 0.48
    dir_ = '/eos/project/h/htozg-dy-privatemc/rzou/bdt/BDT_output_redwood/Output_VBF_rui_redwood_v1_ext_val/'
else:
    sys.exit()
    
hist1_TH1 = ROOT.TH1F(f'{CAT}1_dy_sm', f'{CAT}1_dy_sm', 344, 94, 180)
hist2_TH1 = ROOT.TH1F(f'{CAT}2_dy_sm', f'{CAT}2_dy_sm', 344, 94, 180)
hist3_TH1 = ROOT.TH1F(f'{CAT}3_dy_sm', f'{CAT}3_dy_sm', 344, 94, 180)
hist4_TH1 = ROOT.TH1F(f'{CAT}4_dy_sm', f'{CAT}4_dy_sm', 344, 94, 180)

for y in year:
    for p in prod:
        if p == 'SM1' or p == 'SM2':
            chain.Add(dir_ +p+'_'+y+'_output.root')
        else:
            chain.Add(dir_ +p+'_'+y+'_output_1.root')
if CAT == 'ggf':
    for entry in chain:
        if entry.BDT_score > bdt1 and entry.met < 90 and entry.njet < 2:
            hist1_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)
        elif entry.BDT_score > bdt2 and entry.BDT_score < bdt1 and entry.met < 90 and entry.njet < 2:
            hist2_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)
        elif entry.BDT_score > bdt3 and entry.BDT_score < bdt2 and entry.met < 90 and entry.njet < 2:
            hist3_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)
        elif entry.BDT_score > -1 and entry.BDT_score < bdt3 and entry.met < 90 and entry.njet < 2:
            hist4_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)
            
elif CAT == 'vbf':
    for entry in chain:
        if entry.BDT_score > bdt1:
            hist1_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)
        elif entry.BDT_score > bdt2 and entry.BDT_score < bdt1:
            hist2_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)
        elif entry.BDT_score > bdt3 and entry.BDT_score < bdt2:
            hist3_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)
        elif entry.BDT_score > -1 and entry.BDT_score < bdt3:
            hist4_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)
            
f_out = ROOT.TFile("SST_" + CAT + "_hist_0207.root", "RECREATE")
hist1_TH1.Write()
hist2_TH1.Write()
hist3_TH1.Write()
hist4_TH1.Write()
f_out.Close()
