#!/usr/bin/env python3
import sys
import os
import ROOT
import argparse
import json
sys.path.append(os.path.abspath("../Utilities/"))
sys.path.append(os.path.abspath("../CMS_plotter/"))
import CMS_lumi, tdrstyle
from Xc_Minimizer import *
from plot_utility import *

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--cat', help = 'Category')
parser.add_argument('-con', '--config', help = 'Configuration', default='../Config/chi2_config_xgboost_nodrop.json')
#parser.add_argument('-f', '--func', help = 'Function')
parser.add_argument('-i', '--input', help = 'Input Root File')
#parser.add_argument('-l', '--lep', help = 'Lepton', default='')

args = parser.parse_args()
CAT = args.cat
jfile = open(args.config, 'r')
configs = json.load(jfile)
setting = configs[CAT]

f = ROOT.TFile(args.input)
w = f.Get("w")
lowx = setting["Range"]
nbins = 65
binning = ROOT.RooFit.Binning(nbins,lowx,lowx+65)

x = w.var('CMS_hzg_mass_'+CAT)
x.setRange('left', lowx, 120)
x.setRange('right', 130, lowx+65)


data = w.data("data_obs")
sb_model = w.pdf('model_s').getPdf(CAT)

# Postfit
w.loadSnapshot("MultiDimFit")
#nbkg = w.function('n_exp_final_bin'+CAT+'_proc_bkg').getVal()
#print("post fit nbkg = ", nbkg)
x.setBins(260)
hdata = ROOT.RooDataHist('hdata', 'hdata', x, data)
for i in range(260):
    hdata.get(i)
    print(hdata.weight())
#chi2 = ROOT.RooChi2Var('chi2', 'chi2', sb_model, hdata)
#val_chi2 = chi2.getVal()
#x.setBins(65)

plotClass(x, hdata, sb_model, sb_model, title='Unblind3_'+CAT, fullRatio = True, toy=True)
