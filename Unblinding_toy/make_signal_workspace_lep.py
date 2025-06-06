#!/usr/bin/env python3
import math
import ROOT
import os
import json
import sys
import argparse
sys.path.append(os.path.abspath("../Utilities/"))
from bkg_functions_class import *
from sig_functions_class import *
from Xc_Minimizer import *
from plot_utility import *
from sample_reader import *
from profile_class import *
from sig_functions_class import *
parser = argparse.ArgumentParser(description = "Make signal workspace")
parser.add_argument('-c', '--cat', help="category")
parser.add_argument('-conB', '--configB', help = 'Bkg Configuration', default='../Config/chi2_config_xgboost_nodrop.json')
parser.add_argument('-conS', '--configS', help = 'Sig Configuration', default='../Config/config_DSCB.json')

args = parser.parse_args()
jfile = open(args.configB, 'r')
configs = json.load(jfile)
CAT = args.cat
setting = configs[CAT]
lowx = setting["Range"]

x = ROOT.RooRealVar("CMS_hzg_mass_"+CAT, "CMS_hzg_mass_"+CAT, lowx, lowx + 65.)
x.setBins(260)
MH = ROOT.RooRealVar("MH","MH"       ,125)#, 120., 130.)
sig_model_el = combineSignalLep(x, MH, CAT, args.configS, 'el')
sig_model_mu = combineSignalLep(x, MH, CAT, args.configS, 'mu')
sig_model = combineSignal(x, MH, CAT, args.configS)

print('n tot = ',sig_model.ntot)

f_out1 = ROOT.TFile("workspaces/workspace_sig_toy_unblind_" + CAT+".root", "RECREATE")
w_sig = ROOT.RooWorkspace("workspace_sig","workspace_sig")
getattr(w_sig, "import")(sig_model.pdf)
w_sig.Write()
w_sig_el = ROOT.RooWorkspace("workspace_sig_el","workspace_sig_el")
getattr(w_sig_el, "import")(sig_model_el.pdf)
w_sig_el.Write()
w_sig_mu = ROOT.RooWorkspace("workspace_sig_mu","workspace_sig_mu")
getattr(w_sig_mu, "import")(sig_model_mu.pdf)
w_sig_mu.Write()

f_out1.Close()
