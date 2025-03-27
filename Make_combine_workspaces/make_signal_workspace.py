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
parser.add_argument('-con', '--config', help = 'Configuration')
args = parser.parse_args()
jfile = open('../Config/'+args.config, 'r')
configs = json.load(jfile)
CAT = args.cat
setting = configs[CAT]
lowx = setting["Range"]

x = ROOT.RooRealVar("CMS_hzg_mass_"+CAT, "CMS_hzg_mass_"+CAT, lowx, lowx + 65.)
x.setBins(260)
MH = ROOT.RooRealVar("MH","MH"       ,125)#, 120., 130.)
sig_model = combineSignal(x, MH, CAT, '../Config/config_DSCB.json')
f_out1 = ROOT.TFile("workspaces/workspace_sig_" + CAT + ".root", "RECREATE")
w_sig = ROOT.RooWorkspace("workspace_sig","workspace_sig")
getattr(w_sig, "import")(sig_model.pdf)
w_sig.Print()
w_sig.Write()
f_out1.Close()
