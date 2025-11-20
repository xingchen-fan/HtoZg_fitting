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
parser.add_argument('-conB', '--configB', help = 'Bkg Configuration', default = '../Config/config_xgboost_0924.json')
parser.add_argument('-conS', '--configS', help = 'Sig Configuration', default = '../Config/config_DSCB_double_flat.json')

args = parser.parse_args()
jfile = open(args.configB, 'r')
jfile_s = open(args.configS, 'r')
configs_s = json.load(jfile_s)
configs = json.load(jfile)
CAT = args.cat
setting = configs[CAT]
lowx = setting["Range"][0]
highx = setting["Range"][1]
nbins = int(setting["Bins"])
x = ROOT.RooRealVar("CMS_hzg_mass_"+CAT, "CMS_hzg_mass_"+CAT, lowx, highx)
x.setBins(nbins)
MH = ROOT.RooRealVar("MH","MH"       ,125)#, 120., 130.)
sig_model_el = DSCB_Class(x, MH, CAT+'_el', di_sigma = True)
sig_model_el.assignVal(args.configS, cat=CAT, lep="el")
sig_model_mu = DSCB_Class(x, MH, CAT+'_mu', di_sigma = True)
sig_model_mu.assignVal(args.configS, cat=CAT, lep="mu")
c_el = ROOT.RooRealVar('c_el_'+CAT, 'c_el_'+CAT, sig_model_el.nsig/(sig_model_el.nsig + sig_model_mu.nsig)) 
duo_model = ROOT.RooAddPdf('duo_sig_model_'+CAT, 'duo_sig_model_'+CAT, sig_model_el.pdf, sig_model_mu.pdf, c_el)
#sig_model = combineSignalLep(x, MH, CAT, args.configS, args.lep)
#print('n tot = ',sig_model.ntot)

f_out1 = ROOT.TFile("workspaces/workspace_duo_sig_" + CAT + ".root", "RECREATE")
w_sig = ROOT.RooWorkspace("workspace_sig","workspace_sig")
getattr(w_sig, "import")(duo_model)
w_sig.Print()
w_sig.Write()
f_out1.Close()
