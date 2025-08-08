#!/usr/bin/env python3
import os
import sys
import ROOT
import argparse
import json
import math
sys.path.append(os.path.abspath("../Utilities/"))
sys.path.append(os.path.abspath("../CMS_plotter/"))
import CMS_lumi, tdrstyle
from sig_functions_class import *
from Xc_Minimizer import *
from plot_utility import *
from sample_reader import *

parser = argparse.ArgumentParser(description = "Plot combined signal and single signal")
parser.add_argument('-c', '--cat', help="category")
parser.add_argument('-con', '--config', help = 'Sig Configuration')
args = parser.parse_args()
args = parser.parse_args()
jfile = open(args.config, 'r')
configs = json.load(jfile)
CAT = args.cat
x = ROOT.RooRealVar("CMS_hzg_mass_"+CAT, "CMS_hzg_mass_"+CAT, 115, 130)
x.setBins(60)
x.setRange('signal', 115, 130)
MH = ROOT.RooRealVar("MH","MH"       ,125)
if 'ggf' in CAT:
    direct = '/eos/user/r/rzou/SWAN_projects/Classifier/Output_ggF_rui_commonparam/'
    reader = readRuiROOTggFSignal(x, direct, 0.91,0.82,0.61)
elif 'vbf' in CAT:
    direct = '/eos/user/r/rzou/SWAN_projects/Classifier/Output_ggF_rui_commonparam/'
    reader = readRuiROOTVBFSignal(x, direct, 0.95, 0.91, 0.76)
    
hist = getattr(reader, f"{CAT}")
print('Nsig = ', getattr(reader, f"n{CAT}_el") + getattr(reader, f"n{CAT}_mu"))

SINGLE = False
if SINGLE:
    sig_model = DSCB_Class(x, MH, CAT,  sigmaL_init = 1.2, sigmaR_init = 1.2, nL_init = 4, nL_bond = 100, nR_init = 8, nR_bond = 100, alphaL_init = 1.5, alphaR_init = 1.5, di_sigma = True)
    sig_model.pdf.fitTo(hist, ROOT.RooFit.AsymptoticError(True), ROOT.RooFit.Save(True),ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0))
    sig_model.setStable()
    res = sig_model.pdf.fitTo(hist, ROOT.RooFit.AsymptoticError(True), ROOT.RooFit.Save(True),ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0))
    res.Print("V")
else:
    sig_model = combineSignal(x, MH, CAT, args.config)
    sig_model.pdf.fixCoefNormalization(ROOT.RooArgSet(x))
    print("n tot = ", sig_model.ntot)
"""    
for i in range(60):
    hist.get(i)
    hist.set(hist.weight(), math.sqrt(hist.weight())/32)
"""
chi2 = ROOT.RooChi2Var("chi2", "chi2", sig_model.pdf, hist, ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))
PV = ROOT.TMath.Prob(chi2.getVal(), 60)
sigma_eff_pair = getEffSigma(x, sig_model.pdf)
sigma_eff = (sigma_eff_pair[1] - sigma_eff_pair[0])/2

if SINGLE:
    plotClass(x, hist, sig_model.pdf, sig_model.pdf, CAT+"_single_coarse", "", CMS = "Simulation", fullRatio = True, leftSpace=True, bins = 2, note="#splitline{#Chi^{2} = %.2f"%(chi2.getVal()) + " P-value = %.3f}"%PV + "{#sigma_{eff} = %.2f}"%sigma_eff)
else:
    plotClass(x, hist, sig_model.pdf, sig_model.pdf, CAT+"_combine", "", CMS = "Simulation", fullRatio = True, leftSpace=True, bins = 2, note="#splitline{#Chi^{2} = %.2f"%(chi2.getVal()) + " P-value = %.3f}"%PV + "{#sigma_{eff} = %.2f}"%sigma_eff)
