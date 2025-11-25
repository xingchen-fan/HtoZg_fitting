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
parser.add_argument('-c', '--cat', help="Category")
parser.add_argument('-con', '--config', help="Configuration")
args = parser.parse_args()
CAT = args.cat

x = ROOT.RooRealVar("CMS_hzg_mass_"+CAT, "CMS_hzg_mass_"+CAT, 115, 130)
lowx = x.getMin()
highx = x.getMax()
nbins = int(4* (highx - lowx))
x.setBins(nbins)
nsig_el = 0
nsig_mu = 0
#x.setRange('signal', 115, 130)
MH = ROOT.RooRealVar("MH","MH"       ,125)

if 'ggf' in CAT:
    direct = '/eos/project/h/htozg-dy-privatemc/rzou/bdt/BDT_output_redwood/Output_ggF_rui_redwood_v1_ext_val/'
    reader = readRuiROOTggFSignal(x, direct, 0.94, 0.83, 0.57)
elif 'vbf' in CAT:
    direct = '/eos/project/h/htozg-dy-privatemc/rzou/bdt/BDT_output_redwood/Output_VBF_rui_redwood_v1_ext_val/'
    reader = readRuiROOTVBFSignal(x, direct, 0.91, 0.81, 0.48)
    
hist_el = getattr(reader, f"{CAT}El")
hist_mu = getattr(reader, f"{CAT}Mu")
nsig_el = getattr(reader,f"n{CAT}_el")
nsig_mu = getattr(reader,f"n{CAT}_mu")
hist = getattr(reader, CAT)
    
model_params = {}
with open(args.config, "r") as config_file:
    dscb_config = json.load(config_file)
for lep in ["el", "mu"]:
    model_name = f"Htozg_{lep}_cat_{CAT}_nominal"
    model_params = {}
    sig_model = DSCB_Class(x, MH, CAT+'_'+lep,  sigmaL_init = 1.2, sigmaR_init = 1.2, nL_init = 4, nL_bond = 100, nR_init = 8, nR_bond = 100, alphaL_init = 1.5, alphaR_init = 1.5, di_sigma = True)
    if lep=="el":
        fit_hist = hist_el
    elif lep=="mu":
        fit_hist = hist_mu
    nll = ROOT.RooNLLVar("nll_"+ CAT+"_"+lep, "nll_"+ CAT+"_"+lep, sig_model.pdf, fit_hist, ROOT.RooFit.AsymptoticError(True))
    #sig_model.pdf.fitTo(fit_hist, ROOT.RooFit.AsymptoticError(True), ROOT.RooFit.Save(True),ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0))
    #sig_model.pdf.fitTo(fit_hist, ROOT.RooFit.AsymptoticError(True), ROOT.RooFit.Save(True),ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0))
    Minimizer_NLL(nll, -1, 100, False, 0)
    Minimizer_NLL(nll, -1, 1, True, 0)
    Minimizer_NLL(nll, -1, 0.001, True, 0)
    sig_model.setStable()
    res = Minimizer_NLL(nll, -1, 0.001, True, 0)
    #res = sig_model.pdf.fitTo(fit_hist, ROOT.RooFit.AsymptoticError(True), ROOT.RooFit.Save(True),ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0))
    res.Print('V')
    note = '#splitline{#splitline{#sigmaL = ' + '%.2f'%sig_model.sigmaL.getVal() + '#pm%.2f'%sig_model.sigmaL.getError()+', #sigmaR = ' + '%.2f'%sig_model.sigmaR.getVal() + '#p\
m%.2f'%sig_model.sigmaR.getError() + '}{nL = %.2f'%sig_model.nL.getVal()+', nR = %.2f'%sig_model.nR.getVal()+', #mu = %.2f'%(sig_model.dMH.getVal()+MH.getVal()) +  '}}{#alphaL = %.2f'%sig_model.alphaL.getVal() + ', #alphaR = %.2f'%sig_model.alphaR.getVal()+'}'
    plotClass(x, fit_hist, sig_model.pdf, sig_model.pdf, CAT+"_duo_"+lep, "plots/", note=note, CMS = "Simulation", fullRatio = True, leftSpace=True, bins = 2)
    model_params["disigma"] = sig_model.disigma
    model_params["nexp"] = nsig_el if lep == 'el' else nsig_mu
    for key in ["sigmaL", "sigmaR", "nL", "nR", "alphaL", "alphaR"]:
        model_params[key] = getattr(sig_model,key).getVal()
    model_params["dMH"] = sig_model.dMH.getVal()
    dscb_config[model_name] = model_params
with open(args.config, "w") as config_file:
      config_file.write(json.dumps(dscb_config,indent=2))

"""
c_el = ROOT.RooRealVar('c_el', 'c_el', hist_el.sumEntries()/(hist_el.sumEntries() + hist_mu.sumEntries()))
duo_model = ROOT.RooAddPdf('duo_sig_model_'+CAT, 'duo_sig_model_'+CAT, sig_model_el.pdf, sig_model_mu.pdf, c_el)
duo_model.fixCoefNormalization(ROOT.RooArgSet(x))

chi2 = ROOT.RooChi2Var("chi2", "chi2", duo_model, hist, ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))
PV = ROOT.TMath.Prob(chi2.getVal(), 60)
sigma_eff_pair = getEffSigma(x, duo_model)
sigma_eff = (sigma_eff_pair[1] - sigma_eff_pair[0])/2

plotClass(x, hist, duo_model, duo_model, CAT+"_duo", "", CMS = "Simulation", fullRatio = True, leftSpace=True, bins = 2, note="#splitline{#Chi^{2} = %.2f"%(chi2.getVal()) + " P-value = %.3f}"%PV + "{#sigma_{eff} = %.2f}"%sigma_eff)
"""
