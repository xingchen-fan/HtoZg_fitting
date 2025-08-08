#!/usr/bin/env python3
import os
import sys
import ROOT
import argparse
import json
import math
sys.path.append(os.path.abspath("../../Utilities/"))
sys.path.append(os.path.abspath("../../CMS_plotter/"))
import CMS_lumi, tdrstyle
from sig_functions_class import *
from Xc_Minimizer import *
from plot_utility import *
from sample_reader import *

ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--cat', help = 'Category')
parser.add_argument('-y', '--year', help = 'Year')
parser.add_argument('-p', '--prod', help = 'Production mode')
parser.add_argument('-plot', '--plot', help = 'Plot dir', default="")
args = parser.parse_args()
CAT = args.cat
x = ROOT.RooRealVar("x", "mllg", 115, 130)
MH_el =  ROOT.RooRealVar("MH_el", "MH_el", 125)
MH_mu = ROOT.RooRealVar("MH_mu", "MH_mu", 125)
x.setBins(int(4* (x.getMax() - x.getMin())))
if "ggf" in CAT and args.prod == "ggf":
    sig_sample = readRuiROOTggFSignalggF(x, '/eos/user/r/rzou/SWAN_projects/Classifier/Output_ggF_pinnacles_fix/relpt_peking_run2p3/TreeS/', args.year, 0.81,0.64,0.47)
    if CAT == 'ggf1':
        el_sig_hist = sig_sample.ggf1El
        mu_sig_hist = sig_sample.ggf1Mu
    elif CAT == 'ggf2':
        el_sig_hist = sig_sample.ggf2El
        mu_sig_hist = sig_sample.ggf2Mu
    elif CAT == 'ggf3':
        el_sig_hist = sig_sample.ggf3El
        mu_sig_hist = sig_sample.ggf3Mu
    elif CAT == 'ggf4':
        el_sig_hist = sig_sample.ggf4El
        mu_sig_hist = sig_sample.ggf4Mu

elif "ggf" in CAT and args.prod == "vbf":
    sig_sample = readRuiROOTggFSignalVBF(x, '/eos/user/r/rzou/SWAN_projects/Classifier/Output_ggF_pinnacles_fix/relpt_peking_run2p3/TreeS/', args.year, 0.81,0.64,0.47)
    if CAT == 'ggf1':
        el_sig_hist = sig_sample.ggf1El
        mu_sig_hist = sig_sample.ggf1Mu
    elif CAT == 'ggf2':
        el_sig_hist = sig_sample.ggf2El
        mu_sig_hist = sig_sample.ggf2Mu
    elif CAT == 'ggf3':
        el_sig_hist = sig_sample.ggf3El
        mu_sig_hist = sig_sample.ggf3Mu
    elif CAT == 'ggf4':
        el_sig_hist = sig_sample.ggf4El
        mu_sig_hist = sig_sample.ggf4Mu

elif "vbf" in CAT and args.prod == "ggf":
    sig_sample = readRuiROOTVBFSignalggF(x, '/eos/user/r/rzou/SWAN_projects/Classifier/Output_2JClassic_input_run2p3/', args.year, 0.489,0.286 , 0.083)
    if CAT == 'vbf1':
        el_sig_hist = sig_sample.vbf1El
        mu_sig_hist = sig_sample.vbf1Mu
    elif CAT == 'vbf2':
        el_sig_hist = sig_sample.vbf2El
        mu_sig_hist = sig_sample.vbf2Mu
    elif CAT == 'vbf3':
        el_sig_hist = sig_sample.vbf3El
        mu_sig_hist = sig_sample.vbf3Mu
    elif CAT == 'vbf4':
        el_sig_hist = sig_sample.vbf4El
        mu_sig_hist = sig_sample.vbf4Mu

elif "vbf" in CAT and args.prod == "vbf":
    sig_sample = readRuiROOTVBFSignalVBF(x, '/eos/user/r/rzou/SWAN_projects/Classifier/Output_2JClassic_input_run2p3/', args.year, 0.489,0.286 , 0.083)
    if CAT == 'vbf1':
        el_sig_hist = sig_sample.vbf1El
        mu_sig_hist = sig_sample.vbf1Mu
    elif CAT == 'vbf2':
        el_sig_hist = sig_sample.vbf2El
        mu_sig_hist = sig_sample.vbf2Mu
    elif CAT == 'vbf3':
        el_sig_hist = sig_sample.vbf3El
        mu_sig_hist = sig_sample.vbf3Mu
    elif CAT == 'vbf4':
        el_sig_hist = sig_sample.vbf4El
        mu_sig_hist = sig_sample.vbf4Mu

DISIGMA = True
sig_model_el = DSCB_Class(x, MH_el,  CAT+"_"+args.year+"_el_"+args.prod, sigmaL_init = 1.2, sigmaR_init = 1.2, nL_init = 4, nL_bond = 100, nR_init = 8, nR_bond = 100, alphaL_init = 1.5, alphaR_init = 1.5, di_sigma = DISIGMA)
sig_model_mu = DSCB_Class(x, MH_mu,  CAT+"_"+args.year+"_mu_"+args.prod, sigmaL_init = 1.2, sigmaR_init = 1.2, nL_init = 5.92, nL_bond = 100, nR_init = 10.7, nR_bond = 100, alphaL_init = 0.69, alphaR_init = 1.33, di_sigma = DISIGMA)
#print(el_sig_hist.sumEntries())
#nll_el = ROOT.RooNLLVar("nll_"+ CAT+"_"+args.year+"_el", "nll_"+ CAT+"_"+args.year+"_el", sig_model_el.pdf, el_sig_hist, ROOT.RooFit.AsymptoticError(True))
#nll_mu = ROOT.RooNLLVar("nll_"+ CAT+"_"+args.year+"_mu", "nll_"+ CAT+"_"+args.year+"_mu", sig_model_mu.pdf, mu_sig_hist, ROOT.RooFit.AsymptoticError(True))

#Minimizer_NLL(nll_el, -1, 100, False, 0)
#Minimizer_NLL(nll_el, -1, 1, False, 0)
sig_model_el.pdf.fitTo(el_sig_hist, ROOT.RooFit.AsymptoticError(True), ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0))
sig_model_el.pdf.fitTo(el_sig_hist, ROOT.RooFit.AsymptoticError(True), ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0))
#sig_model_el.pdf.chi2FitTo(el_sig_hist, ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2), ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0))
sig_model_el.setStable()
#Minimizer_NLL(nll_el, -1, 0.1, False, 0)

res_el = sig_model_el.pdf.fitTo(el_sig_hist, ROOT.RooFit.AsymptoticError(True), ROOT.RooFit.Save(True),ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0))
#res_el = sig_model_el.pdf.chi2FitTo(el_sig_hist, ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2), ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0), ROOT.RooFit.Save(True))
res_el.Print('v')
print('nexp_'+CAT+'_el = %.2f'%el_sig_hist.sumEntries())
el_pair = getEffSigma(x, sig_model_el.pdf)
sigma_eff_el = (el_pair[1] - el_pair[0])/2

#Minimizer_NLL(nll_mu, -1, 100, False, 0)
#Minimizer_NLL(nll_mu, -1, 1, False, 0)
sig_model_mu.pdf.fitTo(mu_sig_hist, ROOT.RooFit.AsymptoticError(True), ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0))
sig_model_mu.pdf.fitTo(mu_sig_hist, ROOT.RooFit.AsymptoticError(True), ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0))
#sig_model_mu.pdf.chi2FitTo(mu_sig_hist, ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2), ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0))
sig_model_mu.setStable()
#Minimizer_NLL(nll_mu, -1, 0.1, False, 0)
res_mu = sig_model_mu.pdf.fitTo(mu_sig_hist, ROOT.RooFit.AsymptoticError(True), ROOT.RooFit.Save(True), ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0))
#res_mu = sig_model_mu.pdf.chi2FitTo(mu_sig_hist, ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2), ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0), ROOT.RooFit.Save(True))
res_mu.Print('v')
print('nexp_'+CAT+'_mu = %.2f'%mu_sig_hist.sumEntries())
mu_pair = getEffSigma(x, sig_model_mu.pdf)
sigma_eff_mu = (mu_pair[1] - mu_pair[0])/2

e1_el = sig_model_el.sigmaL.getError()
e2_el = sig_model_el.sigmaR.getError()
e1_mu = sig_model_mu.sigmaL.getError()
e2_mu = sig_model_mu.sigmaR.getError()

if DISIGMA:
    print("sigma_eff_err_el = %.2f"%math.sqrt((e1_el*e1_el + e2_el*e2_el)/2))
    print("sigma_eff_err_el = %.2f"%math.sqrt((e1_mu*e1_mu + e2_mu*e2_mu)/2))
    note_el = '#splitline{#sigma_{eff} = %.2f'%sigma_eff_el +'#pm%.2f'%math.sqrt((e1_el*e1_el + e2_el*e2_el)/2) + '}{#splitline{#splitline{#sigmaL = ' + '%.2f'%sig_model_el.sigmaL.getVal() + '#pm%.2f'%sig_model_el.sigmaL.getError()+', #sigmaR = ' + '%.2f'%sig_model_el.sigmaR.getVal() + '#pm%.2f'%sig_model_el.sigmaR.getError() + '}{nL = %.2f'%sig_model_el.nL.getVal()+', nR = %.2f'%sig_model_el.nR.getVal()+', #mu = %.2f'%(sig_model_el.dMH.getVal()+125) +  '}}{#alphaL = %.2f'%sig_model_el.alphaL.getVal() + ', #alphaR = %.2f'%sig_model_el.alphaR.getVal()+'}}'
    note_mu = '#splitline{#sigma_{eff} = %.2f'%sigma_eff_mu +'#pm%.2f'%math.sqrt((e1_mu*e1_mu + e2_mu*e2_mu)/2) + '}{#splitline{#splitline{#sigmaL = ' + '%.2f'%sig_model_mu.sigmaL.getVal() + '#pm%.2f'%sig_model_mu.sigmaL.getError()+ ', #sigmaR = ' + '%.2f'%sig_model_mu.sigmaR.getVal() + '#pm%.2f'%sig_model_mu.sigmaR.getError() + '}{nL = %.2f'%sig_model_mu.nL.getVal()+', nR = %.2f'%sig_model_mu.nR.getVal()+', #mu = %.2f'%(sig_model_mu.dMH.getVal()+125) +  '}}{#alphaL = %.2f'%sig_model_mu.alphaL.getVal() + ', #alphaR = %.2f'%sig_model_mu.alphaR.getVal()+'}}'

else:
    print("sigma_eff_err_el = %.2f"%e1_el)
    print("sigma_eff_err_el = %.2f"%e1_mu)
    note_el = '#splitline{#splitline{#sigma_{eff} = %.2f'%sigma_eff_el + '#pm%.2f'%e1_el + ' #sigma = ' + '%.2f'%sig_model_el.sigmaL.getVal() + '#pm%.2f'%sig_model_el.sigmaL.getError()+'}{nL = %.2f'%sig_model_el.nL.getVal()+', nR = %.2f'%sig_model_el.nR.getVal()+', #mu = %.2f'%(sig_model_el.dMH.getVal()+125) +  '}}{#alphaL = %.2f'%sig_model_el.alphaL.getVal() + ', #alphaR = %.2f'%sig_model_el.alphaR.getVal()+'}'
    note_mu = '#splitline{#splitline{#sigma_{eff} = %.2f'%sigma_eff_mu + '#pm%.2f'%e1_mu +' #sigma = ' + '%.2f'%sig_model_mu.sigmaL.getVal() + '#pm%.2f'%sig_model_mu.sigmaL.getError()+ '}{nL = %.2f'%sig_model_mu.nL.getVal()+', nR = %.2f'%sig_model_mu.nR.getVal()+', #mu = %.2f'%(sig_model_mu.dMH.getVal()+125) +  '}}{#alphaL = %.2f'%sig_model_mu.alphaL.getVal() + ', #alphaR = %.2f'%sig_model_mu.alphaR.getVal()+'}'
plotClass(x, el_sig_hist, sig_model_el.pdf, sig_model_el.pdf, CAT+"_"+args.year+"_el_"+args.prod, args.plot+"/", note=note_el, CMS = "Simulation", fullRatio = True, leftSpace=True, bins = 2)
plotClass(x, mu_sig_hist, sig_model_mu.pdf, sig_model_mu.pdf, CAT+"_"+args.year+"_mu_"+args.prod, args.plot+"/", note=note_mu, CMS = "Simulation", fullRatio = True, leftSpace=True, bins = 2)

