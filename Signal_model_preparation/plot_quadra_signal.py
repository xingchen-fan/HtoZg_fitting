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
parser = argparse.ArgumentParser(description = "Plot quadra signal model")
parser.add_argument('-c', '--cat', help="category")
parser.add_argument('-con', '--config', help = 'Bkg Configuration', default = '../Config/chi2_config_xgboost_nodrop.json')
args = parser.parse_args()
jfile = open(args.config, 'r')
configs = json.load(jfile)
CAT = args.cat
setting = configs[CAT]
lowx = setting["Range"]

x = ROOT.RooRealVar("CMS_hzg_mass_"+CAT, "CMS_hzg_mass_"+CAT, 115, 130)
x.setBins(60)
MH = ROOT.RooRealVar("MH","MH"       ,125)
year = ['2016', '2016APV', '2017', '2018', '2022', '2022EE', '2023', '2023BPix']
bdt1 = 1.
bdt2 = 0.
if CAT == 'ggf1':
    bdt1 = 1.
    bdt2 = 0.81
elif CAT == 'ggf2':
    bdt1 = 0.81
    bdt2 = 0.64
elif CAT == 'ggf3':
    bdt1 = 0.64
    bdt2 = 0.47
elif CAT == 'ggf4':
    bdt1 = 0.47
    bdt2 = -1
elif CAT == 'vbf1':
    bdt1 = 1.
    bdt2 = 0.489
elif CAT == 'vbf2':
    bdt1 = 0.489
    bdt2 = 0.286
elif CAT == 'vbf3':
    bdt1 = 0.286
    bdt2 = 0.083
elif CAT == 'vbf4':
    bdt1 = 0.083
    bdt2 = -1

hist_TH1_el_ggf = ROOT.TH1F('hist_th1_el_ggf', 'hist_th1_el_ggf', 260, lowx, lowx+65)
hist_TH1_el_vbf = ROOT.TH1F('hist_th1_el_vbf', 'hist_th1_el_vbf', 260, lowx, lowx+65)
hist_TH1_mu_ggf = ROOT.TH1F('hist_th1_mu_ggf', 'hist_th1_mu_ggf', 260, lowx, lowx+65)
hist_TH1_mu_vbf = ROOT.TH1F('hist_th1_mu_vbf', 'hist_th1_mu_vbf', 260, lowx, lowx+65)
hist_TH1 = ROOT.TH1F('hist_th1', 'hist_th1', 260, lowx, lowx+65)
prod1 = ['GGF']
prod2 = ['VBF', 'ZH', 'WH', 'ttH']

if 'ggf' in CAT:
    direct = '/eos/user/r/rzou/SWAN_projects/Classifier/Output_ggF_pinnacles_fix/relpt_peking_run2p3/TreeS/'
    chain1 = ROOT.TChain('TreeS')
    chain2 = ROOT.TChain('TreeS')
    for y in year:
        for p in prod1:
            chain1.Add(direct +p+'_'+y+'_pinnacles_ggf_fixed.root')
        for p in prod2:
            chain2.Add(direct +p+'_'+y+'_pinnacles_ggf_fixed.root')

    for entry in chain1:
        if entry.bdt_score_test > bdt2 and entry.bdt_score_test < bdt1 and entry.met < 90 and entry.weight_corr < 0.5:
            if entry.ll_lepid == 11:
                hist_TH1_el_ggf.Fill(entry.llphoton_refit_m, entry.weight_corr)
                hist_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)
            elif entry.ll_lepid == 13:
                hist_TH1_mu_ggf.Fill(entry.llphoton_refit_m, entry.weight_corr)
                hist_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)
    for entry in chain2:
        if entry.bdt_score_test > bdt2 and entry.bdt_score_test < bdt1 and entry.met < 90 and entry.weight_corr < 0.5:
            if entry.ll_lepid == 11:
                hist_TH1_el_vbf.Fill(entry.llphoton_refit_m, entry.weight_corr)
                hist_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)
            elif entry.ll_lepid == 13:
                hist_TH1_mu_vbf.Fill(entry.llphoton_refit_m, entry.weight_corr)
                hist_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)

elif 'vbf' in CAT:
    direct = '/eos/user/r/rzou/SWAN_projects/Classifier/Output_2JClassic_input_run2p3/'
    chain1 = ROOT.TChain('outtree')
    chain2 = ROOT.TChain('outtree')
    for y in year:
        for p in prod1:
            chain1.Add(direct +p+'_'+y+'_output.root')
        for p in prod2:
            chain2.Add(direct +p+'_'+y+'_output.root')

    for entry in chain1:
        if entry.BDT_score_2j > bdt2 and entry.BDT_score_2j < bdt1 and entry.weight_corr < 0.5:
            if entry.ll_lepid == 11:
                hist_TH1_el_ggf.Fill(entry.llphoton_refit_m, entry.weight_corr)
                hist_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)
            elif entry.ll_lepid == 13:
                hist_TH1_mu_ggf.Fill(entry.llphoton_refit_m, entry.weight_corr)
                hist_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)
    for entry in chain2:
        if entry.BDT_score_2j > bdt2 and entry.BDT_score_2j < bdt1 and entry.weight_corr < 0.5:
            if entry.ll_lepid == 11:
                hist_TH1_el_vbf.Fill(entry.llphoton_refit_m, entry.weight_corr)
                hist_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)
            elif entry.ll_lepid == 13:
                hist_TH1_mu_vbf.Fill(entry.llphoton_refit_m, entry.weight_corr)
                hist_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)

hist_el_ggf = ROOT.RooDataHist('hist_el_ggf', 'hist_el_ggf', x, hist_TH1_el_ggf)
hist_el_vbf = ROOT.RooDataHist('hist_el_vbf', 'hist_el_vbf', x, hist_TH1_el_vbf)
hist_mu_ggf = ROOT.RooDataHist('hist_mu_ggf', 'hist_mu_ggf', x, hist_TH1_mu_ggf)
hist_mu_vbf = ROOT.RooDataHist('hist_mu_vbf', 'hist_mu_vbf', x, hist_TH1_mu_vbf)
hist = ROOT.RooDataHist('hist_ggf', 'hist_ggf', x, hist_TH1)
    
sig_model_el_ggf = DSCB_Class(x, MH, CAT+'_el_ggf',  sigmaL_init = 1.2, sigmaR_init = 1.2, nL_init = 4, nL_bond = 100, nR_init = 8, nR_bond = 100, alphaL_init = 1.5, alphaR_init = 1.5, di_sigma = True)
sig_model_el_vbf = DSCB_Class(x, MH, CAT+'_el_vbf',  sigmaL_init = 1.2, sigmaR_init = 1.2, nL_init = 4, nL_bond = 100, nR_init = 8, nR_bond = 100, alphaL_init = 1.5, alphaR_init = 1.5, di_sigma = True)
sig_model_mu_ggf = DSCB_Class(x, MH, CAT+'_mu_ggf',  sigmaL_init = 1.2, sigmaR_init = 1.2, nL_init = 4, nL_bond = 100, nR_init = 8, nR_bond = 100, alphaL_init = 1.5, alphaR_init = 1.5, di_sigma = True)
sig_model_mu_vbf = DSCB_Class(x, MH, CAT+'_mu_vbf',  sigmaL_init = 1.2, sigmaR_init = 1.2, nL_init = 4, nL_bond = 100, nR_init = 8, nR_bond = 100, alphaL_init = 1.5, alphaR_init = 1.5, di_sigma = True)

sig_model_el_ggf.pdf.fitTo(hist_el_ggf, ROOT.RooFit.AsymptoticError(True), ROOT.RooFit.Save(True),ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0))
sig_model_el_ggf.pdf.fitTo(hist_el_ggf, ROOT.RooFit.AsymptoticError(True), ROOT.RooFit.Save(True),ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0))
sig_model_el_ggf.setStable()
res = sig_model_el_ggf.pdf.fitTo(hist_el_ggf, ROOT.RooFit.AsymptoticError(True), ROOT.RooFit.Save(True),ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0))
#res.Print('V')

sig_model_el_vbf.pdf.fitTo(hist_el_vbf, ROOT.RooFit.AsymptoticError(True), ROOT.RooFit.Save(True),ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0))
sig_model_el_vbf.pdf.fitTo(hist_el_vbf, ROOT.RooFit.AsymptoticError(True), ROOT.RooFit.Save(True),ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0))
sig_model_el_vbf.setStable()
res = sig_model_el_vbf.pdf.fitTo(hist_el_vbf, ROOT.RooFit.AsymptoticError(True), ROOT.RooFit.Save(True),ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0))

sig_model_mu_ggf.pdf.fitTo(hist_mu_ggf, ROOT.RooFit.AsymptoticError(True), ROOT.RooFit.Save(True),ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0))
sig_model_mu_ggf.pdf.fitTo(hist_mu_ggf, ROOT.RooFit.AsymptoticError(True), ROOT.RooFit.Save(True),ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0))
sig_model_mu_ggf.setStable()
res = sig_model_mu_ggf.pdf.fitTo(hist_mu_ggf, ROOT.RooFit.AsymptoticError(True), ROOT.RooFit.Save(True),ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0), ROOT.RooFit.Range('signal'))
#res.Print('V')

sig_model_mu_vbf.pdf.fitTo(hist_mu_vbf, ROOT.RooFit.AsymptoticError(True), ROOT.RooFit.Save(True),ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0))
sig_model_mu_vbf.pdf.fitTo(hist_mu_vbf, ROOT.RooFit.AsymptoticError(True), ROOT.RooFit.Save(True),ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0))
sig_model_mu_vbf.setStable()
res = sig_model_mu_vbf.pdf.fitTo(hist_mu_vbf, ROOT.RooFit.AsymptoticError(True), ROOT.RooFit.Save(True),ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0), ROOT.RooFit.Range('signal'))
#res.Print('V')

#note_el = '#splitline{#splitline{#sigmaL = ' + '%.2f'%sig_model_el.sigmaL.getVal() + '#pm%.2f'%sig_model_el.sigmaL.getError()+', #sigmaR = ' + '%.2f'%sig_model_el.sigmaR.getVal() + '#p\
#m%.2f'%sig_model_el.sigmaR.getError() + '}{nL = %.2f'%sig_model_el.nL.getVal()+', nR = %.2f'%sig_model_el.nR.getVal()+', #mu = %.2f'%MH.getVal() +  '}}{#alphaL = %.2f'%sig_model_el.alphaL.getVal() + ', #alphaR = %.2f'%sig_model_el.alphaR.getVal()+'}'
#plotClass(x, hist_el, sig_model_el.pdf, sig_model_el.pdf, CAT+"_duo_el", "", note=note_el, CMS = "Simulation", fullRatio = False, leftSpace=True, bins = 2)
#note_mu = '#splitline{#splitline{#sigmaL = ' + '%.2f'%sig_model_mu.sigmaL.getVal() + '#pm%.2f'%sig_model_mu.sigmaL.getError()+', #sigmaR = ' + '%.2f'%sig_model_mu.sigmaR.getVal() + '#p\
#m%.2f'%sig_model_mu.sigmaR.getError() + '}{nL = %.2f'%sig_model_mu.nL.getVal()+', nR = %.2f'%sig_model_mu.nR.getVal()+', #mu = %.2f'%MH.getVal() +  '}}{#alphaL = %.2f'%sig_model_mu.alphaL.getVal() + ', #alphaR = %.2f'%sig_model_mu.alphaR.getVal()+'}'
#plotClass(x, hist_mu, sig_model_mu.pdf, sig_model_mu.pdf, CAT+"_duo_mu", "", note=note_mu, CMS = "Simulation", fullRatio = False, leftSpace=True, bins = 2)
tot = hist_el_ggf.sumEntries() + hist_el_vbf.sumEntries() + hist_mu_ggf.sumEntries() + hist_mu_vbf.sumEntries()
c_el_ggf = ROOT.RooRealVar('c_el_ggf', 'c_el_ggf', hist_el_ggf.sumEntries()/tot)
c_el_vbf = ROOT.RooRealVar('c_el_vbf', 'c_el_vbf', hist_el_vbf.sumEntries()/tot)
c_mu_ggf = ROOT.RooRealVar('c_mu_ggf', 'c_mu_ggf', hist_mu_ggf.sumEntries()/tot)
#c_mu_vbf = ROOT.RooRealVar('c_mu_vbf', 'c_mu_vbf', hist_mu_vbf.sumEntries())

quadra_model = ROOT.RooAddPdf('quadra_sig_model_'+CAT, 'quadra_sig_model_'+CAT, ROOT.RooArgList(sig_model_el_ggf.pdf, sig_model_el_vbf.pdf, sig_model_mu_ggf.pdf, sig_model_mu_vbf.pdf), ROOT.RooArgList(c_el_ggf, c_el_vbf, c_mu_ggf))
quadra_model.fixCoefNormalization(ROOT.RooArgSet(x))

chi2 = ROOT.RooChi2Var("chi2", "chi2", quadra_model, hist, ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))
PV = ROOT.TMath.Prob(chi2.getVal(), 60)
sigma_eff_pair = getEffSigma(x, quadra_model)
sigma_eff = (sigma_eff_pair[1] - sigma_eff_pair[0])/2

plotClass(x, hist, quadra_model, quadra_model, CAT+"_quadra", "", CMS = "Simulation", fullRatio = True, leftSpace=True, bins = 2, note="#splitline{#Chi^{2} = %.2f"%(chi2.getVal()) + " P-value = %.3f}"%PV + "{#sigma_{eff} = %.2f}"%sigma_eff)

