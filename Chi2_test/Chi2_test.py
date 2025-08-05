#!/usr/bin/env python3
import os
import sys
import ROOT
import argparse
import json
sys.path.append(os.path.abspath("../Utilities/"))
sys.path.append(os.path.abspath("../CMS_plotter/"))
import CMS_lumi, tdrstyle
from bkg_functions_class import *
from Xc_Minimizer import *
from plot_utility import *
from sample_reader import *
from profile_class import *
#ROOT.gInterpreter.AddIncludePath('../Utilities/HZGRooPdfs.h')
#ROOT.gSystem.Load('../Utilities/HZGRooPdfs_cxx.so')
#ROOT.gInterpreter.AddIncludePath('../Utilities/AsymGenGaussian.h')
#ROOT.gSystem.Load('../Utilities/AsymGenGaussian_cxx.so')
ROOT.gInterpreter.AddIncludePath('../Utilities/RooGaussStepBernstein.h')
ROOT.gSystem.Load('../Utilities/RooGaussStepBernstein_cxx.so')

ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--cat', help = 'Category')
parser.add_argument('-con', '--config', help = 'Configuration')
args = parser.parse_args()
CAT = args.cat
jfile = open(args.config, 'r')
configs = json.load(jfile)
setting = configs[CAT]

# Define vars
lowx = setting["Range"]
x = ROOT.RooRealVar("x", "mllg", lowx, lowx+65.)
y = ROOT.RooRealVar("y", "photon pT", 15., 1000.)
w = ROOT.RooRealVar("w", "w", -40., 40.)
bdt = ROOT.RooRealVar("bdt", "bdt", -1, 1)
year = ROOT.RooRealVar("year", "year", 2015, 21000)
w_year = ROOT.RooRealVar("w_year", "w_year", 0, 100)
lep = ROOT.RooRealVar("lep", "lep", 0, 1) #0 = electron, 1 = muon
ph_eta = ROOT.RooRealVar("ph_eta", "ph_eta", -3, 3)
nlep = ROOT.RooRealVar("nlep", "nlep", 0, 10)
njet = ROOT.RooRealVar("njet", "njet", 0, 10)

x.setBins(260)
x.setRange('left', lowx, 120)
x.setRange('right', 130, lowx+65)
x.setRange('full', lowx, lowx+65)

if CAT == 'ggf1':
    NAMECAT = 'ggF1'
elif CAT == 'ggf2':
    NAMECAT = 'ggF2'
elif CAT == 'ggf3':
    NAMECAT = 'ggF3'
elif CAT == 'ggf4':
    NAMECAT = 'ggF4'
else:
    NAMECAT = 'WRONG'
# Read samples
DAT = False
if DAT:
    read_data = ROOT.RooDataSet.read('../Data/data_'+NAMECAT+'_pinnacles_fix.dat', ROOT.RooArgList(x, y, bdt, w, w_year, year, lep, ph_eta, nlep, njet))
    data = ROOT.RooDataSet('data', 'data', read_data, ROOT.RooArgList(x, y, bdt, w, w_year, year, lep, ph_eta, nlep, njet),'')
    hist_data = ROOT.RooDataHist('hist_data','hist_data', x, data)
else:
    if 'ggf' in CAT:
        read_data = readRuiROOTggFdata(x, '/eos/user/r/rzou/SWAN_projects/Classifier/Output_ggF_xgb_fixedjet/', 0.82, 0.64,0.46)
        if CAT == 'ggf1':
            hist_data = read_data.ggf1
        elif CAT == 'ggf2':
            hist_data = read_data.ggf2
        elif CAT == 'ggf3':
            hist_data = read_data.ggf3
        elif CAT == 'ggf4':
            hist_data = read_data.ggf4
    elif 'vbf' in CAT:
        read_data = readRuiROOTVBFdata(x, '/eos/user/r/rzou/SWAN_projects/Classifier/Output_VBF_xgb_neweight/', 0.95, 0.83, 0.58)
        if CAT == 'vbf1':
            hist_data = read_data.vbf1
        elif CAT == 'vbf2':
            hist_data = read_data.vbf2
        elif CAT == 'vbf3':
            hist_data = read_data.vbf3
        elif CAT == 'vbf4':
            hist_data = read_data.vbf4
print('N data = ', hist_data.sumEntries())
# Bkg funcs
mu_gauss = ROOT.RooRealVar("mu_gauss","always 0"       ,0.)

#AGG_model = AGGClass(x, CAT, kappa_init = -1.27, alpha_init = 14, zeta_init = 105, x_low = lowx, x_high = lowx+65)
profile = profileClass(x, mu_gauss, CAT, args.config)
bkg_list = [profile.bern2_model, profile.bern3_model, profile.bern4_model, profile.bern5_model, profile.pow1_model, profile.pow2_model, profile.pow3_model, profile.exp1_model, profile.exp2_model, profile.exp3_model, profile.lau2_model, profile.lau3_model, profile.lau4_model, profile.modg_model]

# Sideband Chi^2 fit
cuthistogram = hist_data.reduce(ROOT.RooFit.CutRange('left,right'))
th1_cuthist = cuthistogram.createHistogram("cuthist_ks", x, ROOT.RooFit.Binning(260))
n_bins = 220
NLL = True
pass_list = []
for func in bkg_list:
    if NLL:
        chi2 = ROOT.RooNLLVar('chi2_'+func.SBpdf.GetName(), 'chi2_'+func.SBpdf.GetName(), func.SBpdf, cuthistogram)
    else:
        chi2 = ROOT.RooChi2Var('chi2_'+func.SBpdf.GetName(), 'chi2_'+func.SBpdf.GetName(), func.SBpdf, cuthistogram)
    Minimizer_Chi2(chi2, -1, 10, True, 0)
    res=Minimizer_Chi2(chi2, -1, 0.1, True, 0) 
    func.checkBond()
    res.Print('v')
    chi2_ = ROOT.RooChi2Var('chi2_val_'+func.SBpdf.GetName(), 'chi2_val_'+func.SBpdf.GetName(), func.SBpdf, cuthistogram)
    chi2_val = chi2_.getVal()
    dof = n_bins - res.floatParsFinal().getSize()
    chi2_PV = ROOT.Math.chisquared_cdf_c(chi2_.getVal(), dof)
    ks_model_hist = func.SBpdf.generateBinned(x, cuthistogram.sumEntries(), True).createHistogram("hist_"+func.SBpdf.GetName(), x, ROOT.RooFit.Binning(260))
    ks_PV = th1_cuthist.KolmogorovTest(ks_model_hist)
    plotClass(x, hist_data, func.pdf, func.SBpdf, title=func.pdf.GetName(), output_dir="plots/", sideBand = True, fitRange = 'left,right', note='#splitline{#Chi^{2}/dof = ' + '%.2f'%chi2_val + '/%i'%dof +'}{#splitline{#Chi^{2} P-value = ' + '%.2f'%chi2_PV + '}{KS P-value = ' + '%.2f'%ks_PV + '}}', fullRatio = True)
    if chi2_PV>0.1 and ks_PV>0.1 and res.status() == 0:
        pass_list.append(func)
print('List has ', len(pass_list), ' models')
if len(pass_list)>0:
    multiPlotClass(x, hist_data, pass_list, title="Chi2_single_"+CAT, output_dir="plots/",sideBand = True, fitRange = 'left,right',best_index = 0,CMS = "Preliminary", fullRatio = ("vbf" in CAT))
    profile.write_config_file(cuthistogram, "All")
