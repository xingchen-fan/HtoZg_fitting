#!/usr/bin/env python3
import time
import os
import sys
import ROOT
import json
import argparse
sys.path.append(os.path.abspath("../Utilities/"))
sys.path.append(os.path.abspath("../CMS_plotter/"))
import CMS_lumi, tdrstyle
from bkg_functions_fit import *
from bkg_functions_class import *
from Xc_Minimizer import *
from plot_utility import *
from sample_reader import *
ROOT.gInterpreter.AddIncludePath('../Utilities/HZGRooPdfs.h')
ROOT.gSystem.Load('../Utilities/HZGRooPdfs_cxx.so')

ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)
parser = argparse.ArgumentParser()
parser.add_argument('-hi', '--high', help = 'Higher power')
parser.add_argument('-lo', '--low', help = 'Lower power')
args = parser.parse_args()
if int(args.high) <= int(args.low):
    raise Exception('Higher power must be larger than lower power!')
tic = time.perf_counter()
# Define vars
lowx = 105
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

# Read samples
read_data = ROOT.RooDataSet.read('data_ggF3_pinnacles_fix.dat', ROOT.RooArgList(x, y, bdt, w, w_year, year, lep, ph_eta, nlep, njet))
data = ROOT.RooDataSet('data', 'data', read_data, ROOT.RooArgList(x, y, bdt, w, w_year, year, lep, ph_eta, nlep, njet),'')
hist_data = ROOT.RooDataHist('hist_data','hist_data', x, data)

# Bkg funcs
CAT = 'ggf3'
mu_gauss = ROOT.RooRealVar("mu_gauss","always 0"       ,0.)

list_ = []
chi2_list = []
pv_list = []

ONESHOT = True
if ONESHOT:
    lau2_model = Lau2Class(x, mu_gauss, CAT, sigma_init = 3., step_init = 101, p1 = int(args.low), p2 = int(args.high), f_init = 0.8, xmax = lowx+65., const_f1 = True)
    cuthistogram = hist_data.reduce(ROOT.RooFit.CutRange('left,right'))
    n_bins = 220
    chi2 = ROOT.RooChi2Var('chi2_'+lau2_model.SBpdf.GetName(), 'chi2_'+lau2_model.SBpdf.GetName(), lau2_model.SBpdf, cuthistogram)
    Minimizer_Chi2(chi2, -1, 100, False, 0)
    res=Minimizer_Chi2(chi2, -1, 0.1, True, 0)
    lau2_model.checkBond()
    res.Print('v')
    chi2_val = chi2.getVal()
    chi2_pV= ROOT.Math.chisquared_cdf_c(chi2.getVal(), n_bins - res.floatParsFinal().getSize())
    print ('Chi2 = ', chi2_val)
    print ('P-value = ', chi2_pV)
    sys.exit()
for i in range(int(args.low), int(args.high) + 1):
    for j in range(i+1, int(args.high) + 1):
        pow_list = []
        pow_list.append(i)
        pow_list.append(j)
        lau2_model = Lau2Class(x, mu_gauss, CAT, sigma_init = 3., step_init = 101, p1 = i, p2 = j, f_init = 0.8, xmax = lowx+65., const_f1 = True)
        cuthistogram = hist_data.reduce(ROOT.RooFit.CutRange('left,right'))
        n_bins = 220
        chi2 = ROOT.RooChi2Var('chi2_'+lau2_model.SBpdf.GetName(), 'chi2_'+lau2_model.SBpdf.GetName(), lau2_model.SBpdf, cuthistogram)
        Minimizer_Chi2(chi2, -1, 100, False, 0)
        #Minimizer_Chi2(chi2, -1, 0.1, True, 0)
        res=Minimizer_Chi2(chi2, -1, 0.1, True, 0) 
        lau2_model.checkBond()
        res.Print('v')
        chi2_val = chi2.getVal()
        chi2_pV = ROOT.Math.chisquared_cdf_c(chi2.getVal(), n_bins - res.floatParsFinal().getSize())
        chi2_list.append(chi2_val)
        pv_list.append(chi2_pV)
        list_.append(pow_list)
        print ('Finish', pow_list)
min_chi2 = min(chi2_list)
index_ = chi2_list.index(min_chi2)
print ('Power combinations: ', list_)
print ('Chi^2 list: ', chi2_list)
print ('Best lau2 ', list_[index_], ', Chi^2 = ', chi2_list[index_], ', P-value = ', pv_list[index_])
toc = time.perf_counter()
print ('Time spent:', '%.2f'%((toc - tic)/60), 'min')
