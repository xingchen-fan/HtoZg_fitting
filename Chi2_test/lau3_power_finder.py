#!/usr/bin/env python3 
import time
import os
import sys
import ROOT
import json
sys.path.append(os.path.abspath("../Utilities/"))
sys.path.append(os.path.abspath("../CMS_plotter/"))
import CMS_lumi, tdrstyle
from bkg_functions_fit import *
from bkg_functions_class import *
from Xc_Minimizer import *
from plot_utility import *
from sample_reader import *
#ROOT.gInterpreter.AddIncludePath('../Utilities/HZGRooPdfs.h')
#ROOT.gSystem.Load('../Utilities/HZGRooPdfs_cxx.so')
import argparse
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)
parser = argparse.ArgumentParser()
parser.add_argument('-hi', '--high', help = 'Higher power')
parser.add_argument('-lo', '--low', help = 'Lower power')
parser.add_argument('-c', '--cat', help = 'Cat')
parser.add_argument('-r', '--ranges', help = 'Range')
parser.add_argument('-NLL', '--NLL', help = 'Use NLL?', type=int, default=0)

args = parser.parse_args()
if int(args.high) <= int(args.low):
    raise Exception('Higher power must be larger than lower power!')


tic = time.perf_counter()
# Define vars
CAT = args.cat
if args.cat == 'ggf1':
    CATNAME = 'ggF1'
elif args.cat == 'ggf2':
    CATNAME = 'ggF2'
elif args.cat == 'ggf3':
    CATNAME = 'ggF3'
elif args.cat == 'ggf4':
    CATNAME = 'ggF4'

lowx = int(args.ranges)
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
DAT = False
if DAT:
    read_data = ROOT.RooDataSet.read('data_'+CATNAME+'_pinnacles_fix.dat', ROOT.RooArgList(x, y, bdt, w, w_year, year, lep, ph_eta, nlep, njet))
    data = ROOT.RooDataSet('data', 'data', read_data, ROOT.RooArgList(x, y, bdt, w, w_year, year, lep, ph_eta, nlep, njet),'')
    hist_data = ROOT.RooDataHist('hist_data','hist_data', x, data)
else:
    if 'ggf' in CAT:
        read_data = readRuiROOTggFdata(x,  '/eos/user/r/rzou//SWAN_projects/Classifier/Output_ggF_pinnacles_fix/relpt_peking_run2p3/TreeB/', 0.81, 0.64, 0.47)
        if CAT == 'ggf1':
            hist_data = read_data.ggf1
        elif CAT == 'ggf2':
            hist_data = read_data.ggf2
        elif CAT == 'ggf3':
            hist_data = read_data.ggf3
        elif CAT == 'ggf4':
            hist_data = read_data.ggf4
    elif 'vbf' in CAT:
        read_data = readRuiROOTVBFdata(x, '/eos/user/r/rzou/SWAN_projects/Classifier/Output_2JClassic_input_run2p3/', 0.489,0.286, 0.083)
        if CAT == 'vbf1':
            hist_data = read_data.vbf1
        elif CAT == 'vbf2':
            hist_data = read_data.vbf2
        elif CAT == 'vbf3':
            hist_data = read_data.vbf3
        elif CAT == 'vbf4':
            hist_data = read_data.vbf4

# Bkg funcs
mu_gauss = ROOT.RooRealVar("mu_gauss","always 0"       ,0.)
alt_pow = []
if int(args.high) - int(args.low) == 1:
    alt_pow.append(['hi', 1])
    alt_pow.append(['hi', 2])
    alt_pow.append(['lo', -1])
elif int(args.high) - int(args.low) == 2:
    alt_pow.append(['hi', 1])
    alt_pow.append(['hi', -1])
    alt_pow.append(['lo', -1])
else:
    alt_pow.append(['hi', 1])
    alt_pow.append(['hi', -1])
    alt_pow.append(['lo', 1])
    alt_pow.append(['lo', -1])
list_ = []
chi2_list = []
pv_list = []
for entry in alt_pow:
    if entry[0] == 'hi':
        pow_list = [int(args.high), int(args.low), int(args.high) + entry[1]]
    else:
        pow_list = [int(args.high), int(args.low), int(args.low) + entry[1]]

    lau3_model = Lau3Class(x, mu_gauss, CAT, sigma_init = 3.,sigma2_init = 4,  step_init = 101, p1 = pow_list[0], p2 = pow_list[1], p3 = pow_list[2], f1_init = 0.1, f2_init = 0.1, f3_init = 0.1, xmax = lowx+65., const_f1 = True, di_gauss = False, fix_sigma = False, gc_init = 1)
    cuthistogram = hist_data.reduce(ROOT.RooFit.CutRange('left,right'))
    n_bins = 220
    if args.NLL == True:
        chi2 = ROOT.RooNLLVar('NLL_'+lau3_model.SBpdf.GetName(), 'NLL_'+lau3_model.SBpdf.GetName(), lau3_model.SBpdf, cuthistogram)
    else:
        chi2 = ROOT.RooChi2Var('chi2_'+lau3_model.SBpdf.GetName(), 'chi2_'+lau3_model.SBpdf.GetName(), lau3_model.SBpdf, cuthistogram)
    Minimizer_Chi2(chi2, -1, 100, False, 0)
    #Minimizer_Chi2(chi2, -1, 1, True, 0)
    res=Minimizer_Chi2(chi2, -1, 0.1, True, 0) 
    lau3_model.checkBond()
    res.Print('v')
    chi2_ = ROOT.RooChi2Var('chi2_val_'+lau3_model.SBpdf.GetName(), 'chi2_val_'+lau3_model.SBpdf.GetName(), lau3_model.SBpdf, cuthistogram)
    chi2_pV = ROOT.Math.chisquared_cdf_c(chi2_.getVal(), n_bins - res.floatParsFinal().getSize())
    chi2_list.append(chi2.getVal())
    pv_list.append(chi2_pV)
    list_.append(pow_list)
    print ('Finish', pow_list)
min_chi2 = min(chi2_list)
index_ = chi2_list.index(min_chi2)
print ('Power combinations: ', list_)
print ('Chi^2 (NLL) list: ', chi2_list)
print ('Best lau3 ', list_[index_], ', Chi^2 (NLL) = ', chi2_list[index_], ', P-value = ', pv_list[index_])
toc = time.perf_counter()
print ('Time spent:', '%.2f'%((toc - tic)/60), 'min')
