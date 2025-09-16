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
from bkg_functions_class import *
from Xc_Minimizer import *
from plot_utility import *
from sample_reader import *
#ROOT.gInterpreter.AddIncludePath('../Utilities/HZGRooPdfs.h')
#ROOT.gSystem.Load('../Utilities/HZGRooPdfs_cxx.so')

ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--cat', help = 'Cat')
parser.add_argument('-r', '--ranges', help = 'Range')
parser.add_argument('-hi', '--high', help = 'Higher power')
parser.add_argument('-lo', '--low', help = 'Lower power')
parser.add_argument('-NLL', '--NLL', help = 'Use NLL?', type=int, default=1)
parser.add_argument('-d', '--duo', help = 'Double sigma', type=int, default = 0)
parser.add_argument('-f', '--fix', help = 'Fix sigma', type=int, default = 0)

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
#raise Exception('Wrong category name!')
 
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

DAT = False
# Read samples (dat)
if DAT:
    read_data = ROOT.RooDataSet.read('data_' +CATNAME+'_pinnacles_fix.dat', ROOT.RooArgList(x, y, bdt, w, w_year, year, lep, ph_eta, nlep, njet))
    data = ROOT.RooDataSet('data', 'data', read_data, ROOT.RooArgList(x, y, bdt, w, w_year, year, lep, ph_eta, nlep, njet),'')
    hist_data = ROOT.RooDataHist('hist_data','hist_data', x, data)

else:
    if 'ggf' in CAT:
        read_data = readRuiROOTggFdata(x, '/eos/user/r/rzou/SWAN_projects/Classifier/Output_ggF_rui_commonparam/', 0.91,0.82,0.61)
        if CAT == 'ggf1':
            hist_data = read_data.ggf1
        elif CAT == 'ggf2':
            hist_data = read_data.ggf2
        elif CAT == 'ggf3':
            hist_data = read_data.ggf3
        elif CAT == 'ggf4':
            hist_data = read_data.ggf4
    elif 'vbf' in CAT:
        read_data = readRuiROOTVBFdata(x, '/eos/user/r/rzou/SWAN_projects/Classifier/Output_VBF_rui_commonparam/', 0.95, 0.91,0.76)
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

list_ = []
stat_list = []
chi2_list = []
pv_list = []

ONESHOT = False
if ONESHOT:
    lau2_model = Lau2Class(x, mu_gauss, CAT, sigma_init = 3., sigma2_init = 4., step_init = 101, p1 = int(args.low), p2 = int(args.high), f1_init = 0.2, f2_init = 0.1, xmax = lowx+65., const_f1 = True, di_gauss = args.duo, fix_sigma = args.fix, gc_init = 1)
    cuthistogram = hist_data.reduce(ROOT.RooFit.CutRange('left,right'))
    n_bins = 220
    if args.NLL == True:
        chi2 = ROOT.RooNLLVar('NLL_'+lau2_model.SBpdf.GetName(), 'NLL_'+lau2_model.SBpdf.GetName(), lau2_model.SBpdf, cuthistogram)
    else:
        chi2 = ROOT.RooChi2Var('chi2_'+lau2_model.SBpdf.GetName(), 'chi2_'+lau2_model.SBpdf.GetName(), lau2_model.SBpdf, cuthistogram)
        
    Minimizer_Chi2(chi2, -1, 100, False, 0)
    res=Minimizer_Chi2(chi2, -1, 0.1, True, 0)
    lau2_model.checkBond()
    res.Print('v')
    chi2_= ROOT.RooChi2Var('chi2_val_'+lau2_model.SBpdf.GetName(), 'chi2_val_'+lau2_model.SBpdf.GetName(), lau2_model.SBpdf, cuthistogram)
    chi2_val = chi2_.getVal()
    chi2_pV= ROOT.Math.chisquared_cdf_c(chi2_.getVal(), n_bins - res.floatParsFinal().getSize())
    print ('Chi2 = ', chi2_val)
    print ('P-value = ', chi2_pV)
    plotClass(x, hist_data, lau2_model.pdf, lau2_model.SBpdf, "Lau2", "", True, 'left,right')
    sys.exit()
for i in range(int(args.low), int(args.high) + 1):
    for j in range(i+1, int(args.high) + 1):
        pow_list = []
        pow_list.append(i)
        pow_list.append(j)
        lau2_model = Lau2Class(x, mu_gauss, CAT, sigma_init = 3., sigma2_init = 4., step_init = 101, p1 = i, p2 = j, f1_init = 0.2, f2_init = 0.1, xmax = lowx+65., const_f1 = True, di_gauss = args.duo, fix_sigma = False, gc_init = 1)
        cuthistogram = hist_data.reduce(ROOT.RooFit.CutRange('left,right'))
        n_bins = 220
        if args.NLL == True:
            chi2 = ROOT.RooNLLVar('NLL_'+lau2_model.SBpdf.GetName(), 'NLL_'+lau2_model.SBpdf.GetName(), lau2_model.SBpdf, cuthistogram)
        else:
            chi2 = ROOT.RooChi2Var('chi2_'+lau2_model.SBpdf.GetName(), 'chi2_'+lau2_model.SBpdf.GetName(), lau2_model.SBpdf, cuthistogram)
        Minimizer_Chi2(chi2, -1, 100, False, 0)
        #Minimizer_Chi2(chi2, -1, 0.1, True, 0)
        res=Minimizer_Chi2(chi2, -1, 0.1, True, 0) 
        lau2_model.checkBond()
        res.Print('v')
        chi2_ = ROOT.RooChi2Var('chi2_val_'+lau2_model.SBpdf.GetName(), 'chi2_val_'+lau2_model.SBpdf.GetName(), lau2_model.SBpdf, cuthistogram)
        chi2_pV = ROOT.Math.chisquared_cdf_c(chi2_.getVal(), n_bins - res.floatParsFinal().getSize())
        stat_list.append(chi2.getVal())
        chi2_list.append(chi2_.getVal())
        pv_list.append(chi2_pV)
        list_.append(pow_list)
        print ('Finish', pow_list)
min_stat = min(stat_list)
index_ = stat_list.index(min_stat)
print ('Power combinations: ', list_)
print ('Test stat list: ', stat_list)
print ('Chi^2 list: ', chi2_list)
print ('Best lau2 ', list_[index_], ', Chi^2 = ', chi2_list[index_], ', P-value = ', pv_list[index_])
toc = time.perf_counter()
print ('Time spent:', '%.2f'%((toc - tic)/60), 'min')
