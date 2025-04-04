#!/usr/bin/env python3
import os
import sys
import ROOT
import argparse
import json
sys.path.append(os.path.abspath("../Utilities/"))
sys.path.append(os.path.abspath("../CMS_plotter/"))
import CMS_lumi, tdrstyle
from bkg_functions_fit import *
from bkg_functions_class import *
from Xc_Minimizer import *
from plot_utility import *
from sample_reader import *
from profile_class import *
from APFUtilities import *
ROOT.gInterpreter.AddIncludePath('../Utilities/HZGRooPdfs.h')
ROOT.gSystem.Load('../Utilities/HZGRooPdfs_cxx.so')

ROOT.gROOT.SetBatch(True)


#run command:python3 VHttH_Chi2_test.py -c ttH_lep -con ./Config/vhtth_config.json

ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--cat', help = 'Category')
parser.add_argument('-con', '--config', help = 'Configuration')
args = parser.parse_args()
CAT = args.cat
jfile = open('../Config/'+args.config, 'r')
configs = json.load(jfile)
setting = configs[CAT]

# Define vars
lowx = setting["Range"]

#Get m_lly variable with correct range and sidebands
x = define_mlly()

#Define file name
file_name_part1 = "an_int_fit_dist_Data_cat_"
file_name_part2 = "_AllYears__llphoton_refit_m__wgt__lumi_nonorm_lin.root"
file_name = VHttH_file_path + file_name_part1 + CAT + file_name_part2

#Get histogram from input file
hist_data = get_data_rdh(file_name, "rdh_" + CAT, x)

#Print number of events in the sideband
print('N data = ', hist_data.sumEntries())

#Bkg funcs
mu_gauss = ROOT.RooRealVar("mu_gauss", "always 0" ,0.)

#Get Profile from config
profile = profileClass(x, mu_gauss, CAT, '../Config/'+args.config)

bkg_list = [profile.modg_model, profile.exmg_model, profile.agg_model, profile.bern2_model, profile.bern3_model, profile.bern4_model, profile.bern5_model, profile.pow1_model, profile.pow2_model, profile.pow3_model, profile.exp2_model, profile.exp3_model, profile.lau2_model, profile.lau3_model, profile.lau4_model]
#bkg_list = [profile.modg_model, profile.exmg_model, profile.agg_model]
#profile.bern2_model, profile.modg_model, 

# Sideband Chi^2 fit
cuthistogram = hist_data.reduce(ROOT.RooFit.CutRange('left,right'))
th1_cuthist = cuthistogram.createHistogram("cuthist_ks", x, ROOT.RooFit.Binning(160))
n_bins_sideband = 140
HIGHSTAT = CAT == 'untagged'
pass_list = []

#Fit pass criteria
chi2_requirement = 0.1
ks_requirement = 0.2

for func in bkg_list:
    #Only high stat cat is untagged
    if HIGHSTAT:
        chi2 = ROOT.RooChi2Var('chi2_'+func.SBpdf.GetName(), 'chi2_'+func.SBpdf.GetName(), func.SBpdf, cuthistogram)
    else:
        chi2 = ROOT.RooNLLVar('chi2_'+func.SBpdf.GetName(), 'chi2_'+func.SBpdf.GetName(), func.SBpdf, cuthistogram)

    #
    Minimizer_Chi2(chi2, -1, 1, False, 0)
    res=Minimizer_Chi2(chi2, -1, 0.1, True, 0) 
    func.checkBond()
    res.Print('v')
    
    #Get Chi squared value per degrees of freedom
    chi2_ = ROOT.RooChi2Var('chi2_val_'+func.SBpdf.GetName(), 'chi2_val_'+func.SBpdf.GetName(), func.SBpdf, cuthistogram)
    chi2_val = chi2_.getVal()
    dof = n_bins_sideband - res.floatParsFinal().getSize()
    chi2_PV = ROOT.Math.chisquared_cdf_c(chi2_.getVal(), dof)
    
    #Get KS test value
    ks_model_hist = func.SBpdf.generateBinned(x, cuthistogram.sumEntries(), True).createHistogram("hist_"+func.SBpdf.GetName(), x, ROOT.RooFit.Binning(160))
    ks_PV = th1_cuthist.KolmogorovTest(ks_model_hist)
    
    #Plot individual fit on plot
    plotClass(x, hist_data, func.pdf, func.SBpdf, title=func.pdf.GetName(), output_dir="plots/", sideBand = True, fitRange = 'left,right', note='#splitline{#Chi^{2}/dof = ' + '%.2f'%chi2_val + '/%i'%dof +'}{#splitline{#Chi^{2} P-value = ' + '%.2f'%chi2_PV + '}{KS P-value = ' + '%.2f'%ks_PV + '}}', fullRatio = ("vbf" in CAT))
    
    #Checkinging whether or not the fit passes the criteria of the test
    if chi2_PV> chi2_requirement and ks_PV> ks_requirement and res.status() == 0:
        pass_list.append(func)

#Output number of functions passing and the plot with all the functions passing
print('List has ', len(pass_list), ' models')
multiPlotClass(x, hist_data, pass_list, title="Chi2_"+CAT, output_dir="plots/",sideBand = True, fitRange = 'left,right',best_index = 0,CMS = "Preliminary", fullRatio = ("vbf" in CAT), sideband_bins=80, ratio_range=[0.5,1.5])







'''


mlly_low = 100
mlly_high = 180
mlly_sideband_low = 120
mlly_sideband_high = 130

x = ROOT.RooRealVar("x", "mllg", mlly_low, mlly_high)
#y = ROOT.RooRealVar("y", "photon pT", 15., 1000.)
#w = ROOT.RooRealVar("w", "w", -40., 40.)
#bdt = ROOT.RooRealVar("bdt", "bdt", -1, 1)
#year = ROOT.RooRealVar("year", "year", 2015, 21000)
#w_year = ROOT.RooRealVar("w_year", "w_year", 0, 100)
#lep = ROOT.RooRealVar("lep", "lep", 0, 1) #0 = electron, 1 = muon
#ph_eta = ROOT.RooRealVar("ph_eta", "ph_eta", -3, 3)
#nlep = ROOT.RooRealVar("nlep", "nlep", 0, 10)
#njet = ROOT.RooRealVar("njet", "njet", 0, 10)

#x.https://github.com/xingchen-fan/HtoZg_fitting.gitsetBins(160)
x.setRange('left', mlly_low, mlly_sideband_low)
x.setRange('right', mlly_sideband_high, mlly_high)
x.setRange('full', mlly_low, mlly_high)


if DAT:
    read_data = ROOT.RooDataSet.read('../Data/data_'+NAMECAT+'_pinnacles_fix.dat', ROOT.RooArgList(x, y, bdt, w, w_year, year, lep, ph_eta, nlep, njet))
    data = ROOT.RooDataSet('data', 'data', read_data, ROOT.RooArgList(x, y, bdt, w, w_year, year, lep, ph_eta, nlep, njet),'')
    hist_data = ROOT.RooDataHist('hist_data','hist_data', x, data)
else:
    if 'ggf' in CAT:
        read_data = readRuiROOTggFdata(x, '/eos/user/r/rzou//SWAN_projects/Classifier/Output_ggF_pinnacles_fix/relpt_peking_run2p3/TreeB/', 0.81, 0.64, 0.47)
        if CAT == 'ggf1':
            hist_data = read_data.ggf1
        elif CAT == 'ggf2':
            hist_data = read_data.ggf2
        elif CAT == 'ggf3':
            hist_data = read_data.ggf3
        elif CAT == 'ggf4':
            hist_data = read_data.ggf4
    elif 'vbf' in CAT:
        read_data = readRuiROOTVBFdata(x, '/eos/user/r/rzou/SWAN_projects/Classifier/Output_2JClassic_input_run2p3/', 0.489,0.286 , 0.083)
        if CAT == 'vbf1':
            hist_data = read_data.vbf1
        elif CAT == 'vbf2':
            hist_data = read_data.vbf2
        elif CAT == 'vbf3':
            hist_data = read_data.vbf3
        elif CAT == 'vbf4':
            hist_data = read_data.vbf4

#profile.lau2_model, profile.lau3_model, profile.lau4_model,
#vbf1 no lau3, vbf2 no lau4
#AGG_model = AGGClass(x, CAT, kappa_init = -1.27, alpha_init = 14, zeta_init = 105, x_low = lowx, x_high = lowx+65)




'''
