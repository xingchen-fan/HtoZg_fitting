#!/usr/bin/env python3
import ROOT
import os
import sys
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
from wip_profile_class import *
from APFUtilities import *
ROOT.gInterpreter.AddIncludePath('../Utilities/HZGRooPdfs.h')
ROOT.gSystem.Load('../Utilities/HZGRooPdfs_cxx.so')

#ROOT.gROOT.SetBatch(True)

ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)

parser = argparse.ArgumentParser(description = "F test bin and function class")
#parser.add_argument("method")
parser.add_argument('-c', '--cat', help="category")
parser.add_argument('-con', '--config', help = 'Configuration')

args = parser.parse_args()
jfile = open('../Config/'+args.config, 'r')
configs = json.load(jfile)
CAT = args.cat
setting = configs[CAT]

# Define variables
x = define_mlly()

#Define the mean of a gaussian at 0
mu_gauss = ROOT.RooRealVar("mu_gauss","always 0", 0.)

#Define file name
file_name_part1 = "an_int_fit_dist_Data_cat_"
file_name_part2 = "_AllYears__llphoton_refit_m__wgt__lumi_nonorm_lin.root"
file_name = VHttH_file_path + file_name_part1 + CAT + file_name_part2

#Get data histogram as a RooDataHist
data_hist = get_data_rdh(file_name, "rdh_" + CAT, x)

profile = profileClass(x, mu_gauss, CAT, '../Config/'+args.config)

input_key = "PostFT"
envelope_funcs = setting[input_key]

pdfList = []
if CAT=="ttH_lep":
  pdfList = [profile.exmg_model, profile.agg_model, profile.bern2_model, profile.pow1_model,  profile.lau2_model]
  #pdfList = [profile.exmg_model, profile.agg_model]
  #pdfList = [profile.agg_model, profile.bern2_model, profile.pow1_model,  profile.lau2_model]

if CAT=="ttH_had":
  pdfList = [profile.exmg_model, profile.agg_model, profile.bern2_model, profile.pow1_model,  profile.lau2_model]
if CAT=="WH_3l":
  pdfList = [profile.modg_model, profile.exmg_model, profile.agg_model, profile.bern2_model, profile.pow1_model,  profile.lau2_model]
if CAT=="ZH_MET":
  pdfList = [profile.exmg_model, profile.agg_model, profile.bern2_model, profile.pow1_model,  profile.lau2_model]
if CAT=="untagged":
  pdfList = [profile.modg_model, profile.agg_model, profile.bern2_model, profile.bern3_model, profile.pow2_model, profile.exp2_model]

# Sideband Chi^2 fit
cuthistogram = data_hist.reduce(ROOT.RooFit.CutRange('left,right'))
th1_cuthist = cuthistogram.createHistogram("cuthist_ks", x, ROOT.RooFit.Binning(160))
n_bins_sideband = 140
HIGHSTAT = CAT == 'untagged'
pass_list = []

#Fit pass criteria
chi2_requirement = 0.1
ks_requirement = 0.2
HIGHSTAT = (CAT=="untagged")

for func in pdfList:
    #Only high stat cat is untagged
    if HIGHSTAT:
        chi2 = ROOT.RooChi2Var('chi2_'+func.SBpdf.GetName(), 'chi2_'+func.SBpdf.GetName(), func.SBpdf, cuthistogram)
    else:
        chi2 = ROOT.RooNLLVar('chi2_'+func.SBpdf.GetName(), 'chi2_'+func.SBpdf.GetName(), func.SBpdf, cuthistogram)

    Minimizer_Chi2(chi2, -1, 1, False, 0)
    res=Minimizer_Chi2(chi2, -1, 0.1, True, 0) 
    func.checkBond()
    res.Print('v')

#Make Plot
multiPlotClass(x, data_hist, pdfList, title="PostFTest_multi_" + CAT, output_dir="./plots/", sideBand=True, fitRange="left,right", best_index = 1,  sideband_bins=80, ratio_range=[0.5,1.5], fullRatio=False)

profile.write_config_file(data_hist, "PostFT", config_out=CAT + '_postFT_config.json')

