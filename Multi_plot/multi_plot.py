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
ROOT.gInterpreter.AddIncludePath('../Utilities/RooGaussStepBernstein.h')
ROOT.gSystem.Load('../Utilities/RooGaussStepBernstein_cxx.so')

ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--cat', help = 'Category')
parser.add_argument('-con', '--config', help = 'Configuration')
parser.add_argument('-t', '--test', help = 'Test name')
parser.add_argument('-b', '--best', help = 'Best func index', default = -1, type=int)
args = parser.parse_args()
CAT = args.cat
jfile = open(args.config, 'r')
configs = json.load(jfile)
setting = configs[CAT]
lowx = setting["Range"][0]
highx = setting["Range"][1]
nbins = int(setting["Bins"])

x = ROOT.RooRealVar("x", "x", lowx, highx)
x.setBins(nbins)
x.setRange('left', lowx, 120)
x.setRange('right', 130, highx)
x.setRange('full', lowx, highx)
mu_gauss = ROOT.RooRealVar("mu_gauss","always 0"       ,0.)

if 'ggf' in CAT:
        read_data = readRuiROOTggFdata(x, '/eos/project/h/htozg-dy-privatemc/rzou/bdt/BDT_output_redwood/Output_ggF_rui_redwood_v1_ext_val/', 0.94,0.83,0.57)
        if CAT == 'ggf1':
            hist_data = read_data.ggf1
        elif CAT == 'ggf2':
            hist_data = read_data.ggf2
        elif CAT == 'ggf3':
            hist_data = read_data.ggf3
        elif CAT == 'ggf4':
            hist_data = read_data.ggf4
elif 'vbf' in CAT:
        read_data = readRuiROOTVBFdata(x, '/eos/project/h/htozg-dy-privatemc/rzou/bdt/BDT_output_redwood/Output_VBF_rui_redwood_v1_ext_val/', 0.91, 0.81,0.48)
        if CAT == 'vbf1':
            hist_data = read_data.vbf1
        elif CAT == 'vbf2':
            hist_data = read_data.vbf2
        elif CAT == 'vbf3':
            hist_data = read_data.vbf3
        elif CAT == 'vbf4':
            hist_data = read_data.vbf4

print('N data = ', hist_data.sumEntries())
profile = profileClass(x, mu_gauss, CAT, args.config)
bkg_list = profile.testSelection(args.test)
multiPlotClass(x, hist_data, bkg_list, title=args.test+'_'+CAT, output_dir="plots/",sideBand = True, fitRange = 'left,right',best_index = args.best, CMS = "Preliminary", fullRatio = False, ratio_range=[0, 4], bestLabel = args.best > -1, leg_text_size = 0.06)
