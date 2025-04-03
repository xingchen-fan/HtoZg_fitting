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

ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--cat', help = 'Category')
parser.add_argument('-con', '--config', help = 'Configuration')
parser.add_argument('-t', '--test', help = 'Test name')
parser.add_argument('-b', '--best', help = 'Best func index', default = 0)
args = parser.parse_args()
CAT = args.cat
jfile = open('../Config/'+args.config, 'r')
configs = json.load(jfile)
setting = configs[CAT]
lowx = setting["Range"]
x = ROOT.RooRealVar("x", "mllg", lowx, lowx+65.)
x.setBins(260)
x.setRange('left', lowx, 120)
x.setRange('right', 130, lowx+65)
x.setRange('full', lowx, lowx+65)
mu_gauss = ROOT.RooRealVar("mu_gauss","always 0"       ,0.)

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
print('N data = ', hist_data.sumEntries())

BEST = False
if args.test == "FT": BEST = True
profile = profileClass(x, mu_gauss, CAT, '../Config/'+args.config)
bkg_list = profile.testSelection(args.test)

multiPlotClass(x, hist_data, bkg_list, title=args.test+'_'+CAT, output_dir="plots/",sideBand = True, fitRange = 'left,right',best_index = int(args.best), CMS = "Preliminary", fullRatio = True, bestLabel = BEST)
