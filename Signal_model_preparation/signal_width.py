#!/usr/bin/env python3
import os
import sys
import ROOT
import argparse
import json
sys.path.append(os.path.abspath("../Utilities/"))
sys.path.append(os.path.abspath("../CMS_plotter/"))
import CMS_lumi, tdrstyle
from sig_functions_class import *
from Xc_Minimizer import *
from plot_utility import *
from sample_reader import *

ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--cat', help = 'Category')
parser.add_argument('-p', '--prod', help = 'Production mode')
parser.add_argument('-con', '--con', help = 'Config')

args = parser.parse_args()
CAT = args.cat
MH =  ROOT.RooRealVar("MH", "MH", 125)
x = ROOT.RooRealVar("x", "mllg", 110, 135)
years = ["2016", "2016APV", "2017", "2018", "2022", "2022EE", "2023BPix", "2023"]
for lep in ["el", "mu"]:
    print (lep, ":")
    for year in years:
        signal = DSCB_Class(x, MH, di_sigma=True)
        signal.assignVal(MH, args.con, year, CAT, lep, args.prod)
        pair = getEffSigma(x, signal.pdf)
        print('%.2f'%((pair[1] - pair[0])/2))

        
