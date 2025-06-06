#!/usr/bin/env python3
import sys
import os
import ROOT
import argparse
import json
sys.path.append(os.path.abspath("../Utilities/"))
sys.path.append(os.path.abspath("../CMS_plotter/"))
import CMS_lumi, tdrstyle
from Xc_Minimizer import *
from plot_utility import *
from profile_class import *
from sig_functions_class import *
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--cat', help = 'Category')
parser.add_argument('-conB', '--configB', help = 'Bkg Configuration')
parser.add_argument('-conS', '--configS', help = 'Sig Configuration')
parser.add_argument('-i', '--input', help = 'Input Root File')

args = parser.parse_args()
CAT = args.cat
jfile = open(args.configB, 'r')
configs = json.load(jfile)
setting = configs[CAT]
lowx = setting["Range"]

f = ROOT.TFile(args.input)
w = f.Get("w")
data = w.data("data_obs")
x = w.var('CMS_hzg_mass_'+CAT)
mu_gauss = ROOT.RooRealVar("mu_gauss","always 0"       ,0.)
x.setRange('left', lowx, 120)
x.setRange('right', 130, lowx+65)
x.setBins(260)
hdata = ROOT.RooDataHist('hdata', 'hdata', x, data)

MH = ROOT.RooRealVar("MH","MH"       ,125)#, 120., 130.)
combine_model = combineSignal(x, MH, CAT, args.configS)
profile = profileClass(x, mu_gauss, CAT, args.configB)
bkg_list = profile.testSelection("Best")

for bkg_model in bkg_list:
    print ("Starting ", bkg_model.name)
    c01 = ROOT.RooRealVar("c1", "c1", hdata.sumEntries(), 0, 3.* hdata.sumEntries())
    c02 = ROOT.RooRealVar("c2", "c2", 0., -1000., 1000)
    tot_model = ROOT.RooAddPdf(bkg_model.name + "_tot", bkg_model.name + "_tot", ROOT.RooArgList(combine_model.pdf, bkg_model.pdf), ROOT.RooArgList(c02, c01))
    res1 = tot_model.fitTo(hdata, ROOT.RooFit.Save(True), ROOT.RooFit.SumW2Error(False), ROOT.RooFit.Strategy(0), ROOT.RooFit.PrintLevel(-1))
    res1.Print("V")
profile.write_config_file(hdata, "Best")
