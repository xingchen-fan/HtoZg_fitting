#!/usr/bin/env python3
import ROOT
import os
import sys
import json
import argparse
sys.path.append(os.path.abspath("../Utilities/"))
from bkg_functions_fit import *
from bkg_functions_class import *
from sig_functions_class import *
from Xc_Minimizer import *
from plot_utility import *
from sample_reader import *
from profile_class import *
from sig_functions_class import *
#ROOT.gInterpreter.AddIncludePath('../Utilities/HZGRooPdfs.h')
#ROOT.gSystem.Load('../Utilities/HZGRooPdfs_cxx.so')

# ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--cat', help = 'Category')
parser.add_argument('-con', '--config', help = 'Configuration')
args = parser.parse_args() 
jfile = open('../Config/' + args.config, 'r')
configs = json.load(jfile)
CAT = args.cat
setting = configs[CAT]
lowx = setting["Range"]

# Define variables
x = ROOT.RooRealVar("x", "mllg", lowx, lowx + 65.)
mu_gauss = ROOT.RooRealVar("mu_gauss","always 0"       ,0.)
x.setBins(260)

# SST histograms stored in root files
file_open_bkg = ROOT.TFile.Open('../Data/sst_ggf_hist_nodrop.root', 'READ')
#file_open_sig = ROOT.TFile.Open('../Data/sst_ggf_sig_hist_drop.root', 'READ')
hist_bkg_TH1 = file_open_bkg.Get(CAT+'_dy_sm')
#hist_sig_TH1 = file_open_sig.Get('hist_sig_'+CAT)

hist_bkg = ROOT.RooDataHist('hist_bkg', 'hist_bkg', x, hist_bkg_TH1)
#hist_sig = ROOT.RooDataHist('hist_sig', 'hist_sig', x, hist_sig_TH1)

# Define bkg and sig model
MH = ROOT.RooRealVar("MH","MH"       ,124.7, 120., 130.)
profile = profileClass(x, mu_gauss, CAT, '../Config/'+args.config)

bkg_list = profile.testSelection("Chi2")


# Signal model 
#x.setRange("signal",110, 132)
#res = dscb_model.pdf.fitTo(hist_sig, ROOT.RooFit.SumW2Error(True),\
#                         ROOT.RooFit.Save(True),ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0), ROOT.RooFit.Range("signal"))
#res.Print("V")
#plotClass(x, hist_sig, dscb_model.pdf, bkg_model.SBpdf, "Signal_"+cat)


combine_model = combineSignal(x, MH, CAT, '../Config/config_DSCB.json')
MH.setConstant(True)


# Bkg fit
for bkg_model in bkg_list:
    print ("Starting ", bkg_model.name)
    chi2 = ROOT.RooChi2Var("chi2_bkg_" + CAT, "chi2 bkg " + CAT, bkg_model.pdf,  hist_bkg, ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
    Minimizer_Chi2(chi2, -1, 100, False, 0)
    r = Minimizer_Chi2(chi2, -1, 0.1, False, 0)
    bkg_model.checkBond()
    r.Print("V")

# S+B fit
    c01 = ROOT.RooRealVar("c1", "c1", hist_bkg.sumEntries(), 0, 3.* hist_bkg.sumEntries())
    c02 = ROOT.RooRealVar("c2", "c2", 0., -1000., 1000)

    tot_model = ROOT.RooAddPdf(bkg_model.name + "_tot", bkg_model.name + "_tot", ROOT.RooArgList(combine_model.pdf, bkg_model.pdf), ROOT.RooArgList(c02, c01))
    #combine_model.pdf.Print()
    chi2_tot = ROOT.RooChi2Var("chi2_tot", "chi2_tot", tot_model, hist_bkg, ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson), ROOT.RooFit.Extended(True))
    Minimizer_Chi2(chi2_tot, -1, 100, False, 0)
    res1 = Minimizer_Chi2(chi2_tot, -1, 0.1, False, 0)
    bkg_model.checkBond()
    print("Poisson error:")
    res1.Print("V")
    NS = c02.getVal()
    delta_NS = c02.getError()

    chi2_tot = ROOT.RooChi2Var("chi2_tot", "chi2_tot", tot_model, hist_bkg, ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2), ROOT.RooFit.Extended(True))
    Minimizer_Chi2(chi2_tot, -1, 100, False, 0)
    res2 = Minimizer_Chi2(chi2_tot, -1, 0.1, False, 0)
    bkg_model.checkBond()
    print("SumW2 error:")
    res2.Print("V")
    delta_MC = c02.getError()
    pV = ROOT.Math.chisquared_cdf_c(chi2_tot.getVal(), 260 - res2.floatParsFinal().getSize())
    print("p-value = ", pV)
    print("N signal = ", NS, ", Delta Ns = ", delta_NS, ", Delta MC = ", delta_MC)

    plotClass(x, hist_bkg, tot_model, bkg_model.pdf, "S+B_"+tot_model.GetName(), note='#splitline{N_{S} = ' + '%.2f' % NS + ' #zeta_{S} = ' + '%.2f' % max(0., abs(NS) - 2*delta_MC) + '}{#splitline{#Delta_{MC} = '+ '%.2f' % delta_MC + '}{#splitline{#Delta_{N_{S}} = '+ '%.2f' % delta_NS+'}{P-value = %.2f'%pV +'}}}', CMS = 'Simulation')






