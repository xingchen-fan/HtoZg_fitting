#!/usr/bin/env python3
import ROOT
import os
import sys
import json
import argparse
sys.path.append(os.path.abspath("../Utilities/"))
from bkg_functions_class import *
from sig_functions_class import *
from Xc_Minimizer import *
from plot_utility import *
from sample_reader import *
from profile_class import *
from sig_functions_class import *

#ROOT.gInterpreter.AddIncludePath('../Utilities/HZGRooPdfs.h')
#ROOT.gSystem.Load('../Utilities/HZGRooPdfs_cxx.so')
ROOT.gInterpreter.AddIncludePath('../Utilities/RooGaussStepBernstein.h')
ROOT.gSystem.Load('../Utilities/RooGaussStepBernstein_cxx.so')

# ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--cat', help = 'Category')
parser.add_argument('-conB', '--configB', help = 'Bkg configuration')
parser.add_argument('-conS', '--configS', help = 'Sig configuration')
args = parser.parse_args() 
jfile = open(args.configB, 'r')
jfile_s = open(args.configS, 'r')
configs_s = json.load(jfile_s)
configs = json.load(jfile)
CAT = args.cat
setting = configs[CAT]
lowx = setting["Range"]

# Define variables
x = ROOT.RooRealVar("x", "mllg", lowx, lowx + 65.)
mu_gauss = ROOT.RooRealVar("mu_gauss","always 0"       ,0.)
x.setBins(260)

# SST histograms stored in root files
#file_open_bkg = ROOT.TFile.Open('../Data/sst_'+CAT[:3]+'_hist_nodrop.root', 'READ')
file_open_bkg = ROOT.TFile.Open('~/EOS_space/DY_sample_combine/sst_ggf_hist_xgboost_final.root', 'READ')
#file_open_sig = ROOT.TFile.Open('../Data/sst_ggf_sig_hist_drop.root', 'READ')
hist_bkg_TH1 = file_open_bkg.Get(CAT+'_dy_sm')
#hist_sig_TH1 = file_open_sig.Get('hist_sig_'+CAT)

hist_bkg = ROOT.RooDataHist('hist_bkg', 'hist_bkg', x, hist_bkg_TH1)
#hist_sig = ROOT.RooDataHist('hist_sig', 'hist_sig', x, hist_sig_TH1)

# Define bkg and sig model
MH = ROOT.RooRealVar("MH","MH"       ,125)#, 120., 130.)
profile = profileClass(x, mu_gauss, CAT, args.configB)

bkg_list = profile.testSelection("Chi2")


# Signal model 
#x.setRange("signal",110, 132)
#res = dscb_model.pdf.fitTo(hist_sig, ROOT.RooFit.SumW2Error(True),\
#                         ROOT.RooFit.Save(True),ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0), ROOT.RooFit.Range("signal"))
#res.Print("V")
#plotClass(x, hist_sig, dscb_model.pdf, bkg_model.SBpdf, "Signal_"+cat)


#combine_model = combineSignal(x, MH, CAT, args.configS)
sig_model_el = DSCB_Class(x, MH, CAT+'_el', di_sigma = True)
sig_model_el.assignValDoubleModel(MH, args.configS, CAT, 'el')
sig_model_mu = DSCB_Class(x, MH, CAT+'_mu', di_sigma = True)
sig_model_mu.assignValDoubleModel(MH, args.configS, CAT, 'mu')
c_el = ROOT.RooRealVar('c_el', 'c_el', configs_s[CAT]["el"]["nexp"]/(configs_s[CAT]["el"]["nexp"] + configs_s[CAT]["mu"]["nexp"]))
combine_model = ROOT.RooAddPdf('duo_sig_model_'+CAT, 'duo_sig_model_'+CAT, sig_model_el.pdf, sig_model_mu.pdf, c_el)
NLL = True

# Bkg fit
for bkg_model in bkg_list:
    print ("Starting ", bkg_model.name)
    if NLL:
        r = bkg_model.pdf.fitTo(hist_bkg, ROOT.RooFit.Save(True), \
                               ROOT.RooFit.AsymptoticError(True),ROOT.RooFit.PrintLevel(-1))
        th1_hist = hist_bkg.createHistogram("hist_ks", x, ROOT.RooFit.Binning(260))
        ks_model_hist = bkg_model.pdf.generateBinned(x,hist_bkg.sumEntries(), True).createHistogram("model_hist_ks", x, ROOT.RooFit.Binning(260))
        ks_PV = ks_model_hist.KolmogorovTest(th1_hist, 'X')
        note_ = 'KS P-value = %.2f'%ks_PV
    else:    
        chi2 = ROOT.RooChi2Var("chi2_bkg_" + CAT, "chi2 bkg " + CAT, bkg_model.pdf,  hist_bkg, ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))
        Minimizer_Chi2(chi2, -1, 100, False, 0)
        r = Minimizer_Chi2(chi2, -1, 0.1, False, 0)
        note_ = '#Chi^{2}/dof = %.2f'%chi2.getVal()+'/%i'%(260-r.floatParsFinal().getSize())
    bkg_model.checkBond()
    r.Print("V")
    plotClass(x, hist_bkg, bkg_model.pdf, bkg_model.SBpdf, "B_"+bkg_model.name,  note=note_, CMS = 'Simulation', fullRatio = True)
# S+B fit
    c01 = ROOT.RooRealVar("c1", "c1", hist_bkg.sumEntries(), 0, 3.* hist_bkg.sumEntries())
    c02 = ROOT.RooRealVar("c2", "c2", 0., -1000., 1000)
    tot_model = ROOT.RooAddPdf(bkg_model.name + "_tot", bkg_model.name + "_tot", ROOT.RooArgList(combine_model, bkg_model.pdf), ROOT.RooArgList(c02, c01))
    if NLL:
        res1 = tot_model.fitTo(hist_bkg, ROOT.RooFit.Save(True), ROOT.RooFit.SumW2Error(False), ROOT.RooFit.Strategy(0), ROOT.RooFit.PrintLevel(-1))
        note_chi2 = ''
    else:
        chi2_tot = ROOT.RooChi2Var("chi2_tot", "chi2_tot", tot_model, hist_bkg, ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson), ROOT.RooFit.Extended(True))
        chi2_val = chi2_tot.getVal()
        Minimizer_Chi2(chi2_tot, -1, 100, False, 0)
        res1 = Minimizer_Chi2(chi2_tot, -1, 0.1, False, 0)
        pV = ROOT.Math.chisquared_cdf_c(chi2_val, 260 - res1.floatParsFinal().getSize())
        note_chi2 = '#Chi^{2} = %.2f'%chi2_val+',P-value = %.2f'%pV
    bkg_model.checkBond()
    print("Poisson error:")
    res1.Print("V")
    NS = c02.getVal()
    delta_NS = c02.getError()
    #plotClass(x, hist_bkg, tot_model, bkg_model.pdf, "S+B_Poisson_"+tot_model.GetName(), note='', CMS = 'Simulation', fullRatio = True)

    if NLL:
        res2 = tot_model.fitTo(hist_bkg, ROOT.RooFit.Save(True), ROOT.RooFit.AsymptoticError(True), ROOT.RooFit.Strategy(0), ROOT.RooFit.PrintLevel(-1))
    else:
        chi2_tot = ROOT.RooChi2Var("chi2_tot", "chi2_tot", tot_model, hist_bkg, ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2), ROOT.RooFit.Extended(True))
        Minimizer_Chi2(chi2_tot, -1, 100, False, 0)
        res2 = Minimizer_Chi2(chi2_tot, -1, 0.1, False, 0)
    bkg_model.checkBond()
    print("SumW2 error:")
    res2.Print("V")
    delta_MC = c02.getError()
    print("N signal = ", NS, ", Delta Ns = ", delta_NS, ", Delta MC = ", delta_MC)

    plotClass(x, hist_bkg, tot_model, bkg_model.pdf, "S+B_"+tot_model.GetName(), note='#splitline{N_{S} = ' + '%.2f' % NS + ' #zeta_{S} = ' + '%.2f' % max(0., abs(NS) - 2*delta_MC) + '}{#splitline{#Delta_{MC} = '+ '%.2f' % delta_MC + '}{#splitline{#Delta_{N_{S}} = '+ '%.2f' % delta_NS+'}{' + note_chi2 +'}}}', CMS = 'Simulation', fullRatio = True)






