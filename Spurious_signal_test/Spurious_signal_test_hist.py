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
ROOT.gInterpreter.AddIncludePath('../Utilities/HZGRooPdfs.h')
ROOT.gSystem.Load('../Utilities/HZGRooPdfs_cxx.so')

# ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--cat', help = 'Category')
args = parser.parse_args() 
jfile = open('../Config/config.json', 'r')
configs = json.load(jfile)
CAT = args.cat
setting = configs[CAT]
lowx = setting["Range"]

# Define variables
x = ROOT.RooRealVar("x", "mllg", lowx, lowx + 65.)
mu_gauss = ROOT.RooRealVar("mu_gauss","always 0"       ,0.)
x.setBins(260)

# SST histograms stored in root files
file_open_bkg = ROOT.TFile.Open('../Data/sst_ggf_hist_drop.root', 'READ')
file_open_sig = ROOT.TFile.Open('../Data/sst_ggf_sig_hist_drop.root', 'READ')
hist_bkg_TH1 = file_open_bkg.Get(CAT+'_dy_sm')
hist_sig_TH1 = file_open_sig.Get('hist_sig_'+CAT)

hist_bkg = ROOT.RooDataHist('hist_bkg', 'hist_bkg', x, hist_bkg_TH1)
hist_sig = ROOT.RooDataHist('hist_sig', 'hist_sig', x, hist_sig_TH1)

# Define bkg and sig model
MH = ROOT.RooRealVar("MH","MH"       ,124.7, 120., 130.)

bern2_model = Bern2Class(x, mu_gauss, CAT, setting["bern2"]['p0'], setting["bern2"]['p_init'],setting["bern2"]['bond'],setting["bern2"]['sigma_init'],setting["bern2"]['step_init'])
bern3_model = Bern3Class(x, mu_gauss, CAT, setting["bern3"]['p0'], setting["bern3"]['p_init'],setting["bern3"]['bond'],setting["bern3"]['sigma_init'],setting["bern3"]['step_init'])
bern4_model = Bern4Class(x, mu_gauss, CAT, setting["bern4"]['p0'], setting["bern4"]['p_init'],setting["bern4"]['bond'],setting["bern4"]['sigma_init'],setting["bern4"]['step_init'])
bern5_model = Bern5Class(x, mu_gauss, CAT, setting["bern5"]['p0'], setting["bern5"]['p_init'],setting["bern5"]['bond'],setting["bern5"]['sigma_init'],setting["bern5"]['step_init'])

pow1_model = Pow1Class(x, mu_gauss, CAT, setting["pow1"]['sigma_init'], setting["pow1"]['step_init'], setting["pow1"]['p_init'], setting["pow1"]['p_low'], setting["pow1"]['p_high'], setting["pow1"]['di_gauss'])
pow2_model = Pow2Class(x, mu_gauss, CAT, setting["pow2"]['sigma_init'], setting["pow2"]['step_init'], setting["pow2"]['p1_init'], setting["pow2"]['p1_low'], setting["pow2"]['p1_high'], setting["pow2"]['p2_init'], setting["pow2"]['p2_low'], setting["pow2"]['p2_high'], setting["pow2"]['f1_init'], setting["pow2"]['f2_init'], setting["pow2"]['xmax'], setting["pow2"]['const_f1'], setting["pow2"]['di_gauss'])
pow3_model = Pow3Class(x, mu_gauss, CAT,setting["pow3"]['sigma_init'], setting["pow3"]['step_init'], setting["pow3"]['p1_init'], setting["pow3"]['p1_low'], setting["pow3"]['p1_high'], setting["pow3"]['p2_init'], setting["pow3"]['p2_low'], setting["pow3"]['p2_high'], setting["pow3"]['p3_init'], setting["pow3"]['p3_low'], setting["pow3"]['p3_high'], setting["pow3"]['f1_init'], setting["pow3"]['f2_init'], setting["pow3"]['f3_init'], setting["pow3"]['xmax'], setting["pow3"]['const_f1'], setting["pow3"]['di_gauss'])

exp1_model = Exp1Class(x, mu_gauss, CAT, setting["exp1"]['sigma_init'], setting["exp1"]['step_init'], setting["exp1"]['p_init'], setting["exp1"]['p_low'], setting["exp1"]['p_high'], setting["exp1"]['di_gauss'])
exp2_model = Exp2Class(x, mu_gauss, CAT, setting["exp2"]['sigma_init'], setting["exp2"]['step_init'], setting["exp2"]['p1_init'], setting["exp2"]['p1_low'], setting["exp2"]['p1_high'], setting["exp2"]['p2_init'], setting["exp2"]['p2_low'], setting["exp2"]['p2_high'], setting["exp2"]['f1_init'], setting["exp2"]['f2_init'], setting["exp2"]['xmax'], setting["exp2"]['const_f1'], setting["exp2"]['di_gauss'])
exp3_model = Exp3Class(x, mu_gauss, CAT, setting["exp3"]['sigma_init'], setting["exp3"]['step_init'], setting["exp3"]['p1_init'], setting["exp3"]['p1_low'], setting["exp3"]['p1_high'], setting["exp3"]['p2_init'], setting["exp3"]['p2_low'], setting["exp3"]['p2_high'], setting["exp3"]['p3_init'], setting["exp3"]['p3_low'], setting["exp3"]['p3_high'], setting["exp3"]['f1_init'], setting["exp3"]['f2_init'], setting["exp3"]['f3_init'], setting["exp3"]['xmax'], setting["exp3"]['const_f1'], setting["exp3"]['di_gauss'])

lau2_model = Lau2Class(x, mu_gauss, CAT, setting["lau2"]['sigma_init'], setting["lau2"]['step_init'], setting["lau2"]['p1'], setting["lau2"]['p2'], setting["lau2"]['f_init'], setting["lau2"]['xmax'], setting["lau2"]['const_f1'], setting["lau2"]['di_gauss'])
lau3_model = Lau3Class(x, mu_gauss, CAT, setting["lau3"]['sigma_init'], setting["lau3"]['step_init'], setting["lau3"]['p1'], setting["lau3"]['p2'], setting["lau3"]['p3'], setting["lau3"]['f_init'], setting["lau3"]['xmax'], setting["lau3"]['const_f1'], setting["lau3"]['di_gauss'])
lau4_model = Lau4Class(x, mu_gauss, CAT, setting["lau4"]['sigma_init'], setting["lau4"]['step_init'], setting["lau4"]['p1'], setting["lau4"]['p2'], setting["lau4"]['p3'], setting["lau4"]['p4'], setting["lau4"]['f_init'], setting["lau4"]['xmax'], setting["lau4"]['const_f1'], setting["lau4"]['di_gauss'])

modg_model = ModGausClass(x, CAT, lowx, lowx+65, setting["modg"]['m0'], setting["modg"]['sl'], setting["modg"]['sh'], setting["modg"]['vl'], setting["modg"]['vr'])


bkg_list = []
chi2_pass = setting["Chi2"]

if "bern2" in chi2_pass: bkg_list.append(bern2_model)
if "bern3" in chi2_pass: bkg_list.append(bern3_model)
if "bern4" in chi2_pass: bkg_list.append(bern4_model)
if "bern5" in chi2_pass: bkg_list.append(bern5_model)
if "pow1" in chi2_pass: bkg_list.append(pow1_model)
if "pow2" in chi2_pass: bkg_list.append(pow2_model)
if "pow3" in chi2_pass: bkg_list.append(pow3_model)
if "exp1" in chi2_pass: bkg_list.append(exp1_model)
if "exp2" in chi2_pass: bkg_list.append(exp2_model)
if "exp3" in chi2_pass: bkg_list.append(exp3_model)
if "lau2" in chi2_pass: bkg_list.append(lau2_model)
if "lau3" in chi2_pass: bkg_list.append(lau3_model)
if "lau4" in chi2_pass: bkg_list.append(lau4_model)
if "modg" in chi2_pass: bkg_list.append(modg_model)


dscb_model = DSCB_Class(x, MH, CAT)
# Signal fit
x.setRange("signal",110, 132)
res = dscb_model.pdf.fitTo(hist_sig, ROOT.RooFit.SumW2Error(True),\
                         ROOT.RooFit.Save(True),ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0), ROOT.RooFit.Range("signal"))
res.Print("V")
#plotClass(x, hist_sig, dscb_model.pdf, bkg_model.SBpdf, "Signal_"+cat)
dscb_model.setConst(True)
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
    c1 = ROOT.RooRealVar("c1", "c1", hist_bkg.sumEntries(), 0, 3.* hist_bkg.sumEntries())
    c2 = ROOT.RooRealVar("c2", "c2", 0., -1000., 1000)

    tot_model = ROOT.RooAddPdf(bkg_model.name + "_tot", bkg_model.name + "_tot", ROOT.RooArgList(dscb_model.pdf, bkg_model.pdf), ROOT.RooArgList(c2, c1))
    chi2_tot = ROOT.RooChi2Var("chi2_tot", "chi2_tot", tot_model, hist_bkg, ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson), ROOT.RooFit.Extended(True))
    Minimizer_Chi2(chi2_tot, -1, 100, False, 0)
    res1 = Minimizer_Chi2(chi2_tot, -1, 0.1, False, 0)
    bkg_model.checkBond()
    print("Poisson error:")
    res1.Print("V")
    NS = c2.getVal()
    delta_NS = c2.getError()

    chi2_tot = ROOT.RooChi2Var("chi2_tot", "chi2_tot", tot_model, hist_bkg, ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2), ROOT.RooFit.Extended(True))
    Minimizer_Chi2(chi2_tot, -1, 100, False, 0)
    res2 = Minimizer_Chi2(chi2_tot, -1, 0.1, False, 0)
    bkg_model.checkBond()
    print("SumW2 error:")
    res2.Print("V")
    delta_MC = c2.getError()
    pV = ROOT.Math.chisquared_cdf_c(chi2_tot.getVal(), 260 - res2.floatParsFinal().getSize())
    print("p-value = ", pV)
    print("N signal = ", NS, ", Delta Ns = ", delta_NS, ", Delta MC = ", delta_MC)

    plotClass(x, hist_bkg, tot_model, bkg_model.pdf, "S+B_"+tot_model.GetName(), note='#splitline{N_{S} = ' + '%.2f' % NS + ' #zeta_{S} = ' + '%.2f' % max(0., abs(NS) - 2*delta_MC) + '}{#splitline{#Delta_{MC} = '+ '%.2f' % delta_MC + '}{#Delta_{N_{S}} = '+ '%.2f' % delta_NS+'}}', CMS = 'Simulation')






