#!/usr/bin/env python3
import ROOT
import os
import sys
import json
import argparse
sys.path.append(os.path.abspath("../Utilities/"))
sys.path.append(os.path.abspath("../CMS_plotter/"))
import CMS_lumi, tdrstyle
#from bkg_functions_fit import *
from bkg_functions_class import *
from Xc_Minimizer import *
from plot_utility import *
from sample_reader import *
from profile_class import *
#ROOT.gInterpreter.AddIncludePath('../Utilities/HZGRooPdfs.h')
#ROOT.gSystem.Load('../Utilities/HZGRooPdfs_cxx.so')
ROOT.gInterpreter.AddIncludePath('../Utilities/RooGaussStepBernstein.h')
ROOT.gSystem.Load('../Utilities/RooGaussStepBernstein_cxx.so')

parser = argparse.ArgumentParser(description = "F test bin and function class")
#parser.add_argument("method")
parser.add_argument('-c', '--cat', help="category")
parser.add_argument('-con', '--config', help = 'Configuration')
args = parser.parse_args()
jfile = open(args.config, 'r')
configs = json.load(jfile)
CAT = args.cat
setting = configs[CAT]
lowx = setting["Range"][0]
highx = setting["Range"][1]
nbins = int(setting["Bins"])

def goodness(pdfClass, histogram,  e_type = "Poisson", eps = 0.1, n_bins = nbins, className = "Default"):
    if e_type == "Poisson": error = ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson)
    elif e_type == "SumW2": error = ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2)
    good_PV = []
    chi2_=[]
    for entry in pdfClass:
        stat =  ROOT.RooChi2Var("stat_goodness", "goodness test", entry.pdf,  histogram, error)
        Minimizer_Chi2(stat, -1, 100, False, 0)
        r = Minimizer_Chi2(stat, -1, eps, False, 0)
        if r.status() != 0: 
                print(entry.pdf.GetName(), " Minimization fails!")
        entry.checkBond()
        chi2_.append(stat.getVal())
        good_PV.append(ROOT.Math.chisquared_cdf_c(stat.getVal(), n_bins - r.floatParsFinal().getSize()))
    print(className, " goodness = ", good_PV)
    print(className, " Chi2 = ", chi2_)

    

def singleFTestSidebandNLL(x, pdfList, histogram, cat = '', eps = 0.1, offset = True, strategy = 0, range_= "", className = "Default", sideBand = True, FT = True):
    stats = []
    fitres = []
    if sideBand:
        cuthistogram = histogram.reduce(ROOT.RooFit.CutRange(range_))
    else:
        cuthistogram = histogram
    for entry in pdfList:
        entry.reset()
        # nll1_= ROOT.RooNLLVar("stat1_" + entry.pdf.GetName(), "stat1 " + entry.pdf.GetName(), entry.pdf,  histogram, ROOT.RooFit.Range('left'))
        # nll2_= ROOT.RooNLLVar("stat2_" + entry.pdf.GetName(), "stat2 " + entry.pdf.GetName(), entry.pdf,  histogram, ROOT.RooFit.Range('right'))
        if sideBand: nll_= ROOT.RooNLLVar("stat_" + entry.SBpdf.GetName(), "stat " + entry.SBpdf.GetName(), entry.SBpdf,  cuthistogram)
        else: nll_= ROOT.RooNLLVar("stat_" + entry.pdf.GetName(), "stat " + entry.pdf.GetName(), entry.pdf,  cuthistogram)
        #nll2_= ROOT.RooChi2Var("stat2_" + entry.pdf.GetName(), "stat2 " + entry.pdf.GetName(), entry.pdf,  histogram, ROOT.RooFit.Range('right'))

        #nll_ = ROOT.RooAddition("stat", "stat", ROOT.RooArgList(nll1_, nll2_))
        stats.append(nll_)
        Minimizer_NLL(nll_, -1, 100, False, strategy)
        r_=Minimizer_NLL(nll_, -1, eps, offset, strategy)
        r_.Print("V")
        fitres.append(r_)
        entry.checkBond()
        plotClass(x, histogram, entry.pdf, entry.SBpdf, title=entry.pdf.GetName() + "_" + cat, output_dir="plots/", sideBand = sideBand, fitRange = range_)

    print("NLL = ", [ele.getVal() for ele in stats])
    fs = []
    dif = []
    corrNLL = []
    if FT:
        for i in range(len(stats) - 1):
            fs.append(ROOT.Math.chisquared_cdf_c(2*(stats[i].getVal() - stats[i+1].getVal()), fitres[i+1].floatParsFinal().getSize() - fitres[i].floatParsFinal().getSize()))
            dif.append(fitres[i+1].floatParsFinal().getSize() - fitres[i].floatParsFinal().getSize())
    print("DOF diff = ", dif)
    print("P-value = ", fs)
    for i in range(len(stats)):
        corrNLL.append(stats[i].getVal()+ 0.5 * fitres[i].floatParsFinal().getSize())
    print("corr NLL = ", corrNLL)
    highest = 0
    for i, f in enumerate(fs):
        if f<0.05:
            highest = i+1
        else: break
    multiPlotClass(x, histogram, pdfList, title=className+"_multi_" + cat, output_dir="plots/", sideBand=sideBand, fitRange= range_, best_index = highest, bestLabel = True)

# def customNLL(x, histogram, modelClass):
#     x_list = ROOT.RooArgList("x_list")
#     nx_list = 
#     nll_list = ROOT.RooArgList("nll_list")

# ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)

# Define variables
x = ROOT.RooRealVar("CMS_hzg_mass", "CMS_hzg_mass", lowx, highx)
y = ROOT.RooRealVar("y", "photon pT", 15., 1000.)
w = ROOT.RooRealVar("w", "w", -40., 40.)
bdt = ROOT.RooRealVar("bdt", "bdt", -1, 1)
year = ROOT.RooRealVar("year", "year", 2015, 21000)
w_year = ROOT.RooRealVar("w_year", "w_year", 0, 100)
lep = ROOT.RooRealVar("lep", "lep", 0, 1) #0 = electron, 1 = muon
ph_eta = ROOT.RooRealVar("ph_eta", "ph_eta", -3, 3)
nlep = ROOT.RooRealVar("nlep", "nlep", 0, 10)
njet = ROOT.RooRealVar("njet", "njet", 0, 10)

mu_gauss = ROOT.RooRealVar("mu_gauss","always 0"       ,0.)
list = [x, y, w, bdt, year, lep, ph_eta, nlep, njet]

# Expected bkg
# expbkg_u1 = u1_bkg_run2.sumEntries()
# expbkg_u2 = u2_bkg_run2.sumEntries()
# expbkg_u3 = u3_bkg_run2.sumEntries()
# expbkg_u4 =  u4_bkg_run2.sumEntries()

# Cornell MC sample dat reader and make RooDataHist
x.setBins(nbins)
x.setRange('left', lowx, 120)
x.setRange('right', 130, highx)
x.setRange('full', lowx, highx)
#reader = readDat(list, "/afs/cern.ch/user/f/fanx/public/samples/")

#Zebing core func hist
#ZBreader = readWsp(x, '/afs/cern.ch/user/f/fanx/EOS_space/zebing_sample/HZGamma_data_bkg_workspace_cat2.root', 'data_mass_cat0')
DAT = False
if DAT:
    if CAT=='ggf1':
        read_data = ROOT.RooDataSet.read('../Data/data_ggF1_pinnacles_fix.dat', ROOT.RooArgList(x, y, bdt, w, w_year, year, lep, ph_eta, nlep, njet))
        data = ROOT.RooDataSet('data', 'data', read_data, ROOT.RooArgList(x, y, bdt, w, w_year, year, lep, ph_eta, nlep, njet),'')
        hist_data = ROOT.RooDataHist('hist_data','hist_data', x, data)
    elif CAT=='ggf2':
        read_data = ROOT.RooDataSet.read('../Data/data_ggF2_pinnacles_fix.dat', ROOT.RooArgList(x, y, bdt, w, w_year, year, lep, ph_eta, nlep, njet))
        data = ROOT.RooDataSet('data', 'data', read_data, ROOT.RooArgList(x, y, bdt, w, w_year, year, lep, ph_eta, nlep, njet),'')
        hist_data = ROOT.RooDataHist('hist_data','hist_data', x, data)
    elif CAT=='ggf3':
        read_data = ROOT.RooDataSet.read('../Data/data_ggF3_pinnacles_fix.dat', ROOT.RooArgList(x, y, bdt, w, w_year, year, lep, ph_eta, nlep, njet))
        data = ROOT.RooDataSet('data', 'data', read_data, ROOT.RooArgList(x, y, bdt, w, w_year, year, lep, ph_eta, nlep, njet),'')
        hist_data = ROOT.RooDataHist('hist_data','hist_data', x, data)
    elif CAT=='ggf4':
        read_data = ROOT.RooDataSet.read('../Data/data_ggF4_pinnacles_fix.dat', ROOT.RooArgList(x, y, bdt, w, w_year, year, lep, ph_eta, nlep, njet))
        data = ROOT.RooDataSet('data', 'data', read_data, ROOT.RooArgList(x, y, bdt, w, w_year, year, lep, ph_eta, nlep, njet),'')
        hist_data = ROOT.RooDataHist('hist_data','hist_data', x, data)

else:
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

# Define PDF classes
CORE = False
if CORE:
    core_bern2_model = CoreBern2Class(x, CAT, 1, 0.3, 10, fileName="ZGCoreShape_01jet_NAFCorr", shapeName="CoreShape_ZG_NAF_cat2")
    core_bern3_model = CoreBern3Class(x, CAT, 1, 0.3, 10, fileName="ZGCoreShape_01jet_NAFCorr", shapeName="CoreShape_ZG_NAF_cat2")
    core_bern4_model = CoreBern4Class(x, CAT, 1, 0.3, 10, fileName="ZGCoreShape_01jet_NAFCorr", shapeName="CoreShape_ZG_NAF_cat2")
    core_bern5_model = CoreBern5Class(x, CAT, 1, 0.3, 10, fileName="ZGCoreShape_01jet_NAFCorr", shapeName="CoreShape_ZG_NAF_cat2")
    core_bern_list = [core_bern2_model, core_bern3_model, core_bern4_model, core_bern5_model]

    core_pow1_model = CorePow1Class(x, CAT, -6, fileName="ZGCoreShape_01jet_NAFCorr", shapeName="CoreShape_ZG_NAF_cat2")
    core_pow2_model = CorePow2Class(x, CAT, -6, -15, 0.5,  fileName="ZGCoreShape_01jet_NAFCorr", shapeName="CoreShape_ZG_NAF_cat2")
    core_pow_list = [core_pow1_model, core_pow2_model]

profile = profileClass(x, mu_gauss, CAT, args.config)
bern_list = []
pow_list = []
exp_list = []
lau_list = []

if "bern2" in setting["SST"]:
     bern_list.append(profile.bern2_model)
if "bern3" in setting["SST"]:
     bern_list.append(profile.bern3_model)
if "bern4" in setting["SST"]:
     bern_list.append(profile.bern4_model)
if "bern5" in setting["SST"]:
     bern_list.append(profile.bern5_model)
if "pow1" in setting["SST"]:
     pow_list.append(profile.pow1_model)
if "pow2" in setting["SST"]:
     pow_list.append(profile.pow2_model)
if "pow3" in setting["SST"]:
     pow_list.append(profile.pow3_model)
if "exp1" in setting["SST"]:
     exp_list.append(profile.exp1_model)
if "exp2" in setting["SST"]:
     exp_list.append(profile.exp2_model)
if "exp3" in setting["SST"]:
     exp_list.append(profile.exp3_model)
if "lau2" in setting["SST"]:
     lau_list.append(profile.lau2_model)
if "lau3" in setting["SST"]:
     lau_list.append(profile.lau3_model)
if "lau4" in setting["SST"]:
     lau_list.append(profile.lau4_model)


    
best_list = []
# Goodness of fit test (NOT USED)
# goodness(bern_list, reader.hist_data_untagged2_bkg,  e_type = "Poisson", eps = 0.1, n_bins = nbins, className="Bern")
#goodness(pow_list, reader.hist_data_untagged2_bkg,  e_type = "SumW2", eps = 0.1, n_bins = nbins, className="Bern")
#goodness(exp_list, reader.hist_data_untagged1_bkg,  e_type = "SumW2", eps = 0.1, n_bins = nbins, className="Exp")
# goodness(lau_list, reader.hist_data_untagged1_bkg,  e_type = "Poisson", eps = 0.1, n_bins = nbins, className="Bern")

# F Tset
if len(bern_list) > 1:
    singleFTestSidebandNLL(x, bern_list, hist_data, cat = CAT, eps = 0.1, offset = True, strategy = 0, range_= "left,right", className = "Bern", sideBand = True)
if len(pow_list) > 1:
    singleFTestSidebandNLL(x, pow_list, hist_data, cat = CAT, eps = 0.01, offset = True, strategy = 0, range_= "left,right", className = "Pow")
if len(exp_list) > 1:
    singleFTestSidebandNLL(x, exp_list, hist_data, cat = CAT, eps = 0.01, offset = True, strategy = 0, range_= "left,right", className = "Exp")
if len(lau_list) > 1:
    singleFTestSidebandNLL(x, lau_list, hist_data, cat = CAT, eps = 0.1, offset = True, strategy = 0, range_= "left,right", className = "Lau")

# can2 = ROOT.TCanvas("c2","c2", 500, 500)
# can2.cd()
# plot2 = x.frame()
# # x.setBins(65)
# show_hist_data = hist_data.createHistogram("h_hist", x, ROOT.RooFit.Binning(65))
# show_hist_mc = mc_hist.createHistogram("h_hist_mc", x, ROOT.RooFit.Binning(65))

# show_hist_data.Scale(1./show_hist_data.Integral())
# show_hist_mc.Scale(1./show_hist_mc.Integral())
# show_hist_data.SetLineColor(2)
# show_hist_mc.SetLineColor(3)
# show_hist_data.Draw("HIST")
# show_hist_mc.Draw("SAME HIST")
# can2.SaveAs("test.pdf")




