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
ROOT.gInterpreter.AddIncludePath('../Utilities/HZGRooPdfs.h')
ROOT.gSystem.Load('../Utilities/HZGRooPdfs_cxx.so')

parser = argparse.ArgumentParser(description = "F test bin and function class")
#parser.add_argument("method")
parser.add_argument('-c', '--cat', help="category")
args = parser.parse_args()
jfile = open('../Config/config.json', 'r')
configs = json.load(jfile)
CAT = args.cat
setting = configs[CAT]
lowx = setting["Range"]

#if not(args.method == "Chi2" or args.method == "NLL") :
#   print("Please use the correct method.")
#sys.exit(1)


# def singleBernFTest(x, gauss_mu, histogram, cat = "", method = "Chi2", e_type = "Poisson", offset = False):
#     bern2_model = Bern2Minization(x, gauss_mu, histogram, method, e_type, cat, 10, 0.3, 10, 7., 105., \
#                                   printLevel = -1, eps = 1, offSet = offset, strategy = 0)
#     bern2_model[0].Print("v")
#     bern3_model = Bern3Minization(x, gauss_mu, histogram, method, e_type, cat, 10, 0.3, 50, 7., 105.,\
#                                   printLevel = -1, eps = 1, offSet = offset, strategy = 0)
#     bern3_model[0].Print("v")
#     bern4_model = Bern4Minization(x, gauss_mu, histogram, method, e_type, cat, 10, 0.3, 50, 7., 105.,\
#                                   printLevel = -1, eps = 1, offSet = offset, strategy = 0)
#     bern4_model[0].Print("v")
#     bern5_model = Bern5Minization(x, gauss_mu, histogram, method, e_type, cat, 10, 0.3, 50, 7., 105.,\
#                                   printLevel = -1, eps = 1, offSet = offset, strategy = 0)
#     bern5_model[0].Print("v")

    
#     stat = [bern2_model[1], bern3_model[1], bern4_model[1], bern5_model[1]]
#     fs = []
#     if method == "Chi2":
#         for i in range(len(stat) - 1):
#             fs.append(ROOT.Math.fdistribution_cdf_c((stat[i] - stat[i+1])*(260-5-i)/stat[i+1], 1, 260-5-i))
#     elif method == "NLL":
#         for i in range(len(stat) - 1):
#             fs.append(ROOT.Math.chisquared_cdf_c(2*(stat[i] - stat[i+1]), 1))

#     print(method, " = ", stat)
#     print("P-value = ", fs)

def goodness(pdfClass, histogram,  e_type = "Poisson", eps = 0.1, n_bins = 260, className = "Default"):
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

    
'''
def singleBernFTest(x, gauss_mu, histogram, cat = "", method = "Chi2", e_type = "Poisson", eps = 0.1, offset = False, strategy = 0, range_ = "", n_bins = 260):
    bern2_model = Bern2Class(x, gauss_mu, cat, 10, 0.3, 10, 7., 105.)
    bern3_model = Bern3Class(x, gauss_mu, cat, 10, 0.3, 50, 7., 105.)
    bern4_model = Bern4Class(x, gauss_mu, cat, 10, 0.3, 50, 7., 105.)
    bern5_model = Bern5Class(x, gauss_mu, cat, 10, 0.3, 50, 7., 105.)
    if e_type == "Poisson": error = ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson)
    elif e_type == "SumW2": error = ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2)
    #cuthistogram = histogram.reduce(ROOT.RooFit.CutRange("left"))
    if method == "Chi2": 
        stat1 = ROOT.RooChi2Var("stat_bern2_" + cat, "stat bern2 " + cat, bern2_model.pdf,  histogram, error)
        stat2 = ROOT.RooChi2Var("stat_bern3_" + cat, "stat bern3 " + cat, bern3_model.pdf,  histogram, error)
        stat3 = ROOT.RooChi2Var("stat_bern4_" + cat, "stat bern4 " + cat, bern4_model.pdf,  histogram, error)
        stat4 = ROOT.RooChi2Var("stat_bern5_" + cat, "stat bern5 " + cat, bern5_model.pdf,  histogram, error)
    elif method == "NLL":
        stat1 = ROOT.RooNLLVar("stat_bern2_" + cat, "stat bern2 " + cat, bern2_model.pdf,  histogram, ROOT.RooFit.Range(range_))
        stat2 = ROOT.RooNLLVar("stat_bern3_" + cat, "stat bern3 " + cat, bern3_model.pdf,  histogram, ROOT.RooFit.Range(range_))
        stat3 = ROOT.RooNLLVar("stat_bern4_" + cat, "stat bern4 " + cat, bern4_model.pdf,  histogram, ROOT.RooFit.Range(range_))
        stat4 = ROOT.RooNLLVar("stat_bern5_" + cat, "stat bern5 " + cat, bern5_model.pdf,  histogram, ROOT.RooFit.Range(range_))

    stats = [stat1, stat2, stat3, stat4]
    if method == "Chi2": 
        for entry in stats:
            print(entry.GetTitle())
            Minimizer_Chi2(entry, -1, 100, False, strategy)
            r = Minimizer_Chi2(entry, -1, eps, offset, strategy)
            r.Print("V")
            # print("Cov q = ", r.covQual(), " status = ", r.status(), end="\n\n")
    elif method == "NLL": 
        for entry in stats:
            print(entry.GetTitle())
            Minimizer_NLL(entry, -1, 100, False, strategy)
            r = Minimizer_NLL(entry, -1, eps, offset, strategy)
            r.Print("V")
    output = [stat1.getVal(), stat2.getVal(), stat3.getVal(), stat4.getVal()]
    #res = bern4_model.pdf.chi2FitTo(cuthistogram, ROOT.RooFit.Range("left"), error, ROOT.RooFit.PrintLevel(-1), ROOT.RooFit.Save(True))
    #res.Print("V")
    fs = []
    bern2_model.checkBond()
    bern3_model.checkBond()
    bern4_model.checkBond()
    bern5_model.checkBond()
    if method == "Chi2":
        for i in range(len(output) - 1):
            fs.append(ROOT.Math.fdistribution_cdf_c((output[i] - output[i+1])*(260-5-i)/output[i+1], 1, 260-5-i))
    elif method == "NLL":
        for i in range(len(output) - 1):
            fs.append(ROOT.Math.chisquared_cdf_c(2*(output[i] - output[i+1]), 1))

    print(method, " = ", output)
    if method == "Chi2": print("goodness = ", [ROOT.Math.chisquared_cdf_c(output[i], n_bins - 4 -i)for i in range(len(output))])
    print("P-value = ", fs)
    plotClass(x, histogram, bern2_model.pdf, title="Bern2 " + cat, output_dir="plots/", sideBand = True, fitRange = range_)
    plotClass(x, histogram, bern3_model.pdf, title="Bern3 " + cat, output_dir="plots/", sideBand = True, fitRange = range_)
    plotClass(x, histogram, bern4_model.pdf, title="Bern4 " + cat, output_dir="plots/", sideBand = True, fitRange = range_)
    plotClass(x, histogram, bern5_model.pdf, title="Bern5 " + cat, output_dir="plots/", sideBand = True, fitRange = range_)

    multiPlotClass(x, histogram, [bern2_model, bern3_model, bern4_model, bern5_model], title="bern multi " + cat, output_dir="plots/", sideBand=True, fitRange= range_)
'''
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
    multiPlotClass(x, histogram, pdfList, title=className+"_multi_" + cat, output_dir="plots/", sideBand=sideBand, fitRange= range_, best_index = highest)

# def customNLL(x, histogram, modelClass):
#     x_list = ROOT.RooArgList("x_list")
#     nx_list = 
#     nll_list = ROOT.RooArgList("nll_list")

# ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)

# Define variables
x = ROOT.RooRealVar("CMS_hzg_mass", "CMS_hzg_mass", lowx, lowx + 65.)
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
x.setBins(260)
x.setRange('left', lowx, 120)
x.setRange('right', 130, lowx+65)
x.setRange('full', lowx, lowx+65)
#reader = readDat(list, "/afs/cern.ch/user/f/fanx/public/samples/")

#Zebing core func hist
#ZBreader = readWsp(x, '/afs/cern.ch/user/f/fanx/EOS_space/zebing_sample/HZGamma_data_bkg_workspace_cat2.root', 'data_mass_cat0')

if CAT=='ggf1':
    read_data = ROOT.RooDataSet.read('../Data/data_ggF1_pinnacles_fix.dat', ROOT.RooArgList(x, y, bdt, w, w_year, year, lep, ph_eta, nlep, njet))
    data = ROOT.RooDataSet('data', 'data', read_data, ROOT.RooArgList(x, y, bdt, w, w_year, year, lep, ph_eta, nlep, njet),'')
    data_hist = ROOT.RooDataHist('hist_data','hist_data', x, data)
elif CAT=='ggf2':
    read_data = ROOT.RooDataSet.read('../Data/data_ggF2_pinnacles_fix.dat', ROOT.RooArgList(x, y, bdt, w, w_year, year, lep, ph_eta, nlep, njet))
    data = ROOT.RooDataSet('data', 'data', read_data, ROOT.RooArgList(x, y, bdt, w, w_year, year, lep, ph_eta, nlep, njet),'')
    data_hist = ROOT.RooDataHist('hist_data','hist_data', x, data)
elif CAT=='ggf3':
    read_data = ROOT.RooDataSet.read('../Data/data_ggF3_pinnacles_fix.dat', ROOT.RooArgList(x, y, bdt, w, w_year, year, lep, ph_eta, nlep, njet))
    data = ROOT.RooDataSet('data', 'data', read_data, ROOT.RooArgList(x, y, bdt, w, w_year, year, lep, ph_eta, nlep, njet),'')
    data_hist = ROOT.RooDataHist('hist_data','hist_data', x, data)
elif CAT=='ggf4':
    read_data = ROOT.RooDataSet.read('../Data/data_ggF4_pinnacles_fix.dat', ROOT.RooArgList(x, y, bdt, w, w_year, year, lep, ph_eta, nlep, njet))
    data = ROOT.RooDataSet('data', 'data', read_data, ROOT.RooArgList(x, y, bdt, w, w_year, year, lep, ph_eta, nlep, njet),'')
    data_hist = ROOT.RooDataHist('hist_data','hist_data', x, data)

# Define PDF classes
bern2_model = Bern2Class(x, mu_gauss, CAT, setting["bern2"]['p0'], setting["bern2"]['p_init'],setting["bern2"]['bond'],setting["bern2"]['sigma_init'],setting["bern2"]['step_init'])
bern3_model = Bern3Class(x, mu_gauss, CAT, setting["bern3"]['p0'], setting["bern3"]['p_init'],setting["bern3"]['bond'],setting["bern3"]['sigma_init'],setting["bern3"]['step_init'])
bern4_model = Bern4Class(x, mu_gauss, CAT, setting["bern4"]['p0'], setting["bern4"]['p_init'],setting["bern4"]['bond'],setting["bern4"]['sigma_init'],setting["bern4"]['step_init'])
bern5_model = Bern5Class(x, mu_gauss, CAT, setting["bern5"]['p0'], setting["bern5"]['p_init'],setting["bern5"]['bond'],setting["bern5"]['sigma_init'],setting["bern5"]['step_init'])
bern_list = []

core_bern2_model = CoreBern2Class(x, CAT, 1, 0.3, 10, fileName="ZGCoreShape_01jet_NAFCorr", shapeName="CoreShape_ZG_NAF_cat2")
core_bern3_model = CoreBern3Class(x, CAT, 1, 0.3, 10, fileName="ZGCoreShape_01jet_NAFCorr", shapeName="CoreShape_ZG_NAF_cat2")
core_bern4_model = CoreBern4Class(x, CAT, 1, 0.3, 10, fileName="ZGCoreShape_01jet_NAFCorr", shapeName="CoreShape_ZG_NAF_cat2")
core_bern5_model = CoreBern5Class(x, CAT, 1, 0.3, 10, fileName="ZGCoreShape_01jet_NAFCorr", shapeName="CoreShape_ZG_NAF_cat2")
core_bern_list = [core_bern2_model, core_bern3_model, core_bern4_model, core_bern5_model]

pow1_model = Pow1Class(x, mu_gauss, CAT, setting["pow1"]['sigma_init'], setting["pow1"]['step_init'], setting["pow1"]['p_init'], setting["pow1"]['p_low'], setting["pow1"]['p_high'], setting["pow1"]['di_gauss'])
pow2_model = Pow2Class(x, mu_gauss, CAT, setting["pow2"]['sigma_init'], setting["pow2"]['step_init'], setting["pow2"]['p1_init'], setting["pow2"]['p1_low'], setting["pow2"]['p1_high'], setting["pow2"]['p2_init'], setting["pow2"]['p2_low'], setting["pow2"]['p2_high'], setting["pow2"]['f1_init'], setting["pow2"]['f2_init'], setting["pow2"]['xmax'], setting["pow2"]['const_f1'], setting["pow2"]['di_gauss'])
pow3_model = Pow3Class(x, mu_gauss, CAT,setting["pow3"]['sigma_init'], setting["pow3"]['step_init'], setting["pow3"]['p1_init'], setting["pow3"]['p1_low'], setting["pow3"]['p1_high'], setting["pow3"]['p2_init'], setting["pow3"]['p2_low'], setting["pow3"]['p2_high'], setting["pow3"]['p3_init'], setting["pow3"]['p3_low'], setting["pow3"]['p3_high'], setting["pow3"]['f1_init'], setting["pow3"]['f2_init'], setting["pow3"]['f3_init'], setting["pow3"]['xmax'], setting["pow3"]['const_f1'], setting["pow3"]['di_gauss'])
pow_list = []

core_pow1_model = CorePow1Class(x, CAT, -6, fileName="ZGCoreShape_01jet_NAFCorr", shapeName="CoreShape_ZG_NAF_cat2")
core_pow2_model = CorePow2Class(x, CAT, -6, -15, 0.5,  fileName="ZGCoreShape_01jet_NAFCorr", shapeName="CoreShape_ZG_NAF_cat2")
core_pow_list = [core_pow1_model, core_pow2_model]

exp1_model = Exp1Class(x, mu_gauss, CAT, setting["exp1"]['sigma_init'], setting["exp1"]['step_init'], setting["exp1"]['p_init'], setting["exp1"]['p_low'], setting["exp1"]['p_high'], setting["exp1"]['di_gauss'])
exp2_model = Exp2Class(x, mu_gauss, CAT, setting["exp2"]['sigma_init'], setting["exp2"]['step_init'], setting["exp2"]['p1_init'], setting["exp2"]['p1_low'], setting["exp2"]['p1_high'], setting["exp2"]['p2_init'], setting["exp2"]['p2_low'], setting["exp2"]['p2_high'], setting["exp2"]['f1_init'], setting["exp2"]['f2_init'], setting["exp2"]['xmax'], setting["exp2"]['const_f1'], setting["exp2"]['di_gauss'])
exp3_model = Exp3Class(x, mu_gauss, CAT, setting["exp3"]['sigma_init'], setting["exp3"]['step_init'], setting["exp3"]['p1_init'], setting["exp3"]['p1_low'], setting["exp3"]['p1_high'], setting["exp3"]['p2_init'], setting["exp3"]['p2_low'], setting["exp3"]['p2_high'], setting["exp3"]['p3_init'], setting["exp3"]['p3_low'], setting["exp3"]['p3_high'], setting["exp3"]['f1_init'], setting["exp3"]['f2_init'], setting["exp3"]['f3_init'], setting["exp3"]['xmax'], setting["exp3"]['const_f1'], setting["exp3"]['di_gauss'])
exp_list = []

lau2_model = Lau2Class(x, mu_gauss, CAT, setting["lau2"]['sigma_init'], setting["lau2"]['step_init'], setting["lau2"]['p1'], setting["lau2"]['p2'], setting["lau2"]['f_init'], setting["lau2"]['xmax'], setting["lau2"]['const_f1'], setting["lau2"]['di_gauss'])
lau3_model = Lau3Class(x, mu_gauss, CAT, setting["lau3"]['sigma_init'], setting["lau3"]['step_init'], setting["lau3"]['p1'], setting["lau3"]['p2'], setting["lau3"]['p3'], setting["lau3"]['f_init'], setting["lau3"]['xmax'], setting["lau3"]['const_f1'], setting["lau3"]['di_gauss'])
lau4_model = Lau4Class(x, mu_gauss, CAT, setting["lau4"]['sigma_init'], setting["lau4"]['step_init'], setting["lau4"]['p1'], setting["lau4"]['p2'], setting["lau4"]['p3'], setting["lau4"]['p4'], setting["lau4"]['f_init'], setting["lau4"]['xmax'], setting["lau4"]['const_f1'], setting["lau4"]['di_gauss'])
lau_list = []

modg_model = ModGausClass(x, CAT, lowx, lowx+65, setting["modg"]['m0'], setting["modg"]['sl'], setting["modg"]['sh'], setting["modg"]['vl'], setting["modg"]['vr'])

if "bern2" in setting["SST"]:
     bern_list.append(bern2_model)
if "bern3" in setting["SST"]:
     bern_list.append(bern3_model)
if "bern4" in setting["SST"]:
     bern_list.append(bern4_model)
if "bern5" in setting["SST"]:
     bern_list.append(bern5_model)
if "pow1" in setting["SST"]:
     pow_list.append(pow1_model)
if "pow2" in setting["SST"]:
     pow_list.append(pow2_model)
if "pow3" in setting["SST"]:
     pow_list.append(pow3_model)
if "exp1" in setting["SST"]:
     exp_list.append(exp1_model)
if "exp2" in setting["SST"]:
     exp_list.append(exp2_model)
if "exp3" in setting["SST"]:
     exp_list.append(exp3_model)
if "lau2" in setting["SST"]:
     lau_list.append(lau2_model)
if "lau3" in setting["SST"]:
     lau_list.append(lau3_model)
if "lau4" in setting["SST"]:
     lau_list.append(lau4_model)


    
best_list = [bern3_model, bern4_model, pow1_model, pow2_model, exp2_model, modg_model]
# Goodness of fit test
# goodness(bern_list, reader.data_hist_untagged2_bkg,  e_type = "Poisson", eps = 0.1, n_bins = 260, className="Bern")
#goodness(pow_list, reader.data_hist_untagged2_bkg,  e_type = "SumW2", eps = 0.1, n_bins = 260, className="Bern")
#goodness(exp_list, reader.data_hist_untagged1_bkg,  e_type = "SumW2", eps = 0.1, n_bins = 260, className="Exp")
# goodness(lau_list, reader.data_hist_untagged1_bkg,  e_type = "Poisson", eps = 0.1, n_bins = 260, className="Bern")

# F Tset
# singleBernFTest(x, mu_gauss, reader.data_u2, CAT, args.method, "Poisson", eps = 0.1, offset = True, strategy = 0, range_ = "left,right", n_bins = 220)
# singleFTestSidebandNLL(x, pow_list, reader.data_u2, cat = CAT, eps = 0.1, offset = True, strategy = 0, range_= "left,right", className = "Pow")
#singleFTestSidebandNLL(x, exp_list, reader.data_u1, cat = CAT, eps = 0.1, offset = True, strategy = 0, range_= "left,right", className = "Exp")
#singleFTestSidebandNLL(x, lau_list, reader.data_u1, cat = CAT, eps = 0.1, offset = True, strategy = 0, range_= "left,right", calssName = "Lau")
if len(bern_list) > 1:
    singleFTestSidebandNLL(x, bern_list, data_hist, cat = CAT, eps = 0.1, offset = True, strategy = 0, range_= "left,right", className = "Bern", sideBand = True)
if len(pow_list) > 1:
    singleFTestSidebandNLL(x, pow_list, data_hist, cat = CAT, eps = 0.01, offset = True, strategy = 0, range_= "left,right", className = "Pow")
if len(exp_list) > 1:
    singleFTestSidebandNLL(x, exp_list, data_hist, cat = CAT, eps = 0.01, offset = True, strategy = 0, range_= "left,right", className = "Exp")
if len(lau_list) > 1:
    singleFTestSidebandNLL(x, lau_list, data_hist, cat = CAT, eps = 0.1, offset = True, strategy = 0, range_= "left,right", className = "Lau")


# can2 = ROOT.TCanvas("c2","c2", 500, 500)
# can2.cd()
# plot2 = x.frame()
# # x.setBins(65)
# show_hist_data = data_hist.createHistogram("h_hist", x, ROOT.RooFit.Binning(65))
# show_hist_mc = mc_hist.createHistogram("h_hist_mc", x, ROOT.RooFit.Binning(65))

# show_hist_data.Scale(1./show_hist_data.Integral())
# show_hist_mc.Scale(1./show_hist_mc.Integral())
# show_hist_data.SetLineColor(2)
# show_hist_mc.SetLineColor(3)
# show_hist_data.Draw("HIST")
# show_hist_mc.Draw("SAME HIST")
# can2.SaveAs("test.pdf")




