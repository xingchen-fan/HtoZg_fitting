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

ROOT.gROOT.SetBatch(True)

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
            print("Highest d.o.f =" +str(i+1))
        else: break
    multiPlotClass(x, histogram, pdfList, title=className+"_multi_" + cat, output_dir="plots/", sideBand=sideBand, fitRange= range_, best_index = highest,  sideband_bins=70, ratio_range=[0.5,1.5])

ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)

parser = argparse.ArgumentParser(description = "F test bin and function class")
#parser.add_argument("method")
parser.add_argument('-c', '--cat', help="category")
parser.add_argument('-con', '--config', help = 'Configuration')
parser.add_argument('-l', '--poly_list', help = 'Polynomial list')

args = parser.parse_args()
jfile = open('../Config/'+args.config, 'r')
configs = json.load(jfile)
CAT = args.cat
setting = configs[CAT]
poly_list = args.poly_list

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
bern_list = []
pow_list = []
exp_list = []
lau_list = []

input_key = "FT"

if "bern2" in setting[input_key]:
     bern_list.append(profile.bern2_model)
if "bern3" in setting[input_key]:
     bern_list.append(profile.bern3_model)
if "bern4" in setting[input_key]:
     bern_list.append(profile.bern4_model)
if "bern5" in setting[input_key]:
     bern_list.append(profile.bern5_model)
if "pow1" in setting[input_key]:
     pow_list.append(profile.pow1_model)
if "pow2" in setting[input_key]:
     pow_list.append(profile.pow2_model)
if "pow3" in setting[input_key]:
     pow_list.append(profile.pow3_model)
if "exmg" in setting[input_key]:
     exp_list.append(profile.exmg_model)
if "exp2" in setting[input_key]:
     exp_list.append(profile.exp2_model)
if "exp3" in setting[input_key]:
     exp_list.append(profile.exp3_model)
if "lau2" in setting[input_key]:
     lau_list.append(profile.lau2_model)
if "lau3" in setting[input_key]:
     lau_list.append(profile.lau3_model)
if "lau4" in setting[input_key]:
     lau_list.append(profile.lau4_model)
 
best_list = []
if (len(bern_list) > 1 and poly_list=="BERN"):
    print("BERNSTEIN POLYNOMIAL LIST: ")
    singleFTestSidebandNLL(x, bern_list, data_hist, cat = CAT, eps = 0.1, offset = True, strategy = 0, range_= "left,right", className = "Bern", sideBand = True)
    print("-----------------------------------------------------------------------------------------------------------------------------------------------------")
if (len(pow_list) > 1 and poly_list=="POW"):
    print("POWER POLYNOMIAL LIST: ")
    singleFTestSidebandNLL(x, pow_list, data_hist, cat = CAT, eps = 0.01, offset = True, strategy = 0, range_= "left,right", className = "Pow")
    print("-----------------------------------------------------------------------------------------------------------------------------------------------------")
if (len(exp_list) > 1 and poly_list=="EXP"):
    print("EXP. POLYNOMIAL LIST: ")
    singleFTestSidebandNLL(x, exp_list, data_hist, cat = CAT, eps = 0.01, offset = True, strategy = 0, range_= "left,right", className = "Exp")
    print("-----------------------------------------------------------------------------------------------------------------------------------------------------")
if (len(lau_list) > 1 and poly_list=="LAU"):
    print("LANDAU POLYNOMIAL LIST: ")
    singleFTestSidebandNLL(x, lau_list, data_hist, cat = CAT, eps = 0.1, offset = True, strategy = 0, range_= "left,right", className = "Lau")
    print("-----------------------------------------------------------------------------------------------------------------------------------------------------")

# Define PDF classes
#core_bern2_model = CoreBern2Class(x, CAT, 1, 0.3, 10, fileName="ZGCoreShape_01jet_NAFCorr", shapeName="CoreShape_ZG_NAF_cat2")
#core_bern3_model = CoreBern3Class(x, CAT, 1, 0.3, 10, fileName="ZGCoreShape_01jet_NAFCorr", shapeName="CoreShape_ZG_NAF_cat2")
#core_bern4_model = CoreBern4Class(x, CAT, 1, 0.3, 10, fileName="ZGCoreShape_01jet_NAFCorr", shapeName="CoreShape_ZG_NAF_cat2")
#core_bern5_model = CoreBern5Class(x, CAT, 1, 0.3, 10, fileName="ZGCoreShape_01jet_NAFCorr", shapeName="CoreShape_ZG_NAF_cat2")
#core_bern_list = [core_bern2_model, core_bern3_model, core_bern4_model, core_bern5_model]

#core_pow1_model = CorePow1Class(x, CAT, -6, fileName="ZGCoreShape_01jet_NAFCorr", shapeName="CoreShape_ZG_NAF_cat2")
#core_pow2_model = CorePow2Class(x, CAT, -6, -15, 0.5,  fileName="ZGCoreShape_01jet_NAFCorr", shapeName="CoreShape_ZG_NAF_cat2")
#core_pow_list = [core_pow1_model, core_pow2_model]

# def customNLL(x, histogram, modelClass):
#     x_list = ROOT.RooArgList("x_list")
#     nx_list = 
#     nll_list = ROOT.RooArgList("nll_list")

# ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

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

# Goodness of fit test (NOT USED)
# goodness(bern_list, reader.data_hist_untagged2_bkg,  e_type = "Poisson", eps = 0.1, n_bins = 260, className="Bern")
#goodness(pow_list, reader.data_hist_untagged2_bkg,  e_type = "SumW2", eps = 0.1, n_bins = 260, className="Bern")
#goodness(exp_list, reader.data_hist_untagged1_bkg,  e_type = "SumW2", eps = 0.1, n_bins = 260, className="Exp")
# goodness(lau_list, reader.data_hist_untagged1_bkg,  e_type = "Poisson", eps = 0.1, n_bins = 260, className="Bern")

# F Tset
# singleBernFTest(x, mu_gauss, reader.data_u2, CAT, args.method, "Poisson", eps = 0.1, offset = True, strategy = 0, range_ = "left,right", n_bins = 220)
# singleFTestSidebandNLL(x, pow_list, reader.data_u2, cat = CAT, eps = 0.1, offset = True, strategy = 0, range_= "left,right", className = "Pow")
#singleFTestSidebandNLL(x, exp_list, reader.data_u1, cat = CAT, eps = 0.1, offset = True, strategy = 0, range_= "left,right", className = "Exp")
#singleFTestSidebandNLL(x, lau_list, reader.data_u1, cat = CAT, eps = 0.1, offset = True, strategy = 0, range_= "left,right", calssName = "Lau")




