import ROOT
import os
import sys
import argparse
sys.path.append(os.path.abspath("../Utilities/"))
from bkg_functions_fit import *
from bkg_functions_class import *
from Xc_Minimizer import *
from plot_utility import *
from sample_reader import *

parser = argparse.ArgumentParser(description = "F test method (NLL or Chi2)")
parser.add_argument("method")
args = parser.parse_args()
if not(args.method == "Chi2" or args.method == "NLL") :
    print("Please use the correct method.")
    sys.exit(1)


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

def singleBernFTest(x, gauss_mu, histogram, cat = "", method = "Chi2", e_type = "Poisson", eps = 0.1, offset = False, strategy = 0):
    bern2_model = Bern2Class(x, gauss_mu, cat, 10, 0.3, 10, 7., 105.)
    bern3_model = Bern3Class(x, gauss_mu, cat, 10, 0.3, 50, 7., 105.)
    bern4_model = Bern4Class(x, gauss_mu, cat, 10, 0.3, 50, 7., 105.)
    bern5_model = Bern5Class(x, gauss_mu, cat, 10, 0.3, 50, 7., 105.)
    if e_type == "Poisson": error = ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson)
    elif e_type == "SumW2": error = ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2)

    if method == "Chi2": 
        stat1 = ROOT.RooChi2Var("stat_bern2_" + cat, "stat bern2 " + cat, bern2_model.pdf,  histogram, error)
        stat2 = ROOT.RooChi2Var("stat_bern3_" + cat, "stat bern3 " + cat, bern3_model.pdf,  histogram, error)
        stat3 = ROOT.RooChi2Var("stat_bern4_" + cat, "stat bern4 " + cat, bern4_model.pdf,  histogram, error)
        stat4 = ROOT.RooChi2Var("stat_bern5_" + cat, "stat bern5 " + cat, bern5_model.pdf,  histogram, error)
    elif method == "NLL":
        stat1 = ROOT.RooNLLVar("stat_bern2_" + cat, "stat bern2 " + cat, bern2_model.pdf,  histogram)
        stat2 = ROOT.RooNLLVar("stat_bern3_" + cat, "stat bern3 " + cat, bern3_model.pdf,  histogram)
        stat3 = ROOT.RooNLLVar("stat_bern4_" + cat, "stat bern4 " + cat, bern4_model.pdf,  histogram)
        stat4 = ROOT.RooNLLVar("stat_bern5_" + cat, "stat bern5 " + cat, bern5_model.pdf,  histogram)

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
            r = Minimizer_NLL(entry, -1, eps, offset, strategy).Print("V")
            r.Print("V")
    output = [stat1.getVal(), stat2.getVal(), stat3.getVal(), stat4.getVal()]
    fs = []
    if method == "Chi2":
        for i in range(len(output) - 1):
            fs.append(ROOT.Math.fdistribution_cdf_c((output[i] - output[i+1])*(260-5-i)/output[i+1], 1, 260-5-i))
    elif method == "NLL":
        for i in range(len(output) - 1):
            fs.append(ROOT.Math.chisquared_cdf_c(2*(output[i] - output[i+1]), 1))

    print(method, " = ", output)
    print("P-value = ", fs)
    plotClass(x, histogram, bern2_model.pdf)



# ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
lowx = 105.

# Define variables
x = ROOT.RooRealVar("x", "mllg", lowx, lowx + 65.)
y = ROOT.RooRealVar("y", "photon pT", 15., 1000.)
w = ROOT.RooRealVar("w", "w", -40., 40.)
bdt = ROOT.RooRealVar("bdt", "bdt", -1, 1)
year = ROOT.RooRealVar("year", "year", 2015, 2019)
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

# Dat reader
reader = readDat(list)

# Make RooDataHist
x.setBins(260)
data_hist_untagged1 = ROOT.RooDataHist("data_hist_untagged1", "data_hist_untagged1", x, reader.u1_tot_run2)
data_hist_untagged2 = ROOT.RooDataHist("data_hist_untagged2", "data_hist_untagged2", x, reader.u2_tot_run2)
data_hist_untagged3 = ROOT.RooDataHist("data_hist_untagged3", "data_hist_untagged3", x, reader.u3_tot_run2)
data_hist_untagged4 = ROOT.RooDataHist("data_hist_untagged4", "data_hist_untagged4", x, reader.u4_tot_run2)

data_hist_untagged1_bkg = ROOT.RooDataHist("data_hist_untagged1_bkg", "data_hist_untagged1_bkg", x, reader.u1_bkg_run2)
data_hist_untagged2_bkg = ROOT.RooDataHist("data_hist_untagged2_bkg", "data_hist_untagged2_bkg", x, reader.u2_bkg_run2)
data_hist_untagged3_bkg = ROOT.RooDataHist("data_hist_untagged3_bkg", "data_hist_untagged3_bkg", x, reader.u3_bkg_run2)
data_hist_untagged4_bkg = ROOT.RooDataHist("data_hist_untagged4_bkg", "data_hist_untagged4_bkg", x, reader.u4_bkg_run2)

data_hist_untagged1_sig = ROOT.RooDataHist("data_hist_untagged1_sig", "data_hist_untagged1_sig", x, reader.u1_sig_run2)
data_hist_untagged2_sig = ROOT.RooDataHist("data_hist_untagged2_sig", "data_hist_untagged2_sig", x, reader.u2_sig_run2)
data_hist_untagged3_sig = ROOT.RooDataHist("data_hist_untagged3_sig", "data_hist_untagged3_sig", x, reader.u3_sig_run2)
data_hist_untagged4_sig = ROOT.RooDataHist("data_hist_untagged4_sig", "data_hist_untagged4_sig", x, reader.u4_sig_run2)



singleBernFTest(x, mu_gauss, data_hist_untagged1_bkg, "u1", args.method, "Poisson", False)
# singleBernFTest(x, mu_gauss, data_hist_untagged2_bkg, "u2", args.method, "Poisson", False)
# singleBernFTest(x, mu_gauss, data_hist_untagged3_bkg, "u3", args.method, "Poisson", False)
# singleBernFTest(x, mu_gauss, data_hist_untagged4_bkg, "u4", args.method, "Poisson", False)




