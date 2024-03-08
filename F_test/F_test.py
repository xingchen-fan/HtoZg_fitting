import ROOT
import os
import sys
import argparse
sys.path.append(os.path.abspath("../Utilities/"))
from bkg_functions import *
from Xc_Minimizer import *

parser = argparse.ArgumentParser(description = "F test method (NLL or Chi2)")
parser.add_argument("method")
args = parser.parse_args()
if not(args.method == "Chi2" or args.method == "NLL") :
    print("Please use the correct method.")
    sys.exit(1)

def errorType(e_type):
    if e_type == "Poisson":
        return ROOT.RooAbsData.Poisson
    if e_type == "SumW2":
        return ROOT.RooAbsData.SumW2

def singleBernFTest_Chi2(x, gauss_mu, histogram, cat = "", e_type = "Poisson"):
    bern2_model = Bern2(x, gauss_mu, cat, 10, 0.3, 50, 7., 105.)
    bern2chi2 = ROOT.RooChi2Var("bern2chi2_" + cat, "bern2chi2_" + cat, bern2_model, histogram, ROOT.RooFit.DataError(errorType(e_type)))
    bern3_model = Bern3(x, gauss_mu, cat, 10, 0.3, 50, 7., 105.)
    bern3chi2 = ROOT.RooChi2Var("bern2chi2_" + cat, "bern2chi2_" + cat, bern3_model, histogram, ROOT.RooFit.DataError(errorType(e_type)))
    bern4_model = Bern4(x, gauss_mu, cat, 10, 0.3, 50, 7., 105.)
    bern4chi2 = ROOT.RooChi2Var("bern2chi2_" + cat, "bern2chi2_" + cat, bern4_model, histogram, ROOT.RooFit.DataError(errorType(e_type)))
    bern5_model = Bern5(x, gauss_mu, cat, 10, 0.3, 50, 7., 105.)
    bern5chi2 = ROOT.RooChi2Var("bern2chi2_" + cat, "bern2chi2_" + cat, bern5_model, histogram, ROOT.RooFit.DataError(errorType(e_type)))

    Minimizer_Chi2(bern2chi2, -1, 100, False, 0)
    r1 = Minimizer_Chi2(bern2chi2, -1, 0.1, False, 0)
    r1.Print("v")
    Minimizer_Chi2(bern3chi2, -1, 100, False, 0)
    r2 = Minimizer_Chi2(bern3chi2, -1, 0.1, False, 0)
    r2.Print("v")
    Minimizer_Chi2(bern4chi2, -1, 100, False, 0)
    r3 = Minimizer_Chi2(bern4chi2, -1, 0.1, False, 0)
    r3.Print("v")
    Minimizer_Chi2(bern5chi2, -1, 100, False, 0)
    r4 = Minimizer_Chi2(bern5chi2, -1, 0.1, False, 0)
    r4.Print("v")
    chi2s = [bern2chi2.geatVal(), bern3chi2.geatVal(), bern4chi2.geatVal(), bern5chi2.geatVal()]
    fs = []
    for i in range(len(chi2s) - 1):
        fs.append(ROOT.Math.fdistribution_cdf_c((chi2s[i] - chi2s[i+1])*(260-5-i)/chi2s[i+1], 1, 260-5-i))
    
    print("Chi2 = ", chi2s)
    print("P-value = ", fs)

def singleBernFTest_NLL(x, gauss_mu, histogram, cat = ""):
    bern2_model = Bern2(x, gauss_mu, cat, 10, 0.3, 50, 7., 105.)
    bern2NLL = ROOT.RooNLLVar("bern2NLL_" + cat, "bern2NLL_" + cat, bern2_model, histogram)
    bern3_model = Bern3(x, gauss_mu, cat, 10, 0.3, 50, 7., 105.)
    bern3NLL = ROOT.RooNLLVar("bern2NLL_" + cat, "bern2NLL_" + cat, bern3_model, histogram)
    bern4_model = Bern4(x, gauss_mu, cat, 10, 0.3, 50, 7., 105.)
    bern4NLL = ROOT.RooNLLVar("bern2NLL_" + cat, "bern2NLL_" + cat, bern4_model, histogram)
    bern5_model = Bern5(x, gauss_mu, cat, 10, 0.3, 50, 7., 105.)
    bern5NLL = ROOT.RooNLLVar("bern2NLL_" + cat, "bern2NLL_" + cat, bern5_model, histogram)

    Minimizer_NLL(bern2NLL, -1, 100, False, 0)
    r1 = Minimizer_NLL(bern2NLL, -1, 0.1, True, 0)
    r1.Print("v")
    Minimizer_NLL(bern3NLL, -1, 100, False, 0)
    r2 = Minimizer_NLL(bern3NLL, -1, 0.1, True, 0)
    r2.Print("v")
    Minimizer_NLL(bern4NLL, -1, 100, False, 0)
    r3 = Minimizer_NLL(bern4NLL, -1, 0.1, True, 0)
    r3.Print("v")
    Minimizer_NLL(bern5NLL, -1, 100, False, 0)
    r4 = Minimizer_NLL(bern5NLL, -1, 0.1, True, 0)
    r4.Print("v")
    NLLs = [bern2NLL.geatVal(), bern3NLL.geatVal(), bern4NLL.geatVal(), bern5NLL.geatVal()]
    fs = []
    for i in range(len(NLLs) - 1):
        fs.append(ROOT.Math.chisquared_cdf_c(2*(NLLs[i] - NLLs[i+1]), 1))
    
    print("NLL = ", NLLs)
    print("P-value = ", fs)






ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
# ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
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

# Read dat samples
bkg_run2 = ROOT.RooDataSet.read("../SMZg_deathvalley_v3_untagged.dat, ../DY_deathvalley_v3_untagged.dat", ROOT.RooArgList(x, y, bdt, w, year, lep, ph_eta, nlep, njet))
sig_run2 = ROOT.RooDataSet.read("../FullSig_deathvalley_v3_untagged.dat", ROOT.RooArgList(x, y, bdt, w, year, lep, ph_eta, nlep, njet))
tot_run2 = ROOT.RooDataSet.read("../SMZg_deathvalley_v3_untagged.dat, ../DY_deathvalley_v3_untagged.dat, ../FullSig_deathvalley_v3_untagged.dat", ROOT.RooArgList(x, y, bdt, w, year, lep, ph_eta, nlep, njet))

bdt1 = -0.36
bdt2 = -0.06
bdt3 = 0.08
bdt4 = 0.15

# Make RooDataSet
u1_bkg_run2_el = ROOT.RooDataSet("u1_bkg_run2_el", "u1_bkg_run2_el", bkg_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep < 1 && nlep <= 2 && njet < 2 && bdt > " + str(bdt4), "w")
u2_bkg_run2_el = ROOT.RooDataSet("u2_bkg_run2_el", "u2_bkg_run2_el", bkg_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep < 1 && nlep <= 2 && njet < 2 && bdt > " + str(bdt3) +" && bdt < " + str(bdt4), "w")
u3_bkg_run2_el = ROOT.RooDataSet("u3_bkg_run2_el", "u3_bkg_run2_el", bkg_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep < 1 && nlep <= 2 && njet < 2 && bdt > " + str(bdt2) +" && bdt < " + str(bdt3), "w")
u4_bkg_run2_el = ROOT.RooDataSet("u4_bkg_run2_el", "u4_bkg_run2_el", bkg_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep < 1 && nlep <= 2 && njet < 2 && bdt > " + str(bdt1) +" && bdt < " + str(bdt2), "w")

u1_sig_run2_el = ROOT.RooDataSet("u1_sig_run2_el", "u1_sig_run2_el", sig_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep < 1 && nlep <= 2 && njet < 2 && bdt > " + str(bdt4), "w")
u2_sig_run2_el = ROOT.RooDataSet("u2_sig_run2_el", "u2_sig_run2_el", sig_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep < 1 && nlep <= 2 && njet < 2 && bdt > " + str(bdt3) +" && bdt < " + str(bdt4), "w")
u3_sig_run2_el = ROOT.RooDataSet("u3_sig_run2_el", "u3_sig_run2_el", sig_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep < 1 && nlep <= 2 && njet < 2 && bdt > " + str(bdt2) +" && bdt < " + str(bdt3), "w")
u4_sig_run2_el = ROOT.RooDataSet("u4_sig_run2_el", "u4_sig_run2_el", sig_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep < 1 && nlep <= 2 && njet < 2 && bdt > " + str(bdt1) +" && bdt < " + str(bdt2), "w")

u1_tot_run2_el = ROOT.RooDataSet("u1_tot_run2_el", "u1_tot_run2_el", tot_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep < 1 && nlep <= 2 && njet < 2 && bdt > " + str(bdt4), "w")
u2_tot_run2_el = ROOT.RooDataSet("u2_tot_run2_el", "u2_tot_run2_el", tot_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep < 1 && nlep <= 2 && njet < 2 && bdt > " + str(bdt3) +" && bdt < " + str(bdt4), "w")
u3_tot_run2_el = ROOT.RooDataSet("u3_tot_run2_el", "u3_tot_run2_el", tot_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep < 1 && nlep <= 2 && njet < 2 && bdt > " + str(bdt2) +" && bdt < " + str(bdt3), "w")
u4_tot_run2_el = ROOT.RooDataSet("u4_tot_run2_el", "u4_tot_run2_el", tot_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep < 1 && nlep <= 2 && njet < 2 && bdt > " + str(bdt1) +" && bdt < " + str(bdt2), "w")

u1_bkg_run2_mu = ROOT.RooDataSet("u1_bkg_run2_mu", "u1_bkg_run2_mu", bkg_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep > 0 && nlep <= 2 && njet < 2 && bdt > " + str(bdt4), "w")
u2_bkg_run2_mu = ROOT.RooDataSet("u2_bkg_run2_mu", "u2_bkg_run2_mu", bkg_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep > 0 && nlep <= 2 && njet < 2 && bdt > " + str(bdt3) +" && bdt < " + str(bdt4), "w")
u3_bkg_run2_mu = ROOT.RooDataSet("u3_bkg_run2_mu", "u3_bkg_run2_mu", bkg_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep > 0 && nlep <= 2 && njet < 2 && bdt > " + str(bdt2) +" && bdt < " + str(bdt3), "w")
u4_bkg_run2_mu = ROOT.RooDataSet("u4_bkg_run2_mu", "u4_bkg_run2_mu", bkg_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep > 0 && nlep <= 2 && njet < 2 && bdt > " + str(bdt1) +" && bdt < " + str(bdt2), "w")

u1_sig_run2_mu = ROOT.RooDataSet("u1_sig_run2_mu", "u1_sig_run2_mu", sig_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep > 0 && nlep <= 2 && njet < 2 && bdt > " + str(bdt4), "w")
u2_sig_run2_mu = ROOT.RooDataSet("u2_sig_run2_mu", "u2_sig_run2_mu", sig_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep > 0 && nlep <= 2 && njet < 2 && bdt > " + str(bdt3) +" && bdt < " + str(bdt4), "w")
u3_sig_run2_mu = ROOT.RooDataSet("u3_sig_run2_mu", "u3_sig_run2_mu", sig_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep > 0 && nlep <= 2 && njet < 2 && bdt > " + str(bdt2) +" && bdt < " + str(bdt3), "w")
u4_sig_run2_mu = ROOT.RooDataSet("u4_sig_run2_mu", "u4_sig_run2_mu", sig_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep > 0 && nlep <= 2 && njet < 2 && bdt > " + str(bdt1) +" && bdt < " + str(bdt2), "w")

u1_tot_run2_mu = ROOT.RooDataSet("u1_tot_run2_mu", "u1_tot_run2_mu", tot_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep > 0 && nlep <= 2 && njet < 2 && bdt > " + str(bdt4), "w")
u2_tot_run2_mu = ROOT.RooDataSet("u2_tot_run2_mu", "u2_tot_run2_mu", tot_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep > 0 && nlep <= 2 && njet < 2 && bdt > " + str(bdt3) +" && bdt < " + str(bdt4), "w")
u3_tot_run2_mu = ROOT.RooDataSet("u3_tot_run2_mu", "u3_tot_run2_mu", tot_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep > 0 && nlep <= 2 && njet < 2 && bdt > " + str(bdt2) +" && bdt < " + str(bdt3), "w")
u4_tot_run2_mu = ROOT.RooDataSet("u4_tot_run2_mu", "u4_tot_run2_mu", tot_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep > 0 && nlep <= 2 && njet < 2 && bdt > " + str(bdt1) +" && bdt < " + str(bdt2), "w")

u1_bkg_run2 = ROOT.RooDataSet("u1_bkg_run2", "u1_bkg_run2", bkg_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "nlep <= 2 && njet < 2 && bdt > " + str(bdt4), "w")
u2_bkg_run2 = ROOT.RooDataSet("u2_bkg_run2", "u2_bkg_run2", bkg_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "nlep <= 2 && njet < 2 && bdt > " + str(bdt3) +" && bdt < " + str(bdt4), "w")
u3_bkg_run2 = ROOT.RooDataSet("u3_bkg_run2", "u3_bkg_run2", bkg_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "nlep <= 2 && njet < 2 && bdt > " + str(bdt2) +" && bdt < " + str(bdt3), "w")
u4_bkg_run2 = ROOT.RooDataSet("u4_bkg_run2", "u4_bkg_run2", bkg_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "nlep <= 2 && njet < 2 && bdt > " + str(bdt1) +" && bdt < " + str(bdt2), "w")

u1_sig_run2 = ROOT.RooDataSet("u1_sig_run2", "u1_sig_run2", sig_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "nlep <= 2 && njet < 2 && bdt >" + str(bdt4), "w")
u2_sig_run2 = ROOT.RooDataSet("u2_sig_run2", "u2_sig_run2", sig_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "nlep <= 2 && njet < 2 && bdt > " + str(bdt3) +" && bdt < " + str(bdt4), "w")
u3_sig_run2 = ROOT.RooDataSet("u3_sig_run2", "u3_sig_run2", sig_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "nlep <= 2 && njet < 2 && bdt > " + str(bdt2) +" && bdt < " + str(bdt3), "w")
u4_sig_run2 = ROOT.RooDataSet("u4_sig_run2", "u4_sig_run2", sig_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "nlep <= 2 && njet < 2 && bdt > " + str(bdt1) +" && bdt < " + str(bdt2), "w")

u1_tot_run2 = ROOT.RooDataSet("u1_tot_run2", "u1_tot_run2", tot_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "nlep <= 2 && njet < 2 && bdt > " + str(bdt4), "w")
u2_tot_run2 = ROOT.RooDataSet("u2_tot_run2", "u2_tot_run2", tot_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "nlep <= 2 && njet < 2 && bdt > " + str(bdt3) +" && bdt < " + str(bdt4) , "w")
u3_tot_run2 = ROOT.RooDataSet("u3_tot_run2", "u3_tot_run2", tot_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "nlep <= 2 && njet < 2 && bdt > " + str(bdt2) +" && bdt < " + str(bdt3) , "w")
u4_tot_run2 = ROOT.RooDataSet("u4_tot_run2", "u4_tot_run2", tot_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "nlep <= 2 && njet < 2 && bdt > " + str(bdt1) +" && bdt < " + str(bdt2) , "w")

# Expected bkg
expbkg_u1 = u1_bkg_run2.sumEntries()
expbkg_u2 = u2_bkg_run2.sumEntries()
expbkg_u3 = u3_bkg_run2.sumEntries()
expbkg_u4 =  u4_bkg_run2.sumEntries()

# Make RooDataHist
x.setBins(260)
data_hist_untagged1 = ROOT.RooDataHist("data_hist_untagged1", "data_hist_untagged1", x, u1_tot_run2)
data_hist_untagged2 = ROOT.RooDataHist("data_hist_untagged2", "data_hist_untagged2", x, u2_tot_run2)
data_hist_untagged3 = ROOT.RooDataHist("data_hist_untagged3", "data_hist_untagged3", x, u3_tot_run2)
data_hist_untagged4 = ROOT.RooDataHist("data_hist_untagged4", "data_hist_untagged4", x, u4_tot_run2)

data_hist_untagged1_bkg = ROOT.RooDataHist("data_hist_untagged1_bkg", "data_hist_untagged1_bkg", x, u1_bkg_run2)
data_hist_untagged2_bkg = ROOT.RooDataHist("data_hist_untagged2_bkg", "data_hist_untagged2_bkg", x, u2_bkg_run2)
data_hist_untagged3_bkg = ROOT.RooDataHist("data_hist_untagged3_bkg", "data_hist_untagged3_bkg", x, u3_bkg_run2)
data_hist_untagged4_bkg = ROOT.RooDataHist("data_hist_untagged4_bkg", "data_hist_untagged4_bkg", x, u4_bkg_run2)

data_hist_untagged1_sig = ROOT.RooDataHist("data_hist_untagged1_sig", "data_hist_untagged1_sig", x, u1_sig_run2)
data_hist_untagged2_sig = ROOT.RooDataHist("data_hist_untagged2_sig", "data_hist_untagged2_sig", x, u2_sig_run2)
data_hist_untagged3_sig = ROOT.RooDataHist("data_hist_untagged3_sig", "data_hist_untagged3_sig", x, u3_sig_run2)
data_hist_untagged4_sig = ROOT.RooDataHist("data_hist_untagged4_sig", "data_hist_untagged4_sig", x, u4_sig_run2)


if args.method == "Chi2":
    singleBernFTest_Chi2(x, u_gauss, data_hist_untagged1_bkg, "u1", "Poisson")
    singleBernFTest_Chi2(x, u_gauss, data_hist_untagged2_bkg, "u2", "Poisson")
    singleBernFTest_Chi2(x, u_gauss, data_hist_untagged3_bkg, "u3", "Poisson")
    singleBernFTest_Chi2(x, u_gauss, data_hist_untagged4_bkg, "u4", "Poisson")


elif args.method == "NLL":
    singleBernFTest_NLL(x, u_gauss, data_hist_untagged1_bkg, "u1")
    singleBernFTest_NLL(x, u_gauss, data_hist_untagged2_bkg, "u2")
    singleBernFTest_NLL(x, u_gauss, data_hist_untagged3_bkg, "u3")
    singleBernFTest_NLL(x, u_gauss, data_hist_untagged4_bkg, "u4") 



