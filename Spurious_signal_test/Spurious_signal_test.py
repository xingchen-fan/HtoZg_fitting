import ROOT
import os
import sys
import argparse
sys.path.append(os.path.abspath("../Utilities/"))
from bkg_functions_fit import *
from bkg_functions_class import *
from sig_functions_class import *
from Xc_Minimizer import *
from plot_utility import *
from sample_reader import *

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
reader = readDat(list, "../../../../CMSSW_11_3_4/")

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

hist_sig = data_hist_untagged1_sig
hist_bkg = data_hist_untagged1_bkg
# Define bkg and sig model
cat = "u1"
MH = ROOT.RooRealVar("MH","MH"       ,124.7, 120., 130.)
bkg_model = Bern2Class(x, mu_gauss, cat, 10, 0.3, 10, 7., 105.)
dscb_model = DSCB_Class(x, MH, cat)
# Signal fit
x.setRange("signal",110, 132)
res = dscb_model.pdf.fitTo(hist_sig, ROOT.RooFit.SumW2Error(True),\
                         ROOT.RooFit.Save(True),ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0), ROOT.RooFit.Range("signal"))
res.Print("V")
plotClass(x, hist_sig, dscb_model.pdf, "Signal")
dscb_model.setConst(True)

# Bkg fit
chi2 = ROOT.RooChi2Var("chi2_bkg_" + cat, "chi2 bkg " + cat, bkg_model.pdf,  hist_bkg, ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
Minimizer_Chi2(chi2, -1, 100, False, 0)
r = Minimizer_Chi2(chi2, -1, 0.1, False, 0)
r.Print("V")

# S+B fit
c1 = ROOT.RooRealVar("c1", "c1", hist_bkg.sumEntries(), 0, 3.* hist_bkg.sumEntries())
c2 = ROOT.RooRealVar("c2", "c2", 0., -1000., 1000)

tot_model = ROOT.RooAddPdf("tot_model", "tot_model", ROOT.RooArgList(dscb_model.pdf, bkg_model.pdf), ROOT.RooArgList(c2, c1))
chi2_tot = ROOT.RooChi2Var("chi2_tot", "chi2_tot", tot_model, hist_bkg, ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson), ROOT.RooFit.Extended(True))
Minimizer_Chi2(chi2_tot, -1, 100, False, 0)
res1 = Minimizer_Chi2(chi2_tot, -1, 0.1, False, 0)
print("Poisson error:")
res1.Print("V")
NS = c2.getVal()
delta_NS = c2.getError()

chi2_tot = ROOT.RooChi2Var("chi2_tot", "chi2_tot", tot_model, hist_bkg, ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2), ROOT.RooFit.Extended(True))
Minimizer_Chi2(chi2_tot, -1, 100, False, 0)
res2 = Minimizer_Chi2(chi2_tot, -1, 0.1, False, 0)
print("SumW2 error:")
res2.Print("V")
delta_MC = c2.getError()

print("N signal = ", NS, ", Delta Ns = ", delta_NS, ", Delta MC = ", delta_MC)

plotClass(x, hist_bkg, tot_model, "S+B")





