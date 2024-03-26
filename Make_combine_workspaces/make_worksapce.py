import ROOT
import os
import sys
import argparse
sys.path.append(os.path.abspath("../Utilities/"))
from bkg_functions_class import *
from sig_functions_class import *
from Xc_Minimizer import *
from plot_utility import *
from sample_reader import *
# ROOT.gInterpreter.AddIncludePath('../Utilities/HZGRooPdfs.h')
# ROOT.gInterpreter.AddIncludePath('../Utilities/ModGaus.h')
# ROOT.gSystem.Load('../Utilities/HZGRooPdfs_cxx.so')
# ROOT.gSystem.Load('../Utilities/ModGaus_cxx.so')
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)

##############################################
# This script prepares the workspaces as input of combine. 
# Combine is needed and custom classes should be installed in combine instead!!!
# Each function in the profile needs to be fitted before entering the workspace.
# Signal model is fitted and fixed.
##############################################

# Save all the fit results?
LOG = False

# Specify the lower bond of mllg
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
MH = ROOT.RooRealVar("MH","MH"       ,124.7, 120., 130.)
list = [x, y, w, bdt, year, lep, ph_eta, nlep, njet]

# Cornell MC sample dat reader and make RooDataHist
x.setBins(260)
reader = readDat(list, "../../sample/")
reader.numCheck()
cat = "u1"
sig_hist = reader.data_hist_untagged1_sig
bkg_hist = reader.data_hist_untagged1_bkg
tot_hist = reader.data_hist_untagged1

# Assume we have Bern2, Bern3, Pow1, Exp1, Exp2, Lau1, Lau2 in the profile, the signal model is DSCB
bern2 = Bern2Class(x, mu_gauss, cat)
bern3 = Bern2Class(x, mu_gauss, cat)
pow1 = Pow1Class(x, mu_gauss, cat)
exp1 = Exp1Class(x, mu_gauss, cat)
exp2 = Exp2Class(x, mu_gauss, cat)
lau1 = Lau1Class(x, mu_gauss, cat)
lau2 = Lau2Class(x, mu_gauss, cat)
modg = ModGausClass(x,cat, lowx, lowx+65)

sig_model = DSCB_Class(x, MH, cat)

x.setBins(260)
#reader = readDat(list, "/afs/cern.ch/user/f/fanx/public/samples/")
reader = readDat(list, "../../sample/")
reader.numCheck()

x.setRange("signal",115, 132)
sig_model.pdf.fitTo(sig_hist, ROOT.RooFit.SumW2Error(True),  ROOT.RooFit.Save(True),ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0), ROOT.RooFit.Range("signal"))
sig_model.setConst(True)

profile = [bern2, bern3, pow1, exp1, exp2, lau1, lau2, modg]
stat_list = []
error = ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson)
for model in profile:
    stat_list.append(ROOT.RooChi2Var("stat_"+model.pdf.GetName(), "stat_"+model.pdf.GetName(), model.pdf,  bkg_hist, error))

eps = 0.1
strategy = 0
for stat in stat_list:
    Minimizer_Chi2(stat, -1, 100, False, strategy)
    r = Minimizer_Chi2(stat, -1, eps, False, strategy)
    if LOG: r.Print("V")
    else: print(stat.GetName(), ": Cov q = ", r.covQual(), " status = ", r.status())

for model in profile:
    model.checkBond()

# Create signal workspace
f_out1 = ROOT.TFile("workspaces/workspace_sig_" + cat + ".root", "RECREATE")
w_sig = ROOT.RooWorkspace("workspace_sig","workspace_sig")
getattr(w_sig, "import")(sig_model)
w_sig.Print()
w_sig.Write()
f_out1.Close()

# Create background workspace 
cate = ROOT.RooCategory("pdfindex_"+cat, "Index of Pdf which is active for "+cat)
models = ROOT.RooArgList()
for model in profile:
    models.add(model.pdf)
multipdf = ROOT.RooMultiPdf("multipdf_"+cat, "MultiPdf for "+cat, cate, models)

# Penalty term
#multipdf.setCorrectionFactor(0.)

norm = ROOT.RooRealVar("multipdf_"+ cat +"_norm", "Number of background events", tot_hist.sumEntries(), 0, 3*tot_hist.sumEntries())
f_out2 = ROOT.TFile("workspaces/workspace_bkg_profile" + cat + ".root", "RECREATE")
w_bkg = ROOT.RooWorkspace("workspace_bkg","workspace_bkg")
getattr(w_bkg, "import")(tot_hist)
getattr(w_bkg, "import")(cate)
getattr(w_bkg, "import")(norm)
getattr(w_bkg, "import")(multipdf)
w_bkg.Print()
w_bkg.Write()
f_out2.Close()




