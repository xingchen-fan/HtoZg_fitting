#!/usr/bin/env python3
import ROOT
import os
import json
import sys
import argparse
sys.path.append(os.path.abspath("../Utilities/"))
from bkg_functions_class import *
from sig_functions_class import *
from Xc_Minimizer import *
from plot_utility import *
from sample_reader import *
#ROOT.gInterpreter.AddIncludePath('../Utilities/HZGRooPdfs.h')
# ROOT.gInterpreter.AddIncludePath('../Utilities/ModGaus.h')
#ROOT.gSystem.Load('../Utilities/HZGRooPdfs_cxx.so')
# ROOT.gSystem.Load('../Utilities/ModGaus_cxx.so')
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)

##############################################
# This script prepares the workspaces as input of combine. 
# Combine is needed and custom classes should be installed in combine instead!!!
# Each function in the profile needs to be fitted before entering the workspace.
# Signal model is fitted and fixed.
##############################################

parser = argparse.ArgumentParser(description = "Make workspace")
parser.add_argument('-c', '--cat', help="category")
args = parser.parse_args()
jfile = open('../Config/config.json', 'r')
configs = json.load(jfile)
CAT = args.cat
setting = configs[CAT]
lowx = setting["Range"]

# Save all the fit results?
LOG = True

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
MH = ROOT.RooRealVar("MH","MH"       ,124.7, 120., 130.)
list = [x, y, w, w_year, bdt, year, lep, ph_eta, nlep, njet]

# Cornell MC and data sample dat reader and make RooDataHist
# Zebing hist data reader
x.setBins(260)
x.setRange('left', lowx, 120)
x.setRange('right', 130, lowx+65)
x.setRange('full', lowx, lowx+65)
#ZBreader = readWsp(x, '/afs/cern.ch/user/f/fanx/EOS_space/zebing_sample/HZGamma_data_bkg_workspace_cat2.root', 'data_mass_cat0')

file_open_sig = ROOT.TFile.Open('../Data/sst_ggf_sig_hist_drop.root', 'READ')

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

hist_sig_TH1 = file_open_sig.Get('hist_sig_'+CAT)
hist_sig = ROOT.RooDataHist('hist_sig_'+CAT, 'hist_sig_'+CAT, x, hist_sig_TH1)
N_sig = hist_sig.sumEntries()
N = data_hist.sumEntries()
print("N sig = ", N_sig)

# Assume we have ...... in the profile, the signal model is DSCB
'''
core_bern2_model = CoreBern2Class(x, cat, 1, 0.3, 10, fileName="ZGCoreShape_01jet_NAFCorr", shapeName="CoreShape_ZG_NAF_cat2")
core_bern3_model = CoreBern3Class(x, cat, 1, 0.3, 10, fileName="ZGCoreShape_01jet_NAFCorr", shapeName="CoreShape_ZG_NAF_cat2")
core_bern4_model = CoreBern4Class(x, cat, 1, 0.3, 10, fileName="ZGCoreShape_01jet_NAFCorr", shapeName="CoreShape_ZG_NAF_cat2")
core_bern5_model = CoreBern5Class(x, cat, 1, 0.3, 10, fileName="ZGCoreShape_01jet_NAFCorr", shapeName="CoreShape_ZG_NAF_cat2")

core_pow1_model = CorePow1Class(x, cat, -1, fileName="ZGCoreShape_01jet_NAFCorr", shapeName="CoreShape_ZG_NAF_cat2")
core_pow2_model = CorePow2Class(x, cat, -1, -1, 1,  fileName="ZGCoreShape_01jet_NAFCorr", shapeName="CoreShape_ZG_NAF_cat2")
'''
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

sig_model = DSCB_Class(x, MH, CAT)

x.setRange("signal",115, 132)
sig_model.pdf.fitTo(hist_sig, ROOT.RooFit.SumW2Error(True),  ROOT.RooFit.Save(True),ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0), ROOT.RooFit.Range("signal"))
sig_model.setConst(True)
MH.setConstant(True)

profile = []
if "bern2" in setting["FT"]: profile.append(bern2_model)
if "bern3" in setting["FT"]: profile.append(bern3_model)
if "bern4" in setting["FT"]: profile.append(bern4_model)
if "bern5" in setting["FT"]: profile.append(bern5_model)
if "pow1" in setting["FT"]: profile.append(pow1_model)
if "pow2" in setting["FT"]: profile.append(pow2_model)
if "pow3" in setting["FT"]: profile.append(pow3_model)
if "exp1" in setting["FT"]: profile.append(exp1_model)
if "exp2" in setting["FT"]: profile.append(exp2_model)
if "exp3" in setting["FT"]: profile.append(exp3_model)
if "lau2" in setting["FT"]: profile.append(lau2_model)
if "lau3" in setting["FT"]: profile.append(lau3_model)
if "lau4" in setting["FT"]: profile.append(lau4_model)
if "modg" in setting["FT"]: profile.append(modg_model)

stat_list = []
cuthist = data_hist.reduce(ROOT.RooFit.CutRange('left,right'))
error = ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson)
for model in profile:
    stat_list.append(ROOT.RooNLLVar("stat_"+model.pdf.GetName(), "stat_"+model.pdf.GetName(), model.SBpdf,  cuthist))

eps = 0.1
strategy = 0
for stat in stat_list:
    Minimizer_Chi2(stat, -1, 100, False, strategy)
    #Minimizer_Chi2(stat, -1, 1, False, strategy)
    r = Minimizer_Chi2(stat, -1, eps, True, strategy)
    if LOG: r.Print("V")
    else: print(stat.GetName(), ": Cov q = ", r.covQual(), " status = ", r.status())

for model in profile:
    model.checkBond()

# Create signal workspace
f_out1 = ROOT.TFile("workspaces/workspace_sig_" + CAT + ".root", "RECREATE")
w_sig = ROOT.RooWorkspace("workspace_sig","workspace_sig")
getattr(w_sig, "import")(sig_model.pdf)
w_sig.Print()
w_sig.Write()
f_out1.Close()

# Create background workspace 
cate = ROOT.RooCategory("pdfindex_"+CAT, "Index of Pdf which is active for "+CAT)
models = ROOT.RooArgList()
for model in profile:
    models.add(model.pdf)
multipdf = ROOT.RooMultiPdf("multipdf_"+CAT, "MultiPdf for "+CAT, cate, models)

# Penalty term
#multipdf.setCorrectionFactor(0.)

norm = ROOT.RooRealVar("multipdf_"+ CAT +"_norm", "Number of background events", N, 0, 3*N)
f_out2 = ROOT.TFile("workspaces/workspace_bkg_profile_bias_" + CAT + ".root", "RECREATE")
w_bkg = ROOT.RooWorkspace("workspace_bkg","workspace_bkg")
getattr(w_bkg, "import")(cate)
getattr(w_bkg, "import")(norm)
getattr(w_bkg, "import")(multipdf)
getattr(w_bkg, "import")(data_hist)
w_bkg.Print()
w_bkg.Write()
f_out2.Close()

# Plot it
multiPlotClass(x, data_hist, profile, title="Profile_" +CAT, output_dir="", sideBand=True, fitRange= 'left,right')

# Create toy histogram (depreciated)
bias = False
N_toy = 10

if bias:
    f_out3 = ROOT.TFile("~/EOS_space/toy_wsp/workspace_toy_" + CAT + ".root", "RECREATE")
    w_toy = ROOT.RooWorkspace("workspace_toy","workspace_toy")
    for entry in profile:
        x.setBins(260)
        for i in range(N_toy):
            hist_toy = entry.pdf.generateBinned(x, ROOT.RooFit.NumEvents(N))
            hist_toy.SetNameTitle("hist_"+entry.pdf.GetName()+"_"+str(i), "hist_"+entry.pdf.GetName()+"_"+str(i))
            getattr(w_toy, "import")(hist_toy)
    w_toy.Write()
    w_toy.Print()
    f_out3.Close()


