#!/usr/bin/env python3
import math
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
from profile_class import *
from sig_functions_class import *

#ROOT.gInterpreter.AddIncludePath('../Utilities/HZGRooPdfs.h')
#ROOT.gSystem.Load('../Utilities/HZGRooPdfs_cxx.so')
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)

##############################################
# This script prepares the workspaces as input of combine. 
# Combine is needed and custom classes should be installed in combine instead!!!
# Each function in the profile needs to be fitted before entering the workspace.
# Signal model is fitted and fixed.
##############################################

parser = argparse.ArgumentParser(description = "Make workspace")
parser.add_argument('-c', '--cat', help="category")
parser.add_argument('-con', '--config', help = 'Configuration')
parser.add_argument('-b', '--best', help = 'Best', type=int)
args = parser.parse_args()
jfile = open(args.config, 'r')
configs = json.load(jfile)
CAT = args.cat
setting = configs[CAT]
lowx = setting["Range"]

# Save all the fit results? Make signal workspace? Read dat files? Do fit?
LOG = True
SIGNAL = False
DAT = False
FIT = False # 'False' to read the post-fit values from config file

# Define variables
x = ROOT.RooRealVar("CMS_hzg_mass_"+CAT, "CMS_hzg_mass_"+CAT, lowx, lowx + 65.)
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
MH = ROOT.RooRealVar("MH","MH"       ,125, 120., 130.)
list = [x, y, w, w_year, bdt, year, lep, ph_eta, nlep, njet]

# Cornell MC and data sample dat reader and make RooDataHist
# Zebing hist data reader
x.setBins(260)
x.setRange('left', lowx, 120)
x.setRange('right', 130, lowx+65)
x.setRange('full', lowx, lowx+65)
#ZBreader = readWsp(x, '/afs/cern.ch/user/f/fanx/EOS_space/zebing_sample/HZGamma_data_bkg_workspace_cat2.root', 'data_mass_cat0')

#file_open_sig = ROOT.TFile.Open('../Data/sst_ggf_sig_hist_drop.root', 'READ')

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
        read_data = readRuiROOTggFdata(x, '/eos/user/r/rzou//SWAN_projects/Classifier/Output_ggF_pinnacles_fix/relpt_peking_run2p3/TreeB/', 0.81, 0.64, 0.47)
        if CAT == 'ggf1':
            hist_data_el = read_data.ggf1_el
            hist_data_mu = read_data.ggf1_mu
            hist_data = read_data.ggf1
        elif CAT == 'ggf2':
            hist_data_el = read_data.ggf2_el
            hist_data_mu = read_data.ggf2_mu
            hist_data = read_data.ggf2
        elif CAT == 'ggf3':
            hist_data_el = read_data.ggf3_el
            hist_data_mu = read_data.ggf3_mu
            hist_data = read_data.ggf3
        elif CAT == 'ggf4':
            hist_data_el = read_data.ggf4_el
            hist_data_mu = read_data.ggf4_mu
            hist_data = read_data.ggf4
    elif 'vbf' in CAT:
        read_data = readRuiROOTVBFdata(x, '/eos/user/r/rzou/SWAN_projects/Classifier/Output_2JClassic_input_run2p3/', 0.489,0.286 ,0.083)
        if CAT == 'vbf1':
            hist_data_el = read_data.vbf1_el
            hist_data_mu = read_data.vbf1_mu
            hist_data = read_data.vbf1
        elif CAT == 'vbf2':
            hist_data_el = read_data.vbf2_el
            hist_data_mu = read_data.vbf2_mu
            hist_data = read_data.vbf2
        elif CAT == 'vbf3':
            hist_data_el = read_data.vbf3_el
            hist_data_mu = read_data.vbf3_mu
            hist_data = read_data.vbf3
        elif CAT == 'vbf4':
            hist_data_el = read_data.vbf4_el
            hist_data_mu = read_data.vbf4_mu
            hist_data = read_data.vbf4


# Signal model preparation ------------------------
#hist_sig_TH1 = file_open_sig.Get('hist_sig_'+'ggf1')
#hist_sig = ROOT.RooDataHist('hist_sig_'+CAT, 'hist_sig_'+CAT, x, hist_sig_TH1)
sig_model_el = combineSignalLep(x, MH, CAT, '../Config/config_DSCB.json', 'el')
sig_model_mu = combineSignalLep(x, MH, CAT, '../Config/config_DSCB.json', 'mu')
sig_model = combineSignal(x, MH, CAT, '../Config/config_DSCB.json')
#sig_model = DSCB_Class(x, MH, CAT, sigmaL_init = 1.8457, sigmaR_init = 1.3264, nL_init = 3.771, nL_bond = 100, nR_init = 9.962, nR_bond = 100, alphaL_init = 1.154, alphaR_init = 1.231, di_sigma = True)
MH.setVal(125)
MH.setConstant(True)
N_sig_el = sig_model_el.ntot
N_sig_mu = sig_model_mu.ntot
N_el = hist_data_el.sumEntries()
N_mu = hist_data_mu.sumEntries()
N = hist_data.sumEntries()


#N_sig_window =  hist_sig.sumEntries('CMS_hzg_mass_' + CAT + ' > 120 && ' + 'CMS_hzg_mass_' + CAT+ ' < 130')
print("N sig el = ", N_sig_el)
print("N sig mu = ", N_sig_mu)
print ("N bkg = ", N)
#sig_model.setConst(True)
# -------------------------------------------------


#print("N sig window= ", N_sig_window)

# Assume we have ...... in the profile, the signal model is combined DSCB
'''
core_bern2_model = CoreBern2Class(x, cat, 1, 0.3, 10, fileName="ZGCoreShape_01jet_NAFCorr", shapeName="CoreShape_ZG_NAF_cat2")
core_bern3_model = CoreBern3Class(x, cat, 1, 0.3, 10, fileName="ZGCoreShape_01jet_NAFCorr", shapeName="CoreShape_ZG_NAF_cat2")
core_bern4_model = CoreBern4Class(x, cat, 1, 0.3, 10, fileName="ZGCoreShape_01jet_NAFCorr", shapeName="CoreShape_ZG_NAF_cat2")
core_bern5_model = CoreBern5Class(x, cat, 1, 0.3, 10, fileName="ZGCoreShape_01jet_NAFCorr", shapeName="CoreShape_ZG_NAF_cat2")

core_pow1_model = CorePow1Class(x, cat, -1, fileName="ZGCoreShape_01jet_NAFCorr", shapeName="CoreShape_ZG_NAF_cat2")
core_pow2_model = CorePow2Class(x, cat, -1, -1, 1,  fileName="ZGCoreShape_01jet_NAFCorr", shapeName="CoreShape_ZG_NAF_cat2")
'''
'''
x.setRange("signal",115, 132)
res = sig_model.pdf.fitTo(hist_sig, ROOT.RooFit.SumW2Error(True),  ROOT.RooFit.Save(True),ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0), ROOT.RooFit.Range("signal"), ROOT.RooFit.Save(True))
res.Print('v')
sig_model.setConst(True)
MH.setVal(125.)
MH.setConstant(True)
'''


#plotClass(x, hist_sig, sig_model.pdf, sig_model.pdf, "Signal_"+CAT, CMS = 'Simulation', output_dir="")
profile_ = profileClass(x, mu_gauss, CAT, args.config)
profile = profile_.testSelection("FT")


stat_list = []
cuthist = hist_data.reduce(ROOT.RooFit.CutRange('left,right'))
error = ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson)
best_ = args.best
stat = ROOT.RooNLLVar("stat_"+profile[best_].pdf.GetName(), "stat_"+profile[best_].pdf.GetName(), profile[best_].SBpdf,  cuthist)

eps = 0.1
strategy = 0
stat_vals = []
nFloat = 0

Minimizer_NLL(stat, -1, 100, False, strategy)
r = Minimizer_NLL(stat, -1, eps, True, strategy)
r.Print("V")

profile[best_].checkBond()


# Create signal workspace
if SIGNAL:
    f_out1 = ROOT.TFile("workspaces/workspace_sig_toy_unblind_" + CAT + ".root", "RECREATE")
    w_sig = ROOT.RooWorkspace("workspace_sig","workspace_sig")
    getattr(w_sig, "import")(sig_model.pdf)
    getattr(w_sig, "import")(sig_model_el.pdf)
    getattr(w_sig, "import")(sig_model_mu.pdf)
    w_sig.Print()
    w_sig.Write()
    f_out1.Close()

# Create toy samples
c1_el= ROOT.RooRealVar('c1_el','c1_el', N_sig_el)
c2_el = ROOT.RooRealVar('c2_el','c2_el', N_el - N_sig_el)
c1_mu= ROOT.RooRealVar('c1_mu','c1_mu', N_sig_mu)
c2_mu = ROOT.RooRealVar('c2_mu','c2_mu', N_mu - N_sig_mu)

tot_model_el = ROOT.RooAddPdf('tot_pdf_el', 'tot_pdf_el', ROOT.RooArgList(sig_model_el.pdf, profile[best_].pdf), ROOT.RooArgList(c1_el, c2_el))
tot_model_mu = ROOT.RooAddPdf('tot_pdf_mu', 'tot_pdf_mu', ROOT.RooArgList(sig_model_mu.pdf, profile[best_].pdf), ROOT.RooArgList(c1_mu, c2_mu))

toy_hist_el = tot_model_el.generateBinned(ROOT.RooArgSet(x), N_el)
toy_hist_mu = tot_model_mu.generateBinned(ROOT.RooArgSet(x), N_mu)
toy_hist_el.SetNameTitle('hist_' + CAT+ '_el', 'hist_'+ CAT+'_el')
toy_hist_mu.SetNameTitle('hist_' + CAT + '_mu', 'hist_' + CAT + '_mu')

# Create background workspaces
#lep_c = ROOT.RooRealVar('lep_c', 'lep_c', hist_data_el.sumEntries()/(hist_data_el.sumEntries() + hist_data_mu.sumEntries()), 0., 1.)
#lep_model = ROOT.RooAddPdf('lep_combine_model_'+CAT, 'lep_combine_model_'+CAT, profile_el[best_].pdf, profile_mu[best_].pdf, lep_c)
cate_el = ROOT.RooCategory("pdfindex_"+CAT+"_el", "Index of Pdf which is active for "+CAT+" el")
cate_mu = ROOT.RooCategory("pdfindex_"+CAT+"_mu", "Index of Pdf which is active for "+CAT+" mu")
cate = ROOT.RooCategory("pdfindex_"+CAT, "Index of Pdf which is active for "+CAT)

models = ROOT.RooArgList()
for model in profile:
    models.add(model.pdf)

multipdf_el = ROOT.RooMultiPdf("multipdf_"+CAT+"_el", "MultiPdf for "+CAT+" el", cate_el, models)
multipdf_mu = ROOT.RooMultiPdf("multipdf_"+CAT+"_mu", "MultiPdf for "+CAT+" mu", cate_mu, models)
multipdf = ROOT.RooMultiPdf("multipdf_"+CAT, "MultiPdf for "+CAT, cate, models)

# Penalty term
#multipdf.setCorrectionFactor(0.)
print('hist = ', hist_data.sumEntries())

norm_el = ROOT.RooRealVar("multipdf_"+ CAT +"_el_norm", "Number of background events", N_el, 0, 3*N_el)
norm_mu = ROOT.RooRealVar("multipdf_"+ CAT +"_mu_norm", "Number of background events", N_mu, 0, 3*N_mu)
norm = ROOT.RooRealVar("multipdf_"+ CAT +"_norm", "Number of background events", N, 0, 3*N)

f_out2 = ROOT.TFile("workspaces/workspace_bkg_toy_unblind_" + CAT+  ".root", "RECREATE")
w_bkg = ROOT.RooWorkspace("workspace_bkg","workspace_bkg")
w_bkg_el = ROOT.RooWorkspace("workspace_bkg_el","workspace_bkg_el")
w_bkg_mu = ROOT.RooWorkspace("workspace_bkg_mu","workspace_bkg_mu")

getattr(w_bkg_el, "import")(cate_el)
getattr(w_bkg_el, "import")(norm_el)
getattr(w_bkg_el, "import")(multipdf_el)
getattr(w_bkg_el, "import")(toy_hist_el)
#w_bkg.Print()
w_bkg_el.Write()

getattr(w_bkg_mu, "import")(cate_mu)
getattr(w_bkg_mu, "import")(norm_mu)
getattr(w_bkg_mu, "import")(multipdf_mu)
getattr(w_bkg_mu, "import")(toy_hist_mu)
#w_bkg.Print()
w_bkg_mu.Write()

getattr(w_bkg, "import")(cate)
getattr(w_bkg, "import")(norm)
getattr(w_bkg, "import")(multipdf)
toy_hist_el.add(toy_hist_mu)
toy_hist_el.SetNameTitle('hist_' + CAT, 'hist_' + CAT)
getattr(w_bkg, "import")(toy_hist_el)
#w_bkg.Print()
w_bkg.Write()
print ("N toy = ", toy_hist_el.sumEntries())

f_out2.Close()

# Cut and count significance
CDF = False
if CDF:
    cdf = profile[best_].pdf.createCdf(ROOT.RooArgSet(x))
    x.setVal(120.)
    cdf_low = cdf.getVal()
    x.setVal(130.)
    cdf_high = cdf.getVal()
    print ('N bkg (sig window) = ', N*(cdf_high - cdf_low))
    print ('Significance = ', N_sig_window/math.sqrt(N*(cdf_high - cdf_low)))


#plotClass(x, hist_asimov, tot_model, tot_model, title="Asimov", output_dir="", sideBand = False)

# Plot it
#multiPlotClass(x, hist_data, profile, title="Profile_" +CAT, output_dir="", sideBand=True, fitRange= 'left,right', best_index = best_, bestLabel = True)

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


