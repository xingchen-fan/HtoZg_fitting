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
parser.add_argument('-a', '--asimov', help="Asimov", default=0)
parser.add_argument('-con', '--config', help = 'Configuration')
args = parser.parse_args()
jfile = open('../Config/'+args.config, 'r')
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
            hist_data = read_data.ggf1
        elif CAT == 'ggf2':
            hist_data = read_data.ggf2
        elif CAT == 'ggf3':
            hist_data = read_data.ggf3
        elif CAT == 'ggf4':
            hist_data = read_data.ggf4
    elif 'vbf' in CAT:
        read_data = readRuiROOTVBFdata(x, '/eos/user/r/rzou/SWAN_projects/Classifier/Output_2JClassic_input_run2p3/', 0.489,0.286 ,0.083)
        if CAT == 'vbf1':
            hist_data = read_data.vbf1
        elif CAT == 'vbf2':
            hist_data = read_data.vbf2
        elif CAT == 'vbf3':
            hist_data = read_data.vbf3
        elif CAT == 'vbf4':
            hist_data = read_data.vbf4


# Signal model preparation ------------------------
#hist_sig_TH1 = file_open_sig.Get('hist_sig_'+'ggf1')
#hist_sig = ROOT.RooDataHist('hist_sig_'+CAT, 'hist_sig_'+CAT, x, hist_sig_TH1)
sig_model = combineSignal(x, MH, CAT, '../Config/config_DSCB.json')
MH.setVal(125)
MH.setConstant(True)
N_sig = sig_model.ntot
N = hist_data.sumEntries()
#N_sig_window =  hist_sig.sumEntries('CMS_hzg_mass_' + CAT + ' > 120 && ' + 'CMS_hzg_mass_' + CAT+ ' < 130')
print("N sig = ", N_sig)
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
profile_ = profileClass(x, mu_gauss, CAT, '../Config/'+args.config)
profile = profile_.testSelection("FT")

stat_list = []
cuthist = hist_data.reduce(ROOT.RooFit.CutRange('left,right'))
error = ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson)
for model in profile:
    stat_list.append(ROOT.RooNLLVar("stat_"+model.pdf.GetName(), "stat_"+model.pdf.GetName(), model.SBpdf,  cuthist))

eps = 0.1
strategy = 0
stat_vals = []
nFloat = 0
for ind, stat in enumerate(stat_list):
    if FIT:
        Minimizer_NLL(stat, -1, 100, False, strategy)
        #Minimizer_Chi2(stat, -1, 1, False, strategy)
        r = Minimizer_NLL(stat, -1, eps, True, strategy)
        if LOG: r.Print("V")
        else: print(stat.GetName(), ": Cov q = ", r.covQual(), " status = ", r.status())
        nFloat = r.floatParsFinal().getSize()
    else:
        nFloat = profile[ind].pdf.getParameters(hist_data).selectByAttrib('Constant',False).getSize()
    stat_vals.append(stat.getVal() + 0.5*nFloat)
    print (profile[ind].name,' = ', nFloat)
        
for model in profile:
    model.checkBond()

best_ = stat_vals.index(min(stat_vals))
print ("best is ", profile[best_].name)

# Create Asimov (optional)
c1= ROOT.RooRealVar('c1','c1', N_sig)#, 0, 3*N_sig)
c2 = ROOT.RooRealVar('c2','c2', N)#, 0, 3*N)
tot_model = ROOT.RooAddPdf('tot_pdf', 'tot_pdf', ROOT.RooArgList(sig_model.pdf, profile[best_].pdf), ROOT.RooArgList(c1, c2))
hist_asimov = tot_model.generateBinned(ROOT.RooArgSet(x), N, ROOT.RooFit.Asimov(True))
hist_asimov.SetNameTitle('hist_' + CAT + '_asimov', 'hist_' + CAT + '_asimov')
'''
nll1 = ROOT.RooNLLVar('nll1', 'nll2',  profile[best_].pdf, hist_asimov)
r1 = Minimizer_NLL(nll1, -1, eps, True, strategy)
nll1_val = nll1.getVal()
nll2 = ROOT.RooNLLVar('nll1', 'nll2',  tot_model, hist_asimov)
r2 = Minimizer_NLL(nll2, -1, eps, True, strategy)
nll2_val = nll2.getVal()
'''

# Create signal workspace
if SIGNAL:
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

norm = ROOT.RooRealVar("multipdf_"+ CAT +"_norm", "Number of background events", N - N_sig, 0, 3*N)
if int(args.asimov) == 0:
    f_out2 = ROOT.TFile("workspaces/workspace_bkg_profile_bias_" + CAT + ".root", "RECREATE")
else:
    f_out2 = ROOT.TFile("workspaces/workspace_bkg_profile_bias_asimov_" + CAT + ".root", "RECREATE")
w_bkg = ROOT.RooWorkspace("workspace_bkg","workspace_bkg")
getattr(w_bkg, "import")(cate)
getattr(w_bkg, "import")(norm)
getattr(w_bkg, "import")(multipdf)
if int(args.asimov) == 0:
    getattr(w_bkg, "import")(hist_data)
else:
    getattr(w_bkg, "import")(hist_asimov)
w_bkg.Print()
w_bkg.Write()
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
multiPlotClass(x, hist_data, profile, title="Profile_" +CAT, output_dir="", sideBand=True, fitRange= 'left,right', best_index = best_, bestLabel = True)

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


