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
ROOT.gInterpreter.AddIncludePath('../Utilities/RooGaussStepBernstein.h')
ROOT.gSystem.Load('../Utilities/RooGaussStepBernstein_cxx.so')
#ROOT.gInterpreter.AddIncludePath('../Utilities/AsymGenGaussian.h')
#ROOT.gSystem.Load('../Utilities/AsymGenGaussian_cxx.so')
#ROOT.gInterpreter.AddIncludePath('../Utilities/EXModGaus.h')
#ROOT.gSystem.Load('../Utilities/EXModGaus_cxx.so')
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)

#################################################################################
# This script prepares the workspaces as input of combine.                      #
# Combine is needed and custom classes should be installed in combine instead!!!#
# Each function in the profile needs to be fitted before entering the workspace.#
# Signal model is fitted and fixed.                                             #
#################################################################################

parser = argparse.ArgumentParser(description = "Make workspace")
parser.add_argument('-c', '--cat', help="category")
parser.add_argument('-a', '--asimov', help="Asimov", default=0, type=int)
parser.add_argument('-conB', '--config', help = 'Configuration')
parser.add_argument('-conS', '--configS', help = 'Signal Configuration')
parser.add_argument('-t', '--test', help = 'What test to build?')
parser.add_argument('-tag', '--tag', help = 'Tag for compatibility test (i.e. el, mu, r2, r3)', default='')
parser.add_argument('-y', '--type', help = 'Data type', default='')
parser.add_argument('-n', '--name', help = 'Name for output files')
args = parser.parse_args()
jfile = open(args.config, 'r')
configs = json.load(jfile)
CAT = args.cat
CAT_m = CAT
CAT_merge = args.cat[:4]
if CAT[:5] == 'vbf12':
    CAT_merge = 'vbf2'
elif CAT[:5] == 'ggf12':
    CAT_merge = 'ggf2'
else:
    CAT_merge = CAT
    
if CAT == 'untag':
    CAT_m = 'untagged'
elif CAT == 'vhptmiss':
    CAT_m = 'vhmet'

if args.tag == '':
    TAG = ''
    setting = configs[CAT]
else:
    TAG = '_'+args.tag
    setting = configs[CAT_merge]
    
lowx = setting["Range"][0]
highx = setting["Range"][1]
nbins = int(setting["Bins"])
if not CAT[:4] in ["ggf1","ggf2","ggf3","ggf4","vbf1","vbf2","vbf3","vbf4", "vh3l", "tthl", "tthh", "vhpt", "unta"]:
    raise ValueError("Unknown category")
if not args.test in ["FT", "CMSBias", "pruneUnstable"]:
    raise ValueError("Wrong test")

# Save all the fit results? Make signal workspace? Read dat files? Do fit?
LOG = True
FIT = True # 'False' to read the post-fit values from config file

# Defualt penalty term 0.5
PENALTY = 0.5

# Define variables
x = ROOT.RooRealVar("CMS_hzg_mass_"+CAT_merge, "CMS_hzg_mass_"+CAT_merge, lowx, highx)
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
MH = ROOT.RooRealVar("MH","MH"       ,125.38, 120., 130.)
list = [x, y, w, w_year, bdt, year, lep, ph_eta, nlep, njet]

# Cornell MC and data sample dat reader and make RooDataHist
# Zebing hist data reader
x.setBins(nbins)
x.setRange('left', lowx, 120)
x.setRange('right', 130, highx)
x.setRange('full', lowx, highx)

#ZBreader = readWsp(x, '/afs/cern.ch/user/f/fanx/EOS_space/zebing_sample/HZGamma_data_bkg_workspace_cat2.root', 'data_mass_cat0')

#file_open_sig = ROOT.TFile.Open('../Data/sst_ggf_sig_hist_drop.root', 'READ')

if args.type == 'dat':
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
        
elif args.type == 'pico':
    read_data = readPico(x, '~/EOS_space/michael_files/hzg_datacard_v1p4p0_rawdata.root', CAT_m, "data_obs")
    hist_data = getattr(read_data, f"data_obs_cat_{CAT_m}")
else:
    if 'ggf' in CAT:
        read_data = readRuiROOTggFdata(x, '/eos/project/h/htozg-dy-privatemc/rzou/bdt/BDT_output_redwood/Output_ggF_rui_redwood_v1_ext_val/', 0.94,0.83,0.57)
        hist_data = getattr(read_data, CAT)

    elif 'vbf' in CAT:
        read_data = readRuiROOTVBFdata(x, '/eos/project/h/htozg-dy-privatemc/rzou/bdt/BDT_output_redwood/Output_VBF_rui_redwood_v1_ext_val/', 0.91,0.81,0.48)
        hist_data = getattr(read_data, CAT)
#print("N sig window= ", N_sig_window)
print ("done reading the data")
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
profile_ = profileClass(x, mu_gauss, CAT_merge, args.config, TAG)
profile = profile_.testSelection(args.test)
print("done profile")
stat_list = []
cuthist = hist_data.reduce(ROOT.RooFit.CutRange('left,right'))

error = ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson)
for model in profile:
    #model.pdf.SetNameTitle(model.pdf.GetName()+TAG,model.pdf.GetName()+TAG)
    stat_list.append(ROOT.RooNLLVar("stat_"+model.pdf.GetName(), "stat_"+model.pdf.GetName(), model.SBpdf,  cuthist))
print("stats done")

eps = 0.1
strategy = 0
stat_vals = []
nFloat = 0
for ind, stat in enumerate(stat_list):
    if FIT:
        Minimizer_NLL(stat, -1, 100, True, strategy)
        #Minimizer_Chi2(stat, -1, 1, False, strategy)
        r = Minimizer_NLL(stat, -1, eps, True, strategy)
        if LOG: r.Print("V")
        else: print(stat.GetName(), ": Cov q = ", r.covQual(), " status = ", r.status())
        nFloat = r.floatParsFinal().getSize()
    else:
        nFloat = profile[ind].pdf.getParameters(hist_data).selectByAttrib('Constant',False).getSize()
    stat_vals.append(stat.getVal() + 0.5*nFloat)
    print (profile[ind].name,' = ', nFloat)
    print('NLL = ', stat.getVal() + PENALTY*nFloat)
for model in profile:
    model.checkBond()

best_ = stat_vals.index(min(stat_vals))
print ("best is ", profile[best_].name)

# Create Asimov (optional)
N = hist_data.sumEntries()
sig_model_el = DSCB_Class(x, MH, CAT+'_el', di_sigma = False)
sig_model_el.assignVal(args.configS, cat=CAT_merge, lep="el")
sig_model_mu = DSCB_Class(x, MH, CAT+'_mu', di_sigma = False)
sig_model_mu.assignVal(args.configS, cat=CAT_merge, lep="mu")
c_el = ROOT.RooRealVar('c_el', 'c_el', sig_model_el.nsig/(sig_model_el.nsig + sig_model_mu.nsig))
combine_model = ROOT.RooAddPdf('duo_sig_model_'+CAT, 'duo_sig_model_'+CAT, sig_model_el.pdf, sig_model_mu.pdf, c_el)
N_sig = sig_model_el.nsig + sig_model_mu.nsig
c1= ROOT.RooRealVar('c1','c1', N_sig, 0, 300)
c2 = ROOT.RooRealVar('c2','c2', N-N_sig, 0, 3*N)

'''
nll1 = ROOT.RooNLLVar('nll1', 'nll2',  profile[best_].pdf, hist_asimov)
r1 = Minimizer_NLL(nll1, -1, eps, True, strategy)
nll1_val = nll1.getVal()
nll2 = ROOT.RooNLLVar('nll1', 'nll2',  tot_model, hist_asimov)
r2 = Minimizer_NLL(nll2, -1, eps, True, strategy)
nll2_val = nll2.getVal()
'''

# Create background workspace 
cate = ROOT.RooCategory("pdfindex_"+CAT, "Index of Pdf which is active for "+CAT)
models = ROOT.RooArgList()
for model in profile:
    models.add(model.pdf)
multipdf = ROOT.RooMultiPdf("multipdf_"+CAT, "MultiPdf for "+CAT, cate, models)

# Penalty term
multipdf.setCorrectionFactor(PENALTY)

norm = ROOT.RooRealVar("multipdf_"+ CAT + "_norm", "Number of background events", N, 0, 3*N)
if args.asimov == 0:
    f_out2 = ROOT.TFile("workspaces/workspace_bkg_"+ args.name +"_" + CAT + ".root", "RECREATE")
else:
    f_out2 = ROOT.TFile("workspaces/workspace_bkg_profile_bias_asimov_" +args.name +"_" + CAT + ".root", "RECREATE")
w_bkg = ROOT.RooWorkspace("workspace_bkg","workspace_bkg")
getattr(w_bkg, "import")(cate)
getattr(w_bkg, "import")(norm)
getattr(w_bkg, "import")(multipdf)
getattr(w_bkg, "import")(hist_data)
if args.asimov == 1:
    for entry in profile:
        tot_model = ROOT.RooAddPdf('tot_pdf', 'tot_pdf', ROOT.RooArgList(combine_model, entry.pdf), ROOT.RooArgList(c1, c2))
        hist_asimov = tot_model.generateBinned(ROOT.RooArgSet(x), N, ROOT.RooFit.Asimov(True))
        hist_asimov.SetNameTitle('hist_' +  entry.name, 'hist_' + entry.name)
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


# Plot it
multiPlotClass(x, hist_data, profile, title="Profile_"+ args.name + "_" +CAT, output_dir="", sideBand=True, fitRange= 'left,right', best_index = best_, bestLabel = False, fullRatio = False, ratio_range = [0,2])
#profile_.write_config_file(cuthist, "CMSBias")

# Create toy histogram (depreciated)
bias = False
N_toy = 10

if bias:
    f_out3 = ROOT.TFile("~/EOS_space/toy_wsp/workspace_toy_" + CAT + ".root", "RECREATE")
    w_toy = ROOT.RooWorkspace("workspace_toy","workspace_toy")
    for entry in profile:
        x.setBins(nbins)
        for i in range(N_toy):
            hist_toy = entry.pdf.generateBinned(x, ROOT.RooFit.NumEvents(N))
            hist_toy.SetNameTitle("hist_"+entry.pdf.GetName()+"_"+str(i), "hist_"+entry.pdf.GetName()+"_"+str(i))
            getattr(w_toy, "import")(hist_toy)
    w_toy.Write()
    w_toy.Print()
    f_out3.Close()


