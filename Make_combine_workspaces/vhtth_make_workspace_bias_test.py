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
from APFUtilities import *

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
x = define_mlly() 

mu_gauss = ROOT.RooRealVar("mu_gauss", "always 0", 0.)
MH = ROOT.RooRealVar("MH", "MH", 125, 120., 130.)





#plotClass(x, hist_sig, sig_model.pdf, sig_model.pdf, "Signal_"+CAT, CMS = 'Simulation', output_dir="")
profile_ = profileClass(x, mu_gauss, CAT, '../Config/'+args.config)
profile = profile_.testSelection("postFT")

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

# Plot it
multiPlotClass(x, hist_data, profile, title="Profile_" +CAT, output_dir="", sideBand=True, fitRange= 'left,right', best_index = best_, bestLabel = True)

'''
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


# Create signal workspace
if SIGNAL:
    f_out1 = ROOT.TFile("workspaces/workspace_sig_" + CAT + ".root", "RECREATE")
    w_sig = ROOT.RooWorkspace("workspace_sig","workspace_sig")
    getattr(w_sig, "import")(sig_model.pdf)
    w_sig.Print()
    w_sig.Write()
    f_out1.Close()
'''



'''
# Create Asimov (optional)
c1= ROOT.RooRealVar('c1','c1', N_sig)#, 0, 3*N_sig)
c2 = ROOT.RooRealVar('c2','c2', N)#, 0, 3*N)
tot_model = ROOT.RooAddPdf('tot_pdf', 'tot_pdf', ROOT.RooArgList(sig_model.pdf, profile[best_].pdf), ROOT.RooArgList(c1, c2))
hist_asimov = tot_model.generateBinned(ROOT.RooArgSet(x), N, ROOT.RooFit.Asimov(True))
hist_asimov.SetNameTitle('hist_' + CAT + '_asimov', 'hist_' + CAT + '_asimov')
'''
