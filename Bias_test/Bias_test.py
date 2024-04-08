import ROOT
import os
import sys
import argparse
sys.path.append(os.path.abspath("../Utilities/"))
sys.path.append(os.path.abspath("../CMS_plotter/"))
import CMS_lumi, tdrstyle
from bkg_functions_fit import *
from bkg_functions_class import *
from sig_functions_class import *
from Xc_Minimizer import *
from plot_utility import *
from sample_reader import *
from bias_class import *
ROOT.gInterpreter.AddIncludePath('../Utilities/HZGRooPdfs.h')
ROOT.gSystem.Load('../Utilities/HZGRooPdfs_cxx.so')
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)

# Define variables
lowx = 105.
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
"""
N = ROOT.RooRealVar("N", "Extended term", 19000, 10000, 30000);
dummy_N = ROOT.RooRealVar("dummy_N", "dummy N", 10, 0, 1000);
dummy_mu = ROOT.RooRealVar("d_mu", "a dummy mu", 125)
dummy_width = ROOT.RooRealVar("d_sigma", "a dummy width", 0.1)
dummy_sig =  ROOT.RooGaussian("dummy", "dummy for ext", x,dummy_mu, dummy_width)
"""

# Read samples
#reader = readDat(list, "/afs/cern.ch/user/f/fanx/public/samples/")
#reader = readDat(list, "../../sample/")
reader = readRoot(x, "~/beijing_sample/data.root")
#cutHist = reader.data_hist_bin1.reduce(ROOT.RooFit.CutRange('left,right'))
#print ("cut norm = ", cutHist.sumEntries())
N = reader.data_hist_bin1.sumEntries() *10
print("10 times stats = ",N)
x.setBins(260)
x.setRange('left', lowx, 120)
x.setRange('right', 130, lowx+65)
x.setRange('full', lowx, lowx+65)

# Signal model (pre-fit)
MH = ROOT.RooRealVar("MH","MH"       ,124.7, 120., 130.)
dscb_model = DSCB_Class(x, MH, "bin1", 1.78,50, 100, 50, 100,0.845, 2.36)
dscb_model.setConst(True)
N_sig = 10

# Functions to test
bern2_model_seed = Bern2Class(x, mu_gauss, "bin1_seed", 10, 0.3, 10, 7., 105.)
bern3_model_seed = Bern3Class(x, mu_gauss, "bin1_seed", 10, 0.3, 10, 3., 106.)
bern4_model_seed = Bern4Class(x, mu_gauss, "bin1_seed", 10, 0.3, 10, 3., 106.)

bern2_model = Bern2Class(x, mu_gauss, "bin1", 10, 0.3, 10, 7., 105.)
bern3_model = Bern3Class(x, mu_gauss, "bin1", 10, 0.3, 10, 3., 106.)
bern4_model = Bern4Class(x, mu_gauss, "bin1", 10, 0.3, 10, 3., 106.)
#extmodel = ROOT.RooExtendPdf("extmodel", "Extended model", bern2_model.pdf, N, 'full');
#extmodel = ROOT.RooAddPdf("extmodel","extmodel", ROOT.RooArgList(bern3_model.pdf, dummy_sig), ROOT.RooArgList(N, dummy_N))
#r = bern2_model.pdf.fitTo(reader.data_hist_untagged1_bkg,ROOT.RooFit.Save(True), ROOT.RooFit.PrintLevel(-1), ROOT.RooFit.SumW2Error(True))
#r.Print("v")
#ROOT.RooFit.Range('left,right'),
profile_seed = [bern2_model_seed, bern3_model_seed, bern4_model_seed]
profile = [bern2_model, bern3_model, bern4_model]

# Set best-fit values
for entry in profile_seed:
    BiasClass(entry.pdf, reader.data_hist_bin1, True, 'left,right').minimize()
    entry.checkBond()

print("Done seed PDFs")

# Define profile Fit
def profilefFit(profile, sig_model, hist, fix = False, str = 0.):
    min_nll = 0
    ind = 999
    r_sig_ = 0
    r_error_ = 0
    best_=''
    for i, ele in enumerate(profile):
        c1 = ROOT.RooRealVar("c1", "c1", N, 0, 3.* N)
        if (fix): c2 = ROOT.RooRealVar("c2", "c2", str)
        else: c2 = ROOT.RooRealVar("c2", "c2", 0., -100.*N_sig, 100.*N_sig)
        tot_model = ROOT.RooAddPdf("tot_model", "tot_model", ROOT.RooArgList(sig_model.pdf, ele.pdf), ROOT.RooArgList(c2, c1))
        bias = BiasClass(tot_model, hist, False)
        bias.minimize()
        if i ==0: 
            ind = 0
            min_nll=bias.corrNLL
            r_sig_ = c2.getVal()
            r_error_ = c2.getError()
            best_= ele.pdf.GetName()
        elif bias.corrNLL< min_nll: 
            ind = i
            min_nll = bias.corrNLL
            r_sig_ = c2.getVal()
            r_error_ = c2.getError()
            best_= ele.pdf.GetName()
    return [ind, min_nll, r_sig_, r_error_, best_]

# Discrete profiling - Find minimum and (r_down, r_up)
N_toy = 2
N_scan = 30
for entry in profile_seed:
    r_sig = []
    r_error = []
    best_list = []
    for j in range(N_toy):      
        x.setBins(260)
        hist_toy = entry.pdf.generateBinned(x, ROOT.RooFit.NumEvents(N))
        list = profilefFit(profile, dscb_model, hist_toy)
        r_sig.append(list[2]/N_sig)
        #r_error.append(list[3]/N_sig)
        best_list.append(list[4])
        # tot_model_ = ROOT.RooAddPdf("tot_model_", "tot_model_", ROOT.RooArgList(dscb_model.pdf, profile[ind].pdf), ROOT.RooArgList(c2, c1))
        # BiasClass(tot_model_, hist_toy, False).minimize()
        # plotClass(x, hist_toy, tot_model_, title = entry.pdf.GetName(), sideBand = False)
        NLL_list = []
        for k in range(N_scan):
            list_ = profilefFit(profile, dscb_model, hist_toy, True, list[2] * (k - N_scan/2))
            NLL_list.append(list_[1])
        dNLL = [x - list[1] for x in NLL_list]
        left = 0
        right = 0
        for i in range(len(dNLL) - 1):
            if dNLL[i] > 0.5 and dNLL[i+1] < 0.5: left = i
            if dNLL[i] < 0.5 and dNLL[i+1] > 0.5: right = i
        if left == right: print("Scan error! in ", entry.pdf.GetName())
        




    print("r = ", r_sig)
    print("r error = ", r_error)
    print("best func = ", best_list)
    print("Finish ", entry.pdf.GetName()," toy sample")






#print ("nll = ", bern2Bias.nll.getVal())
#print ("corr nll = ", bern2Bias.corrNLL)
#print ("norm = ", extmodel.getNorm())
#i1 = bern2_model.pdf.createIntegral(x, ROOT.RooFit.Range('left,right'))
#i2 = bern2_model.pdf.createIntegral(x, ROOT.RooFit.Range('full'))
#print ("full int = ", i2.getVal(), " side int = ", i1.getVal())

#plotClass(x, reader.data_hist_bin1, bern2_model.pdf, title = "Bern2", sideBand = True, fitRange = 'left,right')
#plotClass(x, reader.data_hist_bin1, bern3_model.pdf, title = "Bern3", sideBand = True, fitRange = 'left,right')
#plotClass(x, reader.data_hist_bin1, bern4_model.pdf, title = "Bern4", sideBand = True, fitRange = 'left,right')
