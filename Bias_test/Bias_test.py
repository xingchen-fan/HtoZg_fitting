import ROOT
import os
import sys
import argparse
sys.path.append(os.path.abspath("../Utilities/"))
sys.path.append(os.path.abspath("../CMS_plotter/"))
import CMS_lumi, tdrstyle
from bkg_functions_fit import *
from bkg_functions_class import *
from Xc_Minimizer import *
from plot_utility import *
from sample_reader import *
from bias_class import *
ROOT.gInterpreter.AddIncludePath('../Utilities/HZGRooPdfs.h')
ROOT.gSystem.Load('../Utilities/HZGRooPdfs_cxx.so')

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

"""
N = ROOT.RooRealVar("N", "Extended term", 19000, 10000, 30000);
dummy_N = ROOT.RooRealVar("dummy_N", "dummy N", 10, 0, 1000);
dummy_mu = ROOT.RooRealVar("d_mu", "a dummy mu", 125)
dummy_width = ROOT.RooRealVar("d_sigma", "a dummy width", 0.1)
dummy_sig =  ROOT.RooGaussian("dummy", "dummy for ext", x,dummy_mu, dummy_width)
"""
x.setBins(260)
#reader = readDat(list, "/afs/cern.ch/user/f/fanx/public/samples/")
#reader = readDat(list, "../../sample/")
reader = readRoot(x, "~/beijing_sample/data.root")
#reader.numCheck()

x.setRange('left', lowx, 120)
x.setRange('right', 130, lowx+65)
x.setRange('full', lowx, lowx+65)




#cutHist = reader.data_hist_bin1.reduce(ROOT.RooFit.CutRange('left,right'))
#print ("cut norm = ", cutHist.sumEntries())
bern2_model = Bern2Class(x, mu_gauss, "bin1", 10, 0.3, 10, 7., 105.)
bern3_model = Bern3Class(x, mu_gauss, "bin1", 10, 0.3, 10, 3., 106.)
bern4_model = Bern4Class(x, mu_gauss, "bin1", 10, 0.3, 10, 3., 106.)



#extmodel = ROOT.RooExtendPdf("extmodel", "Extended model", bern2_model.pdf, N, 'full');
#extmodel = ROOT.RooAddPdf("extmodel","extmodel", ROOT.RooArgList(bern3_model.pdf, dummy_sig), ROOT.RooArgList(N, dummy_N))
#r = bern2_model.pdf.fitTo(reader.data_hist_untagged1_bkg,ROOT.RooFit.Save(True), ROOT.RooFit.PrintLevel(-1), ROOT.RooFit.SumW2Error(True))
#r.Print("v")
#ROOT.RooFit.Range('left,right'),


bern2Bias = BiasClass(bern2_model.pdf, reader.data_hist_bin1, True, 'left,right')
bern3Bias = BiasClass(bern3_model.pdf, reader.data_hist_bin1, True, 'left,right')
bern4Bias = BiasClass(bern4_model.pdf, reader.data_hist_bin1, True, 'left,right')

bern2Bias.minimize()
bern2_model.checkBond()
bern3Bias.minimize()
bern3_model.checkBond()
bern4Bias.minimize()
bern4_model.checkBond()

#print ("nll = ", bern2Bias.nll.getVal())
#print ("corr nll = ", bern2Bias.corrNLL)
#print ("norm = ", extmodel.getNorm())
#i1 = bern2_model.pdf.createIntegral(x, ROOT.RooFit.Range('left,right'))
#i2 = bern2_model.pdf.createIntegral(x, ROOT.RooFit.Range('full'))
#print ("full int = ", i2.getVal(), " side int = ", i1.getVal())

plotClass(x, reader.data_hist_bin1, bern2_model.pdf, title = "Bern2", sideBand = True, fitRange = 'left,right')
plotClass(x, reader.data_hist_bin1, bern3_model.pdf, title = "Bern3", sideBand = True, fitRange = 'left,right')
plotClass(x, reader.data_hist_bin1, bern4_model.pdf, title = "Bern4", sideBand = True, fitRange = 'left,right')
