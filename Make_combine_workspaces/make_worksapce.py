import ROOT
import os
import sys
import argparse
sys.path.append(os.path.abspath("../Utilities/"))
from bkg_functions_class import *
from Xc_Minimizer import *
from plot_utility import *
from sample_reader import *
ROOT.gInterpreter.AddIncludePath('../Utilities/HZGRooPdfs.h')
ROOT.gInterpreter.AddIncludePath('../Utilities/ModGaus.h')
ROOT.gSystem.Load('../Utilities/HZGRooPdfs_cxx.so')
ROOT.gSystem.Load('../Utilities/ModGaus_cxx.so')

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

# Cornell MC sample dat reader and make RooDataHist
x.setBins(260)
reader = readDat(list, "../../sample/")
reader.numCheck()

# Assume we have Bern2, Bern3, Pow1, Exp1, Exp2, Lau1, Lau2 in the profile, the signal model is DSCB

bern2_u1 = Bern2Class(x, gass_mu, "u1")
bern3_u1 = Bern2Class(x, gass_mu, "u1")
pow1_u1 = Pow1Class(x, gauss_mu, "u1")
exp1_u1 = Exp1Class(x, gauss_mu, "u1")
exp2_u1 = Exp2Class(x, gauss_mu, "u1")
lau1_u1 = Lau1Class(x, gass_mu, "u1")
lau2_u1 = Lau2Class(x, gass_mu, "u1")

sig_u1 = DSCB_Class(x, MH, "u1")



