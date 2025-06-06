#!/usr/bin/env python3
import ROOT
import os
import sys
import json
#import matplotlib.pyplot as plt
import numpy as np
import argparse
sys.path.append(os.path.abspath("../Utilities/"))
sys.path.append(os.path.abspath("../CMS_plotter/"))
import CMS_lumi, tdrstyle
from bkg_functions_class import *
from sig_functions_class import *
from Xc_Minimizer import *
from plot_utility import *
from sample_reader import *
from bias_class import *
from profile_class import *
from sig_functions_class import *
import CMS_lumi, tdrstyle

ROOT.gInterpreter.Declare("""
RooDataHist readHistWs(TString filename, string CAT, string lep){
auto file_ = new TFile(filename, "READ");
RooWorkspace *ws = (RooWorkspace *)file_->Get(("workspace_bkg_"+lep).c_str());
std::cout << "Read the ws well" << std::endl;
RooDataHist *hist  = (RooDataHist *)ws->data(("hist_" + CAT +"_"+lep).c_str());
return *hist;
}
""")


file1 = 'workspaces/workspace_bkg_toy_unblind_vbf1.root'
file2 = 'workspaces/workspace_bkg_toy_unblind_vbf2.root'

hist1_el = ROOT.readHistWs(file1, 'vbf1', 'el')
hist1_mu = ROOT.readHistWs(file1, 'vbf1', 'mu')
hist2_el = ROOT.readHistWs(file2, 'vbf2', 'el')
hist2_mu = ROOT.readHistWs(file2, 'vbf2', 'mu')

print('n vbf1 el = ', hist1_el.sumEntries())
print('n vbf2 el = ', hist2_el.sumEntries())

x = ROOT.RooRealVar('CMS_hzg_mass_vbf12','CMS_hzg_mass_vbf12', 96, 161)
mu_gauss = ROOT.RooRealVar("mu_gauss","always 0"       ,0.)
x.setBins(260)
hist_el = ROOT.RooDataHist('hist_vbf12_el', 'hist_vbf12_el', x)
hist_mu = ROOT.RooDataHist('hist_vbf12_mu', 'hist_vbf12_mu', x)

for i in range(260):
    hist2_el.get(i)
    hist2_mu.get(i)
    hist1_el.get(i+4)
    hist1_mu.get(i+4)
    hist_el.get(i)
    hist_mu.get(i)
    if i < 256:
        hist_el.set(hist1_el.weight() + hist2_el.weight())
        hist_mu.set(hist1_mu.weight() + hist2_mu.weight())
    else:
        hist_el.set(hist2_el.weight())
        hist_mu.set(hist2_mu.weight())
print('n vbf12 el = ', hist_el.sumEntries())

'''
frame = x.frame()
hist_el.plotOn(frame)
c = ROOT.TCanvas('c','c', 500, 500)
frame.Draw()
c.SaveAs('test.pdf')
'''

profile_ = profileClass(x, mu_gauss, 'vbf2', '../Config/chi2_config_xgboost_nodrop.json')
profile = profile_.testSelection("FT")

cate_el = ROOT.RooCategory("pdfindex_vbf12_el", "Index of Pdf which is active for vbf12 el")
cate_mu = ROOT.RooCategory("pdfindex_vbf12_mu", "Index of Pdf which is active for vbf12 mu")
models = ROOT.RooArgList()
for model in profile:
    models.add(model.pdf)

multipdf_el = ROOT.RooMultiPdf("multipdf_vbf12_el", "MultiPdf for vbf12 el", cate_el, models)
multipdf_mu = ROOT.RooMultiPdf("multipdf_vbf12_mu", "MultiPdf for vbf12 mu", cate_mu, models)
norm_el = ROOT.RooRealVar("multipdf_vbf12_el_norm", "Number of background events el vbf12", hist_el.sumEntries(), 0, 3*hist_el.sumEntries())
norm_mu = ROOT.RooRealVar("multipdf_vbf12_mu_norm", "Number of background events mu vbf12", hist_mu.sumEntries(), 0, 3*hist_mu.sumEntries())

f_out = ROOT.TFile("workspaces/workspace_bkg_toy_unblind_vbf12.root", "RECREATE")
w_bkg_el =  ROOT.RooWorkspace("workspace_bkg_el","workspace_bkg_el")
w_bkg_mu = ROOT.RooWorkspace("workspace_bkg_mu","workspace_bkg_mu")
getattr(w_bkg_mu, "import")(cate_mu)
getattr(w_bkg_mu, "import")(norm_mu)
getattr(w_bkg_mu, "import")(multipdf_mu)
getattr(w_bkg_mu, "import")(hist_mu)

getattr(w_bkg_el, "import")(cate_el)
getattr(w_bkg_el, "import")(norm_el)
getattr(w_bkg_el, "import")(multipdf_el)
getattr(w_bkg_el, "import")(hist_el)

w_bkg_el.Write()
w_bkg_mu.Write()
f_out.Close()

MH = ROOT.RooRealVar("MH","MH"       ,125)
sig_model1_el = combineSignalLep(x, MH, 'vbf1', '../Config/config_DSCB.json', 'el')
sig_model1_mu = combineSignalLep(x, MH, 'vbf1', '../Config/config_DSCB.json', 'mu')
sig_model2_el = combineSignalLep(x, MH, 'vbf2', '../Config/config_DSCB.json', 'el')
sig_model2_mu = combineSignalLep(x, MH, 'vbf2', '../Config/config_DSCB.json', 'mu')
c_sig_el = ROOT.RooRealVar('c_sig_el', 'c_sig_el', sig_model1_el.ntot/(sig_model1_el.ntot+sig_model2_el.ntot))
c_sig_mu = ROOT.RooRealVar('c_sig_mu', 'c_sig_mu', sig_model1_mu.ntot/(sig_model1_mu.ntot+sig_model2_mu.ntot))
sig_model12_el = ROOT.RooAddPdf('sig_model12_el', 'sig_model12_el', ROOT.RooArgList(sig_model1_el.pdf, sig_model2_el.pdf), c_sig_el)
sig_model12_mu = ROOT.RooAddPdf('sig_model12_mu', 'sig_model12_mu', ROOT.RooArgList(sig_model1_mu.pdf, sig_model2_mu.pdf), c_sig_mu)

f_out1 = ROOT.TFile("workspaces/workspace_sig_toy_unblind_vbf12.root", "RECREATE")
w_sig_el = ROOT.RooWorkspace("workspace_sig_el","workspace_sig_el")
getattr(w_sig_el, "import")(sig_model12_el)
w_sig_el.Write()
w_sig_mu = ROOT.RooWorkspace("workspace_sig_mu","workspace_sig_mu")
getattr(w_sig_mu, "import")(sig_model12_mu)
w_sig_mu.Write()
f_out1.Close()
