#!/usr/bin/env python3
import sys
import os
import ROOT
import argparse
import json
sys.path.append(os.path.abspath("../Utilities/"))
from Xc_Minimizer import *

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--cat', help = 'Category')
parser.add_argument('-con', '--config', help = 'Configuration', default='../Config/chi2_config_xgboost_nodrop.json')
parser.add_argument('-f', '--func', help = 'Function')
parser.add_argument('-i', '--input', help = 'Input Root File')
parser.add_argument('-l', '--lep', help = 'Lepton', default='')

args = parser.parse_args()
CAT = args.cat
jfile = open(args.config, 'r')
configs = json.load(jfile)
setting = configs[CAT]

f = ROOT.TFile(args.input)
w = f.Get("w")
lowx = setting["Range"][0]
highx = setting["Range"][1]
nbins = int(setting["Bins"])
binning = ROOT.RooFit.Binning(nbins,lowx,highx)

x = w.var('CMS_hzg_mass_'+CAT)
x.setRange('left', lowx, 120)
x.setRange('right', 130, highx)

plot = x.frame()
plot1 = x.frame()

plot.SetTitle(CAT + ' Unblind Stage 1 '+args.lep)
plot1.SetTitle(CAT + ' Unblind Stage 2 '+args.lep)

data = w.data("data_obs")
print(type(data))
ntot = data.sumEntries()
dataSB = data.reduce(ROOT.RooFit.CutRange('left,right'))
dataSB.plotOn( plot, binning)
data.plotOn( plot1, binning)
nSB = dataSB.sumEntries()
nbkg_prefit =  w.function('n_exp_final_bin'+CAT+'_proc_bkg').getVal()
print("n tot SB= ", nSB, "ntot = ", ntot)

# Toy sideband fit
b_model = w.pdf(args.func+'_'+CAT+'_model')
SBpdf = ROOT.RooGenericPdf("toy_SB_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, b_model))
nll = ROOT.RooNLLVar('nll_SB', 'nll_SB', SBpdf, dataSB)
Minimizer_NLL(nll, -1, 1, False, 0)
res=Minimizer_NLL(nll, -1, 0.1, True, 0)
res.Print('v')
dof = nbins - res.floatParsFinal().getSize()

# Prefit
b_model.plotOn( plot, ROOT.RooFit.LineColor(2), ROOT.RooFit.Name("prefit"), ROOT.RooFit.Normalization(nbkg_prefit/nSB))
b_model.plotOn( plot1, ROOT.RooFit.LineColor(2), ROOT.RooFit.Name("prefit"), ROOT.RooFit.Normalization(nbkg_prefit/ntot))

# Postfit
w.loadSnapshot("MultiDimFit")
nbkg = w.function('n_exp_final_bin'+CAT+'_proc_bkg').getVal()
print("post fit nbkg = ", nbkg)
x.setBins(nbins)
hdata = ROOT.RooDataHist('hdata', 'hdata', x, data)
chi2 = ROOT.RooChi2Var('chi2', 'chi2', b_model, hdata)
val_chi2 = chi2.getVal()
x.setBins(int(highx - lowx))
latex = ROOT.TLatex()
latex.SetTextSize(0.03)
latex.SetTextAlign(13)

b_model.plotOn( plot, ROOT.RooFit.LineColor(4), ROOT.RooFit.Name("postfit"), ROOT.RooFit.Normalization(nbkg/nSB))
b_model.plotOn( plot1, ROOT.RooFit.LineColor(4), ROOT.RooFit.Name("postfit"), ROOT.RooFit.Normalization(nbkg/ntot))

can = ROOT.TCanvas()
plot.Draw()
leg = ROOT.TLegend(0.55,0.6,0.85,0.85)
leg.AddEntry("prefit", "Prefit B-only model", "L")
leg.AddEntry("postfit", "Postfit B-only model", "L")
leg.Draw("Same")
latex.DrawLatexNDC(.6, .5, '#Chi^{2}/dof = %.2f'%val_chi2+ '/%i'%dof)
can.Update()
can.SaveAs('plots/unblind_stage1_'+CAT+args.lep+'.png')

can1 = ROOT.TCanvas()
can1.cd()
plot1.Draw()
leg.Draw("Same")
can1.Update()
can1.SaveAs('plots/unblind_stage2_'+CAT+args.lep+'.png')

