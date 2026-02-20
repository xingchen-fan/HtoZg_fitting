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
parser.add_argument('-con', '--config', help = 'Configuration', default='')
parser.add_argument('-f1', '--func1', help = 'Function prefit')
parser.add_argument('-f2', '--func2', help = 'Function postfit')
parser.add_argument('-i', '--input', help = 'Input Root File')
parser.add_argument('-t', '--toy', help = 'Toy?', default=False)
parser.add_argument('-d', '--data', help = 'Data/toy histogram input', default='')
parser.add_argument('-l', '--lep', help = 'Lepton', default='')

args = parser.parse_args()
CAT = args.cat
jfile = open(args.config, 'r')
configs = json.load(jfile)
setting = configs[CAT]

if CAT == 'vhptmiss':
    CAT_m = 'vhmet'
elif CAT == 'untag':
    CAT_m = 'untagged'
else:
    CAT_m = CAT
    
f = ROOT.TFile(args.input)
w = f.Get("w")
f_h = ROOT.TFile(args.data)
w_h = f_h.Get("workspace_bkg")
if args.toy:
    data = w_h.data("hist_toy5_"+CAT)
else:
    data = w_h.data("data_obs_cat_"+CAT_m)
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

print(type(data))
ntot = data.sumEntries()
dataSB = data.reduce(ROOT.RooFit.CutRange('left,right'))
dataSB.plotOn( plot, binning)
data.plotOn( plot1, binning)
nSB = dataSB.sumEntries()
nbkg_prefit =  w.function('n_exp_final_bincat_'+CAT+'_proc_background').getVal()
print("Prefit norm = ",nbkg_prefit)
print("n tot SB= ", nSB, "ntot = ", ntot)

# Prefit
b_model_1 = w.pdf(args.func1+'_'+CAT+'_model')
b_model_1.plotOn( plot, ROOT.RooFit.LineColor(2), ROOT.RooFit.Name("prefit"), ROOT.RooFit.Normalization(ntot/nSB))
b_model_1.plotOn( plot1, ROOT.RooFit.LineColor(2), ROOT.RooFit.Name("prefit"), ROOT.RooFit.Normalization(ntot/ntot))

# Postfit
b_model_2 = w.pdf(args.func2+'_'+CAT+'_model')
w.loadSnapshot("MultiDimFit")
nbkg = w.function('n_exp_final_bincat_'+CAT+'_proc_background').getVal()
print("post fit norm = ", nbkg)
x.setBins(nbins)
hdata = ROOT.RooDataHist('hdata', 'hdata', x, data)
chi2 = ROOT.RooChi2Var('chi2', 'chi2', b_model_2, hdata)
val_chi2 = chi2.getVal()
x.setBins(int(highx - lowx))
latex = ROOT.TLatex()
latex.SetTextSize(0.04)
latex.SetTextAlign(13)
pars = b_model_2.getParameters(x)
n_float = 0
for p in pars:
    if isinstance(p, ROOT.RooRealVar) and not p.isConstant():
        n_float += 1
print("nfloat = ", n_float)
dof = nbins - n_float
b_model_2.plotOn( plot, ROOT.RooFit.LineColor(4), ROOT.RooFit.Name("postfit"), ROOT.RooFit.Normalization(nbkg/nSB))
b_model_2.plotOn( plot1, ROOT.RooFit.LineColor(4), ROOT.RooFit.Name("postfit"), ROOT.RooFit.Normalization(nbkg/ntot))

can = ROOT.TCanvas()
plot.Draw()
leg = ROOT.TLegend(0.55,0.6,0.85,0.85)
leg.AddEntry("prefit", "Prefit B-only model", "L")
leg.AddEntry("postfit", "Postfit B-only model", "L")
leg.Draw("Same")
latex.DrawLatexNDC(.6, .55, '#Chi^{2}/dof = %.2f'%val_chi2+ '/%i'%dof)
can.Update()
can.SaveAs('plots/unblind_stage1_'+CAT+args.lep+'.png')

can1 = ROOT.TCanvas()
can1.cd()
plot1.Draw()
leg.Draw("Same")
can1.Update()
can1.SaveAs('plots/unblind_stage2_'+CAT+args.lep+'.png')

