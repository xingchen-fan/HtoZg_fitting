#!/usr/bin/env python3
import sys
import os
import ROOT
import argparse
import json
sys.path.append(os.path.abspath("../Utilities/"))
from Xc_Minimizer import *
ROOT.gStyle.SetOptStat(0)

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--cat', help = 'Category')
parser.add_argument('-con', '--config', help = 'Configuration', default='')
parser.add_argument('-f1', '--func1', help = 'Function prefit')
parser.add_argument('-f2', '--func2', help = 'Function postfit')
parser.add_argument('-i', '--input', help = 'Input Root File')
parser.add_argument('-t', '--toy', help = 'Toy?', default=False)
parser.add_argument('-d', '--data', help = 'Data/toy histogram input', default='')
parser.add_argument('-l', '--lep', help = 'Lepton', default='')
parser.add_argument('-s', '--stage', help = 'Stage', default='1')

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
binning = ROOT.RooFit.Binning(65,lowx,highx)

x = w.var('CMS_hzg_mass_'+CAT)
x.setRange('left', lowx, 120)
x.setRange('right', 130, highx)

plot = x.frame()
plot.SetTitle(CAT + ' Unblind Stage ' +args.stage + ' '+args.lep)
ROOT.gStyle.SetTitleSize(0.08, "t")
plot.GetYaxis().SetTitleSize(0.05)
plot.GetYaxis().SetTitleOffset(0.6)
plot.GetYaxis().SetLabelSize(0.05)
plot.GetYaxis().SetTitle("Events/(1 GeV)")

print(type(data))
ntot = data.sumEntries()
dataSB = data.reduce(ROOT.RooFit.CutRange('left,right'))

nSB = dataSB.sumEntries()
nbkg_prefit =  w.function('n_exp_final_bincat_'+CAT+'_proc_background').getVal()
print("Prefit norm = ",nbkg_prefit)
print("n tot SB= ", nSB, "ntot = ", ntot)
norm_ = 0
if args.stage == '1':
    norm_ = nSB
    show_hist = dataSB
elif args.stage == '2':
    norm_ = ntot
    show_hist = data
else:
    raise Exception("Wrong stage.")
show_hist.plotOn( plot, binning)

# Prefit
b_model_1 = w.pdf(args.func1+'_'+CAT+'_model')
b_model_1.plotOn( plot, ROOT.RooFit.LineColor(2), ROOT.RooFit.Name("prefit"), ROOT.RooFit.Normalization(ntot/norm_))

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
b_model_2.plotOn( plot, ROOT.RooFit.LineColor(4), ROOT.RooFit.Name("postfit"), ROOT.RooFit.Normalization(nbkg/norm_))

can = ROOT.TCanvas("c","c", 700, 500)
can.Divide(1,2)
can.cd(1)
ROOT.gPad.SetPad(0.0, 0.35, 1.0, 0.99)
ROOT.gPad.SetBottomMargin(0)
ROOT.gPad.SetTopMargin(0.15)

plot.Draw()
leg = ROOT.TLegend(0.55,0.6,0.85,0.85)
leg.AddEntry("prefit", "Prefit B-only model", "L")
leg.AddEntry("postfit", "Postfit B-only model", "L")
leg.Draw("Same")
latex.DrawLatexNDC(.6, .55, '#Chi^{2}/dof = %.2f'%val_chi2+ '/%i'%dof)
ROOT.gPad.Update()
if args.stage == '1':
    box = ROOT.TBox(120,0, 130, ROOT.gPad.GetFrame().GetY2())
    box.SetFillColorAlpha(16, 0.4)
    box.Draw("same")

can.cd(2)
ROOT.gPad.SetPad(0.0, 0.1, 1.0, 0.35)
ROOT.gPad.SetTopMargin(0)
ROOT.gPad.SetBottomMargin(0.38)
h_hist = show_hist.createHistogram("h_hist", x, ROOT.RooFit.Binning(65))
model_hist = b_model_2.generateBinned(x, show_hist.sumEntries(), True).createHistogram("model_hist", x, ROOT.RooFit.Binning(65))
for i in range(0, model_hist.GetNbinsX() + 2):
    model_hist.SetBinError(i, 0.0)
h_hist.Add(model_hist, -1)
h_hist.SetMarkerColor(ROOT.kBlack)
h_hist.SetMarkerStyle(8)
h_hist.SetMarkerSize(1)
h_hist.SetLineColor(ROOT.kBlack)
h_hist.SetTitle("")
h_hist.GetXaxis().SetTitleOffset(0.9)
h_hist.GetXaxis().SetTitleSize(0.12)
h_hist.GetXaxis().SetTitle("m_{ll\gamma}/(GeV)")
h_hist.GetXaxis().SetLabelSize(0.12)
h_hist.GetYaxis().SetLabelSize(0.10)
h_hist.GetYaxis().SetTitleOffset(0.2)
h_hist.GetYaxis().SetTitleSize(0.12)
h_hist.GetYaxis().SetNdivisions(4)
h_hist.GetYaxis().SetTitle("Data - Postfit")
minimum = 99999
for i in range(int(x.getMax() - x.getMin())):
    bincont = h_hist.GetBinContent(i+1)
    if bincont != 0 and bincont < minimum:
        minimum = bincont
range_ = 1.1 * max(abs(minimum), abs(h_hist.GetBinContent(h_hist.GetMaximumBin())))
h_hist.GetYaxis().SetRangeUser(-range_, range_)

line = ROOT.TLine( x.getMin(), 0, x.getMax(), 0)
line.SetLineColor(ROOT.kBlack)
line.SetLineStyle(7)
line.SetLineWidth(2)
h_hist.Draw()
line.Draw('same')
ROOT.gPad.Update()
if args.stage == '1':
    box2 = ROOT.TBox(120,ROOT.gPad.GetFrame().GetY1(), 130, ROOT.gPad.GetFrame().GetY2())
    box2.SetFillColorAlpha(16, 0.4)
    box2.Draw("same")

can.Update()
can.SaveAs('plots/unblind_stage'+args.stage+'_'+CAT+args.lep+'.png')

