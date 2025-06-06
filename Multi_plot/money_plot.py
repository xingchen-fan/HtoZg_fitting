#!/usr/bin/env python3
import sys
import os
import ROOT
import argparse
import json
import math
sys.path.append(os.path.abspath("../Utilities/"))
sys.path.append(os.path.abspath("../CMS_plotter/"))
import CMS_lumi, tdrstyle
from Xc_Minimizer import *
from plot_utility import *
from profile_class import *
from sig_functions_class import *
import numpy as np
ROOT.gInterpreter.AddIncludePath('../Utilities/RooGaussStepBernstein.h')
ROOT.gSystem.Load('../Utilities/RooGaussStepBernstein_cxx.so')

parser = argparse.ArgumentParser()
parser.add_argument('-conB', '--configB', help = 'Bkg Configuration')
parser.add_argument('-conS', '--configS', help = 'Sig Configuration')

args = parser.parse_args()
lowx = 110
highx = 160
nbins = 50
binning = ROOT.RooFit.Binning(nbins,lowx,highx)

x = ROOT.RooRealVar('x', 'x', lowx, highx)
mu_gauss = ROOT.RooRealVar("mu_gauss","always 0"       ,0.)
MH = ROOT.RooRealVar("MH","MH"       ,125)

combine_model_ggf1 = combineSignal(x, MH, 'ggf1', args.configS)
combine_model_ggf2 = combineSignal(x, MH, 'ggf2', args.configS)
combine_model_ggf3 = combineSignal(x, MH, 'ggf3', args.configS)
combine_model_ggf4 = combineSignal(x, MH, 'ggf4', args.configS)
combine_model_vbf1 = combineSignal(x, MH, 'vbf1', args.configS)
combine_model_vbf2 = combineSignal(x, MH, 'vbf2', args.configS)
combine_model_vbf3 = combineSignal(x, MH, 'vbf3', args.configS)
combine_model_vbf4 = combineSignal(x, MH, 'vbf4', args.configS)
nggf1 = combine_model_ggf1.ntot
nggf2 = combine_model_ggf2.ntot
nggf3 = combine_model_ggf3.ntot
nggf4 = combine_model_ggf4.ntot
nvbf1 = combine_model_vbf1.ntot
nvbf2 = combine_model_vbf2.ntot
nvbf3 = combine_model_vbf3.ntot
nvbf4 = combine_model_vbf4.ntot
nSigExp = [nggf1, nggf2, nggf3, nggf4, nvbf1, nvbf2, nvbf3, nvbf4]
ntot = sum(nSigExp)
wc1 = ROOT.RooRealVar('wc1', 'wc1', nggf1/ntot)
wc2 = ROOT.RooRealVar('wc2', 'wc2', nggf2/ntot)
wc3 = ROOT.RooRealVar('wc3', 'wc3', nggf3/ntot)
wc4 = ROOT.RooRealVar('wc4', 'wc4', nggf4/ntot)
wc5 = ROOT.RooRealVar('wc5', 'wc5', nvbf1/ntot)
wc6 = ROOT.RooRealVar('wc6', 'wc6', nvbf2/ntot)
wc7 = ROOT.RooRealVar('wc7', 'wc7', nvbf3/ntot)
wc8 = ROOT.RooRealVar('wc8', 'wc8', nvbf4/ntot)

bkg_ggf1 = profileClass(x, mu_gauss, 'ggf1', args.configB).testSelection("Best")[0]
bkg_ggf2 = profileClass(x, mu_gauss, 'ggf2', args.configB).testSelection("Best")[0]
bkg_ggf3 = profileClass(x, mu_gauss, 'ggf3', args.configB).testSelection("Best")[0]
bkg_ggf4 = profileClass(x, mu_gauss, 'ggf4', args.configB).testSelection("Best")[0]
bkg_vbf1 = profileClass(x, mu_gauss, 'vbf1', args.configB).testSelection("Best")[0]
bkg_vbf2 = profileClass(x, mu_gauss, 'vbf2', args.configB).testSelection("Best")[0]
bkg_vbf3 = profileClass(x, mu_gauss, 'vbf3', args.configB).testSelection("Best")[0]
bkg_vbf4 = profileClass(x, mu_gauss, 'vbf4', args.configB).testSelection("Best")[0]

hist_content = np.zeros(int(4*(highx-lowx)))
hist_errorsq = np.zeros(int(4*(highx-lowx)))

nSigFit = [12.992, 78.563, 205.64, 250.77, 4.5228, 9.6859, -7.0575, 113.58]
norm_fit = np.zeros(8)
ntot_weighted = 0
x.setBins(int(4*(highx-lowx)))
combine_hist = ROOT.RooDataHist('hist', 'hist', x)
sig_hist = ROOT.RooDataHist('sig_hist', 'sig_hist', x)
CATS = ['ggf1', 'ggf2', 'ggf3', 'ggf4', 'vbf1', 'vbf2', 'vbf3', 'vbf4']
for j, CAT in enumerate(CATS):
    f = ROOT.TFile('../Unblinding_toy/higgsCombine.'+CAT+'.MultiDimFit.mH125.root')
    w = f.Get("w")
    x_ = w.var('CMS_hzg_mass_'+CAT)
    lowx_ = x_.getMin()
    x_.setBins(260)
    data = w.data("data_obs")
    #norm_plot = data.sumEntries('CMS_hzg_mass_'+CAT+'> '+str(lowx)+ '&& CMS_hzg_mass_'+CAT+' < '+str(highx))
    norm_fit[j] = data.sumEntries()
    hdata = ROOT.RooDataHist('hdata', 'hdata', x_, data)
    norm_plot = 0
    for i in range(int(4*(highx-lowx))):
        hdata.get(int(i+4*(lowx - lowx_)))
        norm_plot += hdata.weight()
    ratio = nSigExp[j]/norm_plot
    ntot_weighted += nSigExp[j]*ratio
    print('ratio = ', ratio)
    for i in range(4*(highx-lowx)):
        hdata.get(int(i+4*(lowx - lowx_)))
        hist_content[i] += hdata.weight()*ratio
        hist_errorsq[i] += hdata.weight()*ratio*ratio
for i in range(int(4*(highx-lowx))):
    combine_hist.get(i)
    combine_hist.set(hist_content[i], math.sqrt(hist_errorsq[i]))
    #print (math.sqrt(hist_errorsq[i]))    
c2_ggf1 = ROOT.RooRealVar('c2_ggf1', 'c2_ggf1', nSigFit[0])
c2_ggf2 = ROOT.RooRealVar('c2_ggf2', 'c2_ggf2', nSigFit[1])
c2_ggf3 = ROOT.RooRealVar('c2_ggf3', 'c2_ggf3', nSigFit[2])
c2_ggf4 = ROOT.RooRealVar('c2_ggf4', 'c2_ggf4', nSigFit[3])
c2_vbf1 = ROOT.RooRealVar('c2_vbf1', 'c2_vbf1', nSigFit[4])
c2_vbf2 = ROOT.RooRealVar('c2_vbf2', 'c2_vbf2', nSigFit[5])
c2_vbf3 = ROOT.RooRealVar('c2_vbf3', 'c2_vbf3', nSigFit[6])
c2_vbf4 = ROOT.RooRealVar('c2_vbf4', 'c2_vbf4', nSigFit[7])
c1_ggf1 = ROOT.RooRealVar('c1_ggf1', 'c1_ggf1', norm_fit[0])
c1_ggf2 = ROOT.RooRealVar('c1_ggf2', 'c1_ggf2', norm_fit[1])
c1_ggf3 = ROOT.RooRealVar('c1_ggf3', 'c1_ggf3', norm_fit[2])
c1_ggf4 = ROOT.RooRealVar('c1_ggf4', 'c1_ggf4', norm_fit[3])
c1_vbf1 = ROOT.RooRealVar('c1_vbf1', 'c1_vbf1', norm_fit[4])
c1_vbf2 = ROOT.RooRealVar('c1_vbf2', 'c1_vbf2', norm_fit[5])
c1_vbf3 = ROOT.RooRealVar('c1_vbf3', 'c1_vbf3', norm_fit[6])
c1_vbf4 = ROOT.RooRealVar('c1_vbf4', 'c1_vbf4', norm_fit[7])

model_ggf1 = ROOT.RooAddPdf('model_ggf1', 'model_ggf1', ROOT.RooArgList(combine_model_ggf1.pdf, bkg_ggf1.pdf), ROOT.RooArgList(c2_ggf1, c1_ggf1))
model_ggf2 = ROOT.RooAddPdf('model_ggf2', 'model_ggf2', ROOT.RooArgList(combine_model_ggf2.pdf, bkg_ggf2.pdf), ROOT.RooArgList(c2_ggf2, c1_ggf2))
model_ggf3 = ROOT.RooAddPdf('model_ggf3', 'model_ggf3',  ROOT.RooArgList(combine_model_ggf3.pdf, bkg_ggf3.pdf), ROOT.RooArgList(c2_ggf3, c1_ggf3))
model_ggf4 = ROOT.RooAddPdf('model_ggf4', 'model_ggf4',  ROOT.RooArgList(combine_model_ggf4.pdf, bkg_ggf4.pdf), ROOT.RooArgList(c2_ggf4, c1_ggf4))
model_vbf1 = ROOT.RooAddPdf('model_vbf1', 'model_vbf1',  ROOT.RooArgList(combine_model_vbf1.pdf, bkg_vbf1.pdf), ROOT.RooArgList(c2_vbf1, c1_vbf1))
model_vbf2 = ROOT.RooAddPdf('model_vbf2', 'model_vbf2',  ROOT.RooArgList(combine_model_vbf2.pdf, bkg_vbf2.pdf), ROOT.RooArgList(c2_vbf2, c1_vbf2))
model_vbf3 = ROOT.RooAddPdf('model_vbf3', 'model_vbf3',  ROOT.RooArgList(combine_model_vbf3.pdf, bkg_vbf3.pdf), ROOT.RooArgList(c2_vbf3, c1_vbf3))
model_vbf4 = ROOT.RooAddPdf('model_vbf4', 'model_vbf4',  ROOT.RooArgList(combine_model_vbf4.pdf, bkg_vbf4.pdf), ROOT.RooArgList(c2_vbf4, c1_vbf4))

tot_model = ROOT.RooAddPdf('S+B', 'S+B', ROOT.RooArgList(model_ggf1, model_ggf2, model_ggf3, model_ggf4, model_vbf1, model_vbf2, model_vbf3, model_vbf4), ROOT.RooArgList(wc1, wc2, wc3, wc4, wc5, wc6, wc7))
b_only = ROOT.RooAddPdf('B-Only', 'B-Only', ROOT.RooArgList(bkg_ggf1.pdf, bkg_ggf2.pdf, bkg_ggf3.pdf, bkg_ggf4.pdf, bkg_vbf1.pdf, bkg_vbf2.pdf, bkg_vbf3.pdf, bkg_vbf4.pdf), ROOT.RooArgList(wc1, wc2, wc3, wc4, wc5, wc6, wc7))
s_only = ROOT.RooAddPdf('S-Only', 'S-Only', ROOT.RooArgList(combine_model_ggf1.pdf, combine_model_ggf2.pdf, combine_model_ggf3.pdf, combine_model_ggf4.pdf, combine_model_vbf1.pdf, combine_model_vbf2.pdf, combine_model_vbf3.pdf, combine_model_vbf4.pdf), ROOT.RooArgList(wc1, wc2, wc3, wc4, wc5, wc6, wc7))
b_hist = b_only.generateBinned(x, combine_hist.sumEntries() - ntot_weighted, True)
for i in range(int(4*(highx-lowx))):
    sig_hist.get(i)
    b_hist.get(i)
    sig_hist.set(hist_content[i] - b_hist.weight(), math.sqrt(hist_errorsq[i]))
#plotClass(x, combine_hist, tot_model, tot_model, title='Money_plot', fullRatio = True)
#plot = x.frame()
#combine_hist.plotOn(plot, ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))
#can = ROOT.TCanvas('can','can', 600, 500)
#plot.Draw()
#can.SaveAs('test.pdf')

# Plotting starts from here
ROOT.gStyle.SetOptStat(0)
CMS = "Internal"
CMS_lumi.lumi_sqrtS = "137.61 fb^{-1} (13 TeV) + 62.32 fb^{-1} (13.6 TeV)"
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "      " + CMS
CMS_lumi.cmsTextSize = 0.5
CMS_lumi.lumiTextSize = 0.3

print (ntot_weighted)
x.setBins(int(highx-lowx))
plot = x.frame()
plot.SetTitle("")
plot2 = x.frame()
plot2.SetTitle("")
show_hist = ROOT.RooDataHist("show_hist", "show_hist", ROOT.RooArgSet(x), combine_hist)
show_sig_hist = ROOT.RooDataHist("show_sig_hist", "show_sig_hist", ROOT.RooArgSet(x), sig_hist)
show_hist.plotOn( plot, ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2), ROOT.RooFit.Name('hist'))
b_only.plotOn(plot,  ROOT.RooFit.LineColor(6), ROOT.RooFit.LineWidth(2), ROOT.RooFit.Name('bonly_model'), ROOT.RooFit.Normalization((combine_hist.sumEntries() - ntot_weighted)/combine_hist.sumEntries()), ROOT.RooFit.LineStyle(2))
tot_model.plotOn(plot,  ROOT.RooFit.LineColor(2), ROOT.RooFit.LineWidth(2), ROOT.RooFit.Name('model'))
show_sig_hist.plotOn( plot2, ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))
print (combine_hist.sumEntries())
print (show_sig_hist.sumEntries())
s_only.plotOn(plot2,  ROOT.RooFit.LineColor(2), ROOT.RooFit.LineWidth(2))#, ROOT.RooFit.Normalization(ntot_weighted, ROOT.RooAbsReal.NumEvent))

can = ROOT.TCanvas("c1","c1", 700, 500)
can.Divide(1,2)
can.cd(1)
ROOT.gPad.SetPad(0.0, 0.4, 1.0, 0.98)
ROOT.gPad.SetBottomMargin(0)
ROOT.gPad.SetTopMargin(0.12)
plot.GetYaxis().SetTitleSize(0.05)
plot.GetYaxis().SetTitleOffset(0.5)
plot.GetYaxis().SetLabelSize(0.05)
plot.GetYaxis().SetTitle('S/(S+B) Weighted Events/(1)')
plot.Draw()
leg = ROOT.TLegend(.7,.5,.9,.8)
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetFillStyle(0)
#leg.SetTextFont(50)
leg.SetTextSize(0.05)
leg.AddEntry('hist','Toy',"LP")
leg.AddEntry('model',tot_model.GetName(),"L")
leg.AddEntry('bonly_model',b_only.GetName(),"L")
leg.Draw("SAME")
latex = ROOT.TLatex()
latex.SetTextSize(0.05)
latex.SetTextAlign(13)
#latex.SetTextFont(50)
latex.DrawLatexNDC(.65, .45, '#mu = 1.557 + 0.547(-0.606)')
CMS_lumi.CMS_lumi(can, 0, 0)
can.cd(2)
ROOT.gPad.SetPad(0.0, 0.1, 1.0, 0.4)
ROOT.gPad.SetTopMargin(0)
ROOT.gPad.SetBottomMargin(0.22)
plot2.GetYaxis().SetTitleSize(0.1)
plot2.GetYaxis().SetTitle('Data - B')
plot2.GetYaxis().SetTitleOffset(0.25)
plot2.GetYaxis().SetLabelSize(0.1)
plot2.GetXaxis().SetLabelSize(0.1)
plot2.GetXaxis().SetTitleSize(0.1)
plot2.GetXaxis().SetTitle('m_{#font[12]{ll}\gamma}(GeV)')
plot2.Draw()
line = ROOT.TLine( x.getMin(), 0, x.getMax(), 0)
line.SetLineColor(ROOT.kBlack)
line.SetLineStyle(7)
line.SetLineWidth(2)
line.Draw("same")
can.SaveAs('plots/Money_plot.pdf')
