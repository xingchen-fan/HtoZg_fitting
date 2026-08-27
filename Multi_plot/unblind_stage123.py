#!/usr/bin/env python3
import sys
import os
import ROOT
import argparse
import json
sys.path.append(os.path.abspath("../Utilities/"))
sys.path.append(os.path.abspath("../CMS_plotter/"))
from Xc_Minimizer import *
ROOT.gStyle.SetOptStat(0)
import numpy as np
import pandas as pd
import glob
import math
import ctypes
import CMS_lumi, tdrstyle
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)

def set_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--cat', help = 'Category')
    parser.add_argument('-con', '--config', help = 'Configuration', default='')
    parser.add_argument('-f', '--func', help = 'Function postfit')
    parser.add_argument('-i', '--input', help = 'Input Root File')
    parser.add_argument('-t', '--toy', help = 'Toy file?', default='')
    parser.add_argument('-l', '--lep', help = 'Lepton', default='')
    parser.add_argument('-s', '--stage', help = 'Stage', default='1')
    parser.add_argument('-eb', '--errBand', help = 'Error band?', default=False)
    parser.add_argument('-ed', '--errDir', help = 'Error band toy dir', default='')
    parser.add_argument('-cl', '--catLabel', help = 'Cat label', default='')
    return parser.parse_args()

def extractBandProperties(data,category,bidx):
  props = {}
  if category == 'all': c = 'sum'
  elif category == 'wall': c = 'wsum'
  else: c = category
  props['median'] = np.median(data['%s_%g'%(c,bidx)].values)
  props['up1sigma'] = np.percentile(data['%s_%g'%(c,bidx)].values,50*(1+math.erf(1./math.sqrt(2))))
  props['down1sigma'] = np.percentile(data['%s_%g'%(c,bidx)].values,50*(1+math.erf(-1./math.sqrt(2))))
  props['up2sigma'] = np.percentile(data['%s_%g'%(c,bidx)].values,50*(1+math.erf(2./math.sqrt(2))))
  props['down2sigma'] = np.percentile(data['%s_%g'%(c,bidx)].values,50*(1+math.erf(-2./math.sqrt(2))))
  return props

def findPrefitNorm(x_, pdf_, nSB_):
    x_.setRange("full", x_.getMin(), x_.getMax())
    x_.setRange("signal", 120, 130)
    total = pdf_.createIntegral(ROOT.RooArgSet(x_), ROOT.RooFit.Range("full"))
    signal_window = pdf_.createIntegral(ROOT.RooArgSet(x_), ROOT.RooFit.Range("signal"))
    return nSB_ / (1 - signal_window.getVal()/total.getVal())
    
args = set_args()
CAT = args.cat
jfile = open(args.config, 'r')
configs = json.load(jfile)
setting = configs[CAT]

cats = ['ggf4', 'ggf3', 'ggf2', 'ggf1', 'vbf4', 'vbf3', 'vbf2', 'vbf1', 'vh3l', 'vhptmiss', 'tthhad', 'tthlep', 'untag']
cat_index = cats.index(CAT)

f = ROOT.TFile(args.input)
w = f.Get("w")

if args.toy != '':
    data = w_h.data("hist_toy5_"+CAT)
else:
    data = w.data("data_obs").reduce(f'CMS_channel=={cat_index}')
lowx = setting["Range"][0]
highx = setting["Range"][1]
nbins = int(setting["Bins"])
binning = ROOT.RooFit.Binning(65,lowx,highx)

x = w.var('CMS_hzg_mass_'+CAT)
x.setRange('left', lowx, 120)
x.setRange('right', 130, highx)
x.setBins(65)

plot = x.frame()
plot.SetTitle("")#(CAT + ' Unblind Stage ' +args.stage + ' '+args.lep)
ROOT.gStyle.SetTitleSize(0.08, "t")
plot.GetYaxis().SetTitleSize(0.05)
plot.GetYaxis().SetTitleOffset(0.6)
plot.GetYaxis().SetLabelSize(0.05)
plot.GetYaxis().SetTitle("Events/GeV")


print(type(data))
ntot = data.sumEntries()
#dataSB = data.reduce(ROOT.RooFit.CutRange('left,right'))
dataSB = data.reduce(f"CMS_hzg_mass_{CAT} < 120 || CMS_hzg_mass_{CAT} > 130")

nSB = dataSB.sumEntries()
nbkg_prefit = findPrefitNorm(x, w.pdf('model_b').getPdf("cat_"+CAT), nSB) #w.function('n_exp_final_bincat_'+CAT+'_proc_background').getVal()
print("Prefit norm = ",nbkg_prefit, ", Expect bkg in signal window = ", nbkg_prefit - nSB)
print("n tot SB= ", nSB, "ntot = ", ntot)
norm_ = 0
if args.stage == '1':
    norm_ = nSB
    show_hist = dataSB
elif args.stage == '2':
    norm_ = ntot
    show_hist = data
elif args.stage == '3':
    norm_ = ntot
    show_hist = data
else:
    raise Exception("Wrong stage.")
show_hist.plotOn( plot, binning, ROOT.RooFit.Name('hist'), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))



# Error band
cats = ['ggf4', 'ggf3', 'ggf2', 'ggf1', 'vbf4', 'vbf3', 'vbf2', 'vbf1', 'vh3l', 'vhptmiss', 'tthhad', 'tthlep', 'untag']
cat_index = cats.index(CAT)
cols = []
nbad = 0
for i in range(1, 66):
    cols.append(f"{CAT}_{i}")
df_bands = pd.DataFrame(columns=cols)
if args.errBand:
    toyFiles = glob.glob("%s/toy_*.root"%args.errDir)
    for it in range(len(toyFiles)):
        ftoy = ROOT.TFile(toyFiles[it])
        toy = ftoy.Get("toys/toy_asimov")
        vetoToy = False
        values = {}
        dtoy = toy.reduce(f'CMS_channel=={cat_index}')
        htoy = x.createHistogram(f"h_{CAT}",ROOT.RooFit.Binning(65,x.getMin(),x.getMax()))
        dtoy.fillHistogram(htoy,ROOT.RooArgList(x))
        #htoy.Rebin(4)
        for ibin in range(1,htoy.GetNbinsX()+1):
            v = htoy.GetBinContent(ibin)
            if v!=v: vetoToy = True # Veto toys which have a NaN
            values['%s_%g'%(CAT,ibin)] = v
        if values['%s_1'%CAT] == 0: vetoToy = True # Veto toys with 1st bin is empty
        htoy.Delete()
        dtoy.Delete()
        toy.Delete()
        ftoy.Close()
        if not vetoToy: df_bands.loc[len(df_bands)] = values
        if vetoToy: nbad  += 1
print('nbad = ', nbad)
b_hist = ROOT.TH1F()
sb_hist = ROOT.TH1F()

CMS = "" #"Supplementary"
CMS_lumi.lumi_sqrtS = "138 fb^{-1} (13 TeV) + 62 fb^{-1} (13.6 TeV)"
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "        " + CMS
CMS_lumi.cmsTextSize = 0.7
CMS_lumi.lumiTextSize = 0.4
latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextSize(0.08)
latex.SetTextFont(42)

can = ROOT.TCanvas("c","c", 700, 500)
can.Divide(1,2)
leg = ROOT.TLegend(0.7 if args.stage == '3' else 0.55,0.5,0.85,0.8)
leg.SetLineColor(0)
leg.SetTextSize(0.06)
if args.stage == '1' or args.stage == '2':
# Prefit
    b_model_1 = w.pdf('model_b').getPdf("cat_"+CAT)#w.pdf(args.func1+'_'+CAT+'_model')
    b_model_1.plotOn( plot, ROOT.RooFit.LineColor(2), ROOT.RooFit.Name("prefit"), ROOT.RooFit.Normalization(nbkg_prefit/norm_), ROOT.RooFit.LineWidth(2))

# Postfit
    b_model_2 = w.pdf('model_b').getPdf("cat_"+CAT)
    w.loadSnapshot("MultiDimFit")
    nbkg = w.function('n_exp_final_bincat_'+CAT+'_proc_background').getVal()
    print("post fit norm = ", nbkg)
    x.setBins(nbins)
    hdata = ROOT.RooDataHist('hdata', 'hdata', x, data)
    chi2 = ROOT.RooChi2Var('chi2', 'chi2', b_model_2, hdata)
    val_chi2 = chi2.getVal()
    x.setBins(int(highx - lowx))
    latex = ROOT.TLatex()
    latex.SetTextSize(0.03)
    latex.SetTextAlign(13)
    pars = w.pdf(args.func+'_'+CAT+'_model').getParameters(x)
    n_float = 0
    for p in pars:
        if isinstance(p, ROOT.RooRealVar) and not p.isConstant(): n_float += 1
    print("nfloat = ", n_float)
    dof = nbins - n_float
    PV = ROOT.TMath.Prob(val_chi2, dof)
    b_hist = b_model_2.generateBinned(x,  nbkg, True).createHistogram("b_hist", x, ROOT.RooFit.Binning(65))
    b_model_2.plotOn( plot, ROOT.RooFit.LineColor(4), ROOT.RooFit.Name("postfit"), ROOT.RooFit.Normalization(nbkg/norm_), ROOT.RooFit.LineWidth(2))
    latex.DrawLatexNDC(.7, .65, '#splitline{#Chi^{2}/dof = %i'%val_chi2+ '/%i'%dof+'}{P-value = %.2f}'%PV)
else:
    b_model_2 = w.pdf('model_b').getPdf("cat_"+CAT)
    sb_model = w.pdf('model_s').getPdf("cat_"+CAT)
    w.loadSnapshot("MultiDimFit")
    nbkg = w.function('n_exp_final_bincat_'+CAT+'_proc_background').getVal()
    x.setBins(nbins)
    b_hist = b_model_2.generateBinned(x, nbkg, True).createHistogram("b_hist", x, ROOT.RooFit.Binning(65))
    b_model_2.plotOn( plot, ROOT.RooFit.LineColor(2), ROOT.RooFit.Name("b_only"), ROOT.RooFit.LineWidth(2), ROOT.RooFit.LineStyle(2), ROOT.RooFit.Normalization(nbkg/norm_))
    sb_model.plotOn( plot, ROOT.RooFit.LineColor(2), ROOT.RooFit.Name("s+b"),ROOT.RooFit.LineWidth(2))
    sb_hist = sb_model.generateBinned(x, ntot, True).createHistogram("sb_hist", x, ROOT.RooFit.Binning(65))
    sb_hist.Add(b_hist, -1)
    sb_hist.SetLineWidth(2)
can.cd(1)
ROOT.gPad.SetPad(0.0, 0.35, 1.0, 0.99)
ROOT.gPad.SetBottomMargin(0)
ROOT.gPad.SetTopMargin(0.15)
plot.GetYaxis().SetTitleOffset(0.75)
plot.GetYaxis().SetTitleSize(0.065)
plot.Draw()
ROOT.gPad.Update()   # Make sure the axis exists
plot.GetYaxis().ChangeLabel(1, -1, 0, -1, -1, -1, "")
latex.DrawLatex(0.15, 0.7, args.catLabel)

h_hist = show_hist.createHistogram("h_hist", x, ROOT.RooFit.Binning(65))

# Error band
gr_1sig_r, gr_2sig_r = ROOT.TGraphAsymmErrors(), ROOT.TGraphAsymmErrors()
gr_i = 0
if args.errBand:
    for ibin in range(1, 66):
        xval = h_hist.GetXaxis().GetBinCenter(ibin)
        xerr = 0.5*(h_hist.GetXaxis().GetBinWidth(ibin))
        properties = extractBandProperties(df_bands,CAT,ibin)
        bkgval = b_hist.GetBinContent(ibin)
        gr_1sig_r.SetPoint(gr_i,xval,properties['median']-bkgval)
        gr_2sig_r.SetPoint(gr_i,xval,properties['median']-bkgval)
        gr_1sig_r.SetPointError(gr_i,xerr,xerr,properties['median']-properties['down1sigma'],properties['up1sigma']-properties['median'])
        gr_2sig_r.SetPointError(gr_i,xerr,xerr,properties['median']-properties['down2sigma'],properties['up2sigma']-properties['median'])
        gr_i += 1
        
gr_1sig_r.SetFillColorAlpha(ROOT.kGreen+1, 0.4)
gr_2sig_r.SetFillColorAlpha(ROOT.kOrange, 0.4)
gr_1sig_r.SetLineWidth(0)
gr_2sig_r.SetLineWidth(0)

leg.AddEntry("hist", "Data", "LP")
if args.stage == '3':
    leg.AddEntry("b_only", "B model", "L")
    leg.AddEntry("s+b", "S+B model", "L")
else:
    leg.AddEntry("prefit", "Prefit B-only model", "L")
    leg.AddEntry("postfit", "Postfit B-only model", "L")
if args.errBand:
    leg.AddEntry(gr_1sig_r,"95% CI","F")
    leg.AddEntry(gr_2sig_r,"68% CI","F")
leg.Draw("same")
ROOT.gPad.Update()
if args.stage == '1':
    box = ROOT.TBox(120,0, 130, ROOT.gPad.GetFrame().GetY2())
    box.SetFillColorAlpha(16, 0.4)
    box.Draw("same")
    
CMS_lumi.CMS_lumi(can, 0, 0)
                  
can.cd(2)
ROOT.gPad.SetPad(0.0, 0.1, 1.0, 0.35)
ROOT.gPad.SetTopMargin(0)
ROOT.gPad.SetBottomMargin(0.42) #38
#h_hist.Sumw2(False)

for i in range(0, b_hist.GetNbinsX() + 2):
    b_hist.SetBinError(i, 0.0)
h_hist.SetBinErrorOption(ROOT.TH1.kPoisson)
"""
h_hist.Add(b_hist, -1)
h_hist.SetMarkerColor(ROOT.kBlack)
h_hist.SetMarkerStyle(8)
h_hist.SetMarkerSize(1)
h_hist.SetLineColor(ROOT.kBlack)
"""
h_ratio = ROOT.TGraphAsymmErrors()
h_ratio.SetMarkerColor(ROOT.kBlack)
h_ratio.SetMarkerStyle(8)
h_ratio.SetMarkerSize(1)
h_ratio.SetLineColor(ROOT.kBlack)
rooHist = plot.getObject(0)
for i in range(rooHist.GetN()):
    x_ = ctypes.c_double()
    y_ = ctypes.c_double()
    rooHist.GetPoint(i, x_, y_)

    y_res = y_.value - b_hist.GetBinContent(i+1)
    h_ratio.SetPoint(i, float(x_.value), y_res)

    h_ratio.SetPointError(
        i,
        rooHist.GetErrorXlow(i),
        rooHist.GetErrorXhigh(i),
        rooHist.GetErrorYlow(i),
        rooHist.GetErrorYhigh(i)
    )

"""    
for i in range(1, h_hist.GetNbinsX() + 1):
    x_ = h_hist.GetBinCenter(i)
    y_ = h_hist.GetBinContent(i) - b_hist.GetBinContent(i)

    ex = h_hist.GetBinWidth(i) / 2.0
    eyl = h_hist.GetBinErrorLow(i)
    eyh = h_hist.GetBinErrorUp(i)
    print(eyl, eyh)
    point = i - 1
    h_ratio.SetPoint(point, x_, y_)
    h_ratio.SetPointError(point, ex, ex, eyl, eyh)
"""
h_range = x.createHistogram(f"h_range_{CAT}",ROOT.RooFit.Binning(65,x.getMin(),x.getMax()))
h_range.SetTitle("")
h_range.GetXaxis().SetTitleOffset(1)
h_range.GetXaxis().SetTitleSize(0.2)#18
h_range.GetXaxis().SetTitle("m_{#font[12]{ll}\gamma}(GeV)")
h_range.GetXaxis().SetLabelSize(0.15)#13
h_range.GetYaxis().SetLabelSize(0.10)
h_range.GetYaxis().SetTitleOffset(0.3)
h_range.GetYaxis().SetTitleSize(0.15)
h_range.GetYaxis().SetNdivisions(4)
h_range.GetYaxis().SetTitle("Data - B model")
minimum = 99999
h_hist.Add(b_hist, -1)
for i in range(int(x.getMax() - x.getMin())):
    bincont = h_hist.GetBinContent(i+1)
    if bincont != 0 and bincont < minimum and (x.getMin() + i < 120 or x.getMin() + i > 130) :
        minimum = bincont
range_ = 1.1 * max(abs(minimum), abs(h_hist.GetBinContent(h_hist.GetMaximumBin())))

line = ROOT.TLine( x.getMin(), 0, x.getMax(), 0)
line.SetLineColor(ROOT.kBlack)
line.SetLineStyle(7)
line.SetLineWidth(2)
h_range.GetYaxis().SetRangeUser(-range_, range_)
h_range.Draw()
line.Draw("SAME")
if args.errBand:
    gr_2sig_r.Draw("LE3SAME")
    gr_1sig_r.Draw("LE3SAME")
h_ratio.Draw("P SAME")
if args.stage == '3':
    sb_hist.SetLineColor(2)
    sb_hist.Draw("HIST C SAME")
ROOT.gPad.Update()
if args.stage == '1':
    box2 = ROOT.TBox(120,ROOT.gPad.GetFrame().GetY1(), 130, ROOT.gPad.GetFrame().GetY2())
    box2.SetFillColorAlpha(16, 0.4)
    box2.Draw("same")

can.Update()
can.SaveAs('plots/unblind_stage'+args.stage+'_'+CAT+args.lep+'.pdf')
#can.SaveAs('plots/unblind_stage'+args.stage+'_'+CAT+args.lep+'.root')
