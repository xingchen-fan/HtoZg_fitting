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
import pandas as pd
import glob
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)
ROOT.gInterpreter.AddIncludePath('../Utilities/RooGaussStepBernstein.h')
ROOT.gSystem.Load('../Utilities/RooGaussStepBernstein_cxx.so')

def findCatIndex(cat):
    cats_ = ['ggf4', 'ggf3', 'ggf2', 'ggf1', 'vbf4', 'vbf3', 'vbf2', 'vbf1', 'vh3l', 'vhptmiss', 'tthhad', 'tthlep', 'untag']
    return cats_.index(cat)

def findYRange(x, hist):
    maxi = -1
    for i in range(int(x.getMax() - x.getMin())):
        hist.get(i)
        bincont = hist.weight()
        if abs(bincont) > maxi:
            maxi = abs(bincont)
    return 1.1 * maxi

def set_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-conS', '--configS', help = 'Signal Configuration')
    parser.add_argument('-i', '--input', help = 'Input of postfit combine output file')
    parser.add_argument('-c', '--cat', help = 'Category', default='all')
    parser.add_argument('-t', '--errToy', help = 'Error band toy dir')
    return parser.parse_args()

def readErrorBandToys(cats, weights, lowx, highx, range_l, range_h):
    cols = []
    toyFiles = glob.glob("%s/toy_*.root"%args.errToy)
    for i in range(1, int(highx - lowx + 1)):
        cols.append(f"wsum_{i}")
    df_bands = pd.DataFrame(columns=cols)
    nbad = 0
    for it in range(len(toyFiles)):
        ftoy = ROOT.TFile(toyFiles[it])
        toy = ftoy.Get("toys/toy_asimov")
        values = {}
        for ibin in range(1,int(highx - lowx + 1)):
            values['wsum_%g'%ibin] = 0
        vetoToy = False
        for cat in cats:
            c_ind = findCatIndex(cat)
            #if cat == 'untag': print('ID = ', c_ind)
            dtoy = toy.reduce(f'CMS_channel=={c_ind}')
            x = ROOT.RooRealVar('CMS_hzg_mass_'+cat, 'CMS_hzg_mass_'+cat, int(range_l[cat]), int(range_h[cat]))
            htoy = dtoy.createHistogram(f"h_{cat}", x, ROOT.RooFit.Binning(65,x.getMin(),x.getMax()))
            for ibin in range(1,int(highx - lowx + 1)):
                ind_sel = int(ibin + lowx - range_l[cat])
                v = htoy.GetBinContent(ind_sel)
                if v!=v: vetoToy = True # Veto toys which have a NaN
                values['wsum_%g'%ibin] += weights[cat]*v
                #if v == 0 and ibin == 1: vetoToy = True # Veto toys with 1st bin is empty
            htoy.Delete()
            dtoy.Delete()
        toy.Delete()
        ftoy.Close()
        if not vetoToy: df_bands.loc[len(df_bands)] = values
        if vetoToy: nbad  += 1
    print('n bad toys = ', nbad)
    return df_bands

def extractBandProperties(data,bidx):
  props = {}
  c = 'wsum'
  props['median'] = np.median(data['%s_%g'%(c,bidx)].values)
  props['up1sigma'] = np.percentile(data['%s_%g'%(c,bidx)].values,50*(1+math.erf(1./math.sqrt(2))))
  props['down1sigma'] = np.percentile(data['%s_%g'%(c,bidx)].values,50*(1+math.erf(-1./math.sqrt(2))))
  props['up2sigma'] = np.percentile(data['%s_%g'%(c,bidx)].values,50*(1+math.erf(2./math.sqrt(2))))
  props['down2sigma'] = np.percentile(data['%s_%g'%(c,bidx)].values,50*(1+math.erf(-2./math.sqrt(2))))
  return props

args = set_args()
lowx = 106
highx = 159
nbins = int(highx - lowx)
binning = ROOT.RooFit.Binning(nbins,lowx,highx)

jfileS = open(args.configS, 'r')
config = json.load(jfileS)

x = ROOT.RooRealVar('x', 'x', lowx, highx)
x_c = ROOT.RooRealVar('x_c', 'x_c', lowx, highx)

DEBUG = False
hist_content = np.zeros(int(highx-lowx))
hist_errorsq = np.zeros(int(highx-lowx))
hist_content_B = np.zeros(int(4 *(highx-lowx)))
hist_content_SB = np.zeros(int(4 * (highx-lowx)))

ntot_weighted = 0
x.setBins(int(highx-lowx))
x_c.setBins(int(4*(highx-lowx)))

combine_hist = ROOT.RooDataHist('hist', 'hist', x) # Weighted data
sig_hist = ROOT.RooDataHist('sig_hist', 'sig_hist', x) # Data - B only
combine_B_hist = ROOT.RooDataHist('hist_B', 'hist_B', x_c) # B only model
combine_SB_hist = ROOT.RooDataHist('hist_SB', 'hist_SB', x_c) # S+B model
combine_S_hist = ROOT.RooDataHist('hist_S', 'hist_S', x_c) # S model

if args.cat == 'all':
    CATS = ['ggf1', 'ggf2', 'ggf3', 'ggf4', 'vbf1', 'vbf2', 'vbf3', 'vbf4', 'vh3l', 'vhptmiss', 'untag', 'tthlep', 'tthhad']
    #funcs = ['lau2', 'bern3', 'pow2', 'lau3', 'exp1', 'lau2', 'bern3', 'bern2', 'bern2', 'pow1', 'exp2', 'bern2', 'bern2']
elif args.cat == 'ggf':
    CATS = ['ggf1', 'ggf2', 'ggf3', 'ggf4']
elif args.cat == 'vbf':
    CATS = ['vbf1', 'vbf2', 'vbf3', 'vbf4']
elif args.cat == 'vhtth':
    CATS = ['vh3l', 'vhptmiss', 'untag', 'tthlep', 'tthhad']
else:
    CATS = [args.cat]
    
gr_1sig_r, gr_2sig_r = ROOT.TGraphAsymmErrors(), ROOT.TGraphAsymmErrors()
gr_i = 0
weights = {}
range_l = {}
range_h = {}
for j, CAT in enumerate(CATS):
    cat_index = findCatIndex(CAT)
    
    # Read combine output file
    f_ = ROOT.TFile(args.input)
    w_ = f_.Get("w")
    x_ = w_.var('CMS_hzg_mass_'+CAT)
    b_model = w_.pdf('model_b').getPdf("cat_"+CAT)
    sb_model = w_.pdf('model_s').getPdf("cat_"+CAT)
    w_.loadSnapshot("MultiDimFit")
    
    lowx_ = x_.getMin()
    range_l[CAT] = lowx_
    range_h[CAT] = x_.getMax()
    x_.setBins(260)

    data = w_.data("data_obs").reduce(f'CMS_channel=={cat_index}')
    data_hist = ROOT.RooDataHist(f'h_{CAT}_data', f'h_{CAT}_data', x_, data)
    ntot = data_hist.sumEntries()
    print(f"n {CAT} = ", ntot)

    # Generate b-only model hist and calculate exp S/(S+B) ratio
    nbkg = w_.function(f'n_exp_final_bincat_{CAT}_proc_background').getVal()
    b_hist = b_model.generateBinned(x_,  nbkg, True)
    nbkg_sig = b_hist.reduce(f"CMS_hzg_mass_{CAT} > 120 && CMS_hzg_mass_{CAT} < 130").sumEntries()
    nsig = config[f'Htozg_el_cat_{CAT}_nominal']['nexp'] + config[f'Htozg_mu_cat_{CAT}_nominal']['nexp']
    ratio = nsig/(nsig+nbkg_sig)
    weights[CAT] = ratio

    if DEBUG:
        print('ratio = ', ratio)
        print('Post-fit nbkg = ', nbkg)

    # Generate s+b model hist
    sb_hist = sb_model.generateBinned(x_, ntot, True)
    #sb_hist = sb_Rhist.createHistogram(f"{CAT}_sb_hist", x_, ROOT.RooFit.Binning(65))
    if DEBUG:
        print('N b_hist  = ', b_hist.sumEntries())
        print('N sb_hist  = ', sb_hist.sumEntries())
        
    # Store weighted data, b model and s+b model hist to array
    for i in range(int(4 * (highx-lowx))):
        #data_hist.get(int(i+4*(lowx - lowx_))))
        b_hist.get(int(i+4*(lowx - lowx_)))
        sb_hist.get(int(i+4*(lowx - lowx_)))
        #hist_content[i] += data_hist.weight()*ratio
        #hist_errorsq[i] += data_hist.weight()*ratio*ratio
        hist_content_B[i] += b_hist.weight()*ratio
        hist_content_SB[i] += sb_hist.weight()*ratio
    for i in range(int(highx-lowx)):
        bin_w = 0
        for j in range(4):
            data_hist.get(int(j + 4*(lowx - lowx_ +i)))
            bin_w += data_hist.weight()
        hist_content[i] += bin_w*ratio
        hist_errorsq[i] += bin_w*ratio*ratio
        
# Fill the histograms
for i in range(int(4 * (highx-lowx))):
    combine_B_hist.get(i)
    combine_SB_hist.get(i)
    combine_S_hist.get(i)
    combine_B_hist.set(4 * hist_content_B[i], 0)
    combine_SB_hist.set(4 * hist_content_SB[i], 0)
    combine_S_hist.set(4 * (hist_content_SB[i] - hist_content_B[i]), 0)

for i in range(int(highx-lowx)):
    combine_hist.get(i)
    sig_hist.get(i)
    combine_hist.set(hist_content[i], math.sqrt(hist_errorsq[i]))
    bin_w = 0
    for j in range(4):
        bin_w += hist_content_B[int(4*i + j)]
    sig_hist.set(hist_content[i] - bin_w, math.sqrt(hist_errorsq[i]))

if DEBUG:
    print ('S+B hist = ',combine_SB_hist.sumEntries())
    print ('B hist = ', combine_B_hist.sumEntries())
    print ('Combined hist = ',combine_hist.sumEntries())
    print ('Combined - B = ', sig_hist.sumEntries())
    
print(weights)

# Error band
h_hist = combine_hist.createHistogram("h_hist", x)
df_bands = readErrorBandToys(CATS, weights, lowx, highx, range_l, range_h)
for ibin in range(highx - lowx):
    xval = h_hist.GetXaxis().GetBinCenter(ibin+1)
    xerr = 0.5*(h_hist.GetXaxis().GetBinWidth(ibin+1))
    properties = extractBandProperties(df_bands,ibin+1)
    bkgval = 0
    for j in range(4):
        bkgval += hist_content_B[int(4*ibin + j)]
    if DEBUG: print(properties['median'], ' bhist = ', bkgval)
    gr_1sig_r.SetPoint(gr_i,xval,properties['median']-bkgval)
    gr_2sig_r.SetPoint(gr_i,xval,properties['median']-bkgval)
    gr_1sig_r.SetPointError(gr_i,xerr,xerr,properties['median']-properties['down1sigma'],properties['up1sigma']-properties['median'])
    gr_2sig_r.SetPointError(gr_i,xerr,xerr,properties['median']-properties['down2sigma'],properties['up2sigma']-properties['median'])
    gr_i += 1
gr_1sig_r.SetFillColorAlpha(ROOT.kGreen+1, 0.4)
gr_2sig_r.SetFillColorAlpha(ROOT.kOrange, 0.4)
gr_1sig_r.SetLineWidth(0)
gr_2sig_r.SetLineWidth(0)

# Plotting starts from here
ROOT.gStyle.SetOptStat(0)
CMS = ""
CMS_lumi.lumi_sqrtS = "138 fb^{-1} (13 TeV) + 62 fb^{-1} (13.6 TeV)"
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "        " + CMS
CMS_lumi.cmsTextSize = 0.7
CMS_lumi.lumiTextSize = 0.4

x.setBins(int(highx-lowx))
plot = x.frame()
plot.SetTitle("")
plot_c = x_c.frame()
plot2 = x.frame()
plot2_c = x_c.frame()
plot2.SetTitle("")

combine_hist.plotOn( plot, ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2), ROOT.RooFit.Name('hist'), ROOT.RooFit.MarkerSize(0.9))
combine_B_hist.plotOn(plot_c,  ROOT.RooFit.LineColor(6), ROOT.RooFit.LineWidth(2), ROOT.RooFit.Name('bonly_model'), ROOT.RooFit.LineStyle(2), ROOT.RooFit.DrawOption("C"), ROOT.RooFit.DataError(2), ROOT.RooFit.XErrorSize(0))
combine_SB_hist.plotOn(plot_c,  ROOT.RooFit.LineColor(2), ROOT.RooFit.LineWidth(2), ROOT.RooFit.Name('model'), ROOT.RooFit.DrawOption("C"), ROOT.RooFit.DataError(2), ROOT.RooFit.XErrorSize(0))

sig_hist.plotOn( plot2, ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2), ROOT.RooFit.MarkerSize(0.9))
combine_S_hist.plotOn(plot2_c,  ROOT.RooFit.LineColor(2), ROOT.RooFit.LineWidth(2), ROOT.RooFit.DrawOption("C"), ROOT.RooFit.DataError(2), ROOT.RooFit.XErrorSize(0))#, ROOT.RooFit.Normalization(ntot_weighted, ROOT.RooAbsReal.NumEvent))

can = ROOT.TCanvas("c1","c1", 700, 500)
can.Divide(1,2)
can.cd(1)
ROOT.gPad.SetPad(0.0, 0.4, 1.0, 0.99)
ROOT.gPad.SetBottomMargin(0)
ROOT.gPad.SetTopMargin(0.15)
plot.GetYaxis().SetTitleSize(0.06)
plot.GetYaxis().SetTitleOffset(0.5)
plot.GetYaxis().SetLabelSize(0.06)
plot.GetYaxis().SetTitle('S/(S+B) weighted events/GeV')
plot.SetMinimum(0.0)
plot.GetYaxis().ChangeLabel(1, -1, 0, -1, -1, -1, "")
plot.Draw()
plot_c.Draw("same")
ROOT.gPad.Update()

leg = ROOT.TLegend(0.6,.35,.9,.8)
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetTextFont(42)
leg.SetTextSize(0.06)
leg.AddEntry('hist',' Data',"LP")
leg.AddEntry('model',' S+B model',"L")
leg.AddEntry('bonly_model',' B only model',"L")
leg.AddEntry(gr_1sig_r," 95% CI","F")
leg.AddEntry(gr_2sig_r," 68% CI","F")
leg.AddEntry(0, "", "")
leg.AddEntry(0, "#mu = %.2f ^{#plus 0.52}_{#minus 0.61}"%1.10, "")
leg.Draw("SAME")
#latex = ROOT.TLatex()
#latex.SetTextSize(0.06)
#latex.SetTextAlign(13)
#latex.SetTextFont(50)
#latex.DrawLatexNDC(.65, .45, "#mu = %.2f ^{#plus 0.52}_{#minus 0.61}"%1.10)
#latex.SetTextFont(42)
#latex.SetTextSize(0.06)
CMS_lumi.CMS_lumi(can, 0, 0)

can.cd(2)
ROOT.gPad.SetPad(0.0, 0.1, 1.0, 0.4)
ROOT.gPad.SetTopMargin(0)
ROOT.gPad.SetBottomMargin(0.33)
h_range = x.createHistogram(f"h_range",ROOT.RooFit.Binning(65,x.getMin(),x.getMax()))
h_range.SetTitle("")
h_range.GetXaxis().SetTitleOffset(0.9)
h_range.GetXaxis().SetTitleSize(0.17) #15
h_range.GetXaxis().SetTitle('m_{#font[12]{ll}\gamma}(GeV)')
h_range.GetXaxis().SetLabelSize(0.14) #12
h_range.GetYaxis().SetLabelSize(0.10)
h_range.GetYaxis().SetTitleOffset(0.2)
h_range.GetYaxis().SetTitleSize(0.12)
h_range.GetYaxis().SetNdivisions(4)
h_range.GetYaxis().SetTitle("Data - B model ")
yrange = findYRange(x, sig_hist)
h_range.GetYaxis().SetRangeUser(-yrange, yrange)
h_range.Draw()
gr_2sig_r.Draw("LE3SAME")
gr_1sig_r.Draw("LE3SAME")
plot2.Draw("same")
line = ROOT.TLine( x.getMin(), 0, x.getMax(), 0)
line.SetLineColor(ROOT.kBlack)
line.SetLineStyle(7)
line.SetLineWidth(2)
line.Draw("same")
plot2_c.Draw("same")

can.SaveAs('plots/Money_plot_update_'+args.cat + '.pdf')
can.SaveAs('plots/Money_plot_update_'+args.cat + '.root')
