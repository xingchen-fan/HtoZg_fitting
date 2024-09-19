import ROOT
import sys
import os
sys.path.append(os.path.abspath("../CMS_plotter/"))
from bkg_functions_class import *
import CMS_lumi

def plotClass (x, datahist, pdf, SBpdf, title="Histogram", output_dir="plots/", sideBand = False, fitRange = ''):
    ROOT.gStyle.SetOptStat(0)
    CMS_lumi.lumi_sqrtS = "13 TeV"
    CMS_lumi.writeExtraText = 1
    CMS_lumi.extraText = "          Preliminary"
    CMS_lumi.cmsTextSize = 0.4
    CMS_lumi.lumiTextSize = 0.5
    if sideBand:
        datahist = datahist.reduce(ROOT.RooFit.CutRange(fitRange))
    tot = datahist.sumEntries()
    h_hist = datahist.createHistogram("h_hist", x, ROOT.RooFit.Binning(65))
    model_hist = SBpdf.generateBinned(x, tot, True).createHistogram("model_hist", x, ROOT.RooFit.Binning(65))
    norm_hist = pdf.generateBinned(x, 100000, True)
    if sideBand: sum_SB = norm_hist.sumEntries('x', fitRange)
    else: sum_SB = norm_hist.sumEntries()
    
    for i in range(65):
        model_hist.SetBinError(i+1, 0)
    ratio = ROOT.TH1D("ratio", "ratio", 65, x.getMin(), x.getMax())
    ratio.Divide(h_hist, model_hist)
    ratio.SetMarkerColor(ROOT.kBlack)
    ratio.SetMarkerStyle(8)
    ratio.SetMarkerSize(1)
    ratio.SetTitle("")
    ratio.GetYaxis().SetRangeUser(0.5, 1.5)
    ratio.GetXaxis().SetTitleOffset(0.4)
    ratio.GetXaxis().SetTitleSize(0.15)
    ratio.GetXaxis().SetTitle("m_{#font[12]{ll}\gamma}")
    ratio.GetXaxis().SetLabelSize(0.13)
    ratio.GetYaxis().SetTitleSize(0.1)
    ratio.GetYaxis().SetLabelSize(0.1)
    ratio.GetYaxis().SetTitleOffset(0.2)
    ratio.GetYaxis().SetTitleSize(0.15)
    ratio.GetYaxis().SetTitle("Ratio")
    
    line = ROOT.TLine( x.getMin(), 1, x.getMax(), 1)
    line.SetLineColor(ROOT.kBlack)
    line.SetLineStyle(7)
    line.SetLineWidth(2)

    leg = ROOT.TLegend(.5,.7,.9,.9)
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(50)
    leg.SetTextSize(0.035)
    leg.AddEntry(pdf,pdf.GetName(),"")

    x.setBins(65)
    show_hist = ROOT.RooDataHist("show_hist", "show_hist", x, datahist)
    #show_normhist = ROOT.RooDataHist("show_normhist", "show_normhist", x, normhist)
    ROOT.gStyle.SetOptStat(0)
    can = ROOT.TCanvas("","", 500, 500)
    #pad1 = ROOT.TPad("pad1","",0,0.15,1,1)
    #pad2 = ROOT.TPad("pad2","",0,0,1,0.2)
    #pad1.Draw()
    #pad2.Draw()
    #pad1.cd()
    plot = x.frame()
    plot.SetTitle("")
    can.Divide(1,2)
    can.cd(1)
    ROOT.gPad.SetPad(0.0, 0.3, 1.0, 0.99)
    ROOT.gPad.SetBottomMargin(0)
    ROOT.gPad.SetTopMargin(0.15)
    #show_normhist.plotOn( plot,ROOT.RooFit.LineColor(0), ROOT.RooFit.MarkerColor(0))
    show_hist.plotOn( plot, ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))
    if sideBand: pdf.plotOn(plot,  ROOT.RooFit.LineColor(2), ROOT.RooFit.LineWidth(2), ROOT.RooFit.Range('full'), ROOT.RooFit.Normalization(100000./sum_SB, ROOT.RooAbsReal.Relative))#ROOT.RooFit.NormRange(fitRange))
    else: pdf.plotOn(plot,  ROOT.RooFit.LineColor(2), ROOT.RooFit.LineWidth(2))
    plot.Draw()
    leg.Draw("SAME")
    CMS_lumi.CMS_lumi(can, 0, 0)
    can.cd(2)
    ROOT.gPad.SetPad(0.0, 0.1, 1.0, 0.3)
    ROOT.gPad.SetTopMargin(0)
    ROOT.gPad.SetBottomMargin(0.2)
    ratio.Draw()
    line.Draw("same")
    can.SaveAs(output_dir + title+".pdf")
    x.setBins(260)

def multiPlotClass(x, datahist, classList, title="Histogram", output_dir="plots/", sideBand = False, fitRange = ''):
    ROOT.gStyle.SetOptStat(0)
    CMS_lumi.lumi_sqrtS = "13 TeV"
    CMS_lumi.writeExtraText = 1
    CMS_lumi.extraText = "          Preliminary"
    can2 = ROOT.TCanvas("c2","c2", 500, 500)
    can2.cd()
    plot2 = x.frame()
    x.setBins(65)
    show_hist = ROOT.RooDataHist("show_multi_hist", "show_multi_hist", x, datahist)
    if sideBand: show_hist = show_hist.reduce(ROOT.RooFit.CutRange(fitRange))
    show_hist.plotOn(plot2,ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))
    for i,entry in enumerate(classList):
        norm_hist = entry.pdf.generateBinned(x, 100000, True)
        if sideBand: sum_SB = norm_hist.sumEntries('x', fitRange)
        else: sum_SB = norm_hist.sumEntries()

        if sideBand: entry.pdf.plotOn(plot2, ROOT.RooFit.LineColor(i+2), ROOT.RooFit.LineWidth(3), ROOT.RooFit.Name(entry.pdf.GetName()), ROOT.RooFit.Range('full'), ROOT.RooFit.Normalization(100000./sum_SB, ROOT.RooAbsReal.Relative))#ROOT.RooFit.NormRange(fitRange))
        else: entry.pdf.plotOn(plot2, ROOT.RooFit.LineColor(i+2), ROOT.RooFit.LineWidth(3), ROOT.RooFit.Name(entry.pdf.GetName()))
    plot2.Draw()
    plot2.GetXaxis().SetTitle("m_{#font[12]{ll}\gamma} (GeV)")
    plot2.GetYaxis().SetTitleOffset(1)
    plot2.SetTitle("")
    leg = ROOT.TLegend(.55,.7,.9,.9)
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.030)
    for entry in classList:
        leg.AddEntry(entry.pdf.GetName(), entry.pdf.GetName(),"L")
    leg.Draw("same")
    CMS_lumi.CMS_lumi(can2, 0, 0)
    #    can2.Draw()
    can2.SaveAs(output_dir + title + ".pdf")
    x.setBins(260)
    
