import ROOT
import sys
import os
sys.path.append(os.path.abspath("../CMS_plotter/"))
from bkg_functions_class import *
import CMS_lumi

def plotClass (x, datahist, pdf, SBpdf, title="Histogram", output_dir="plots/", sideBand = False, fitRange = '', note = "", CMS = "Preliminary", fullRatio = False, leftSpace=False, bins = 1, toy = False, bkg_normSF = 1):
    #-------------------------------------------------------------
    # In the case of toy plotting, SBpdf is the background-only pdf
    #--------------------------------------------------------------
    
    ROOT.gStyle.SetOptStat(0)
    CMS_lumi.lumi_sqrtS = "137.61 fb^{-1} (13 TeV) + 62.32 fb^{-1} (13.6 TeV)"
    CMS_lumi.writeExtraText = 1
    CMS_lumi.extraText = "      " + CMS
    CMS_lumi.cmsTextSize = 0.5
    CMS_lumi.lumiTextSize = 0.3
    if sideBand:
        datahist = datahist.reduce(ROOT.RooFit.CutRange(fitRange))
        model_hist = SBpdf.generateBinned(x, datahist.sumEntries(), True).createHistogram("model_hist", x, ROOT.RooFit.Binning(int(bins * (x.getMax() - x.getMin()))))
    else:
        model_hist = pdf.generateBinned(x, datahist.sumEntries(), True).createHistogram("model_hist", x, ROOT.RooFit.Binning(int(bins * (x.getMax() - x.getMin()))))
        
    h_hist = datahist.createHistogram("h_hist", x, ROOT.RooFit.Binning(int(bins * (x.getMax() - x.getMin()))))
    if toy: h_hist.Sumw2(False)
    norm_hist = pdf.generateBinned(x, 100000, True)
    if sideBand: sum_SB = norm_hist.sumEntries('x', fitRange)
    else: sum_SB = norm_hist.sumEntries()
    
    for i in range(int(bins * (x.getMax() - x.getMin()))):
        model_hist.SetBinError(i+1, 0)
    ratio = ROOT.TH1D("ratio", "ratio", int(bins * (x.getMax() - x.getMin())), x.getMin(), x.getMax())
    ratio.Divide(h_hist, model_hist)
    ratio.SetMarkerColor(ROOT.kBlack)
    ratio.SetMarkerStyle(8)
    ratio.SetMarkerSize(1)
    ratio.SetLineColor(ROOT.kBlack)
    ratio.SetTitle("")
    if fullRatio:
        minimum = 99999
        for i in range(int(bins * (x.getMax() - x.getMin()))):
            bincont = ratio.GetBinContent(i+1)
            if bincont != 0 and bincont < minimum:
                minimum = bincont
        ratio.GetYaxis().SetRangeUser(minimum, ratio.GetBinContent(ratio.GetMaximumBin()))
    else:
        ratio.GetYaxis().SetRangeUser(0.9, 1.1)
    ratio.GetXaxis().SetTitleOffset(0.9)
    ratio.GetXaxis().SetTitleSize(0.17)
    ratio.GetXaxis().SetTitle("m_{#font[12]{ll}\gamma}(GeV)")
    ratio.GetXaxis().SetLabelSize(0.17)
    ratio.GetYaxis().SetLabelSize(0.12)
    ratio.GetYaxis().SetTitleOffset(0.2)
    ratio.GetYaxis().SetTitleSize(0.17)
    ratio.GetYaxis().SetTitle("Ratio")
    ratio.GetYaxis().SetNdivisions(4)
    line = ROOT.TLine( x.getMin(), 1, x.getMax(), 1)
    line.SetLineColor(ROOT.kBlack)
    line.SetLineStyle(7)
    line.SetLineWidth(2)

    if CMS == "Simulation":
        hist_name = "MC Sample"
    else:
        hist_name = "Data"

    
    latex = ROOT.TLatex()
    latex.SetTextSize(0.03)
    latex.SetTextAlign(13)
    x.setBins(int(bins * (x.getMax() - x.getMin())))
    show_hist = ROOT.RooDataHist("show_hist", "show_hist", ROOT.RooArgSet(x), datahist)
    #show_normhist = ROOT.RooDataHist("show_normhist", "show_normhist", x, normhist)
    ROOT.gStyle.SetOptStat(0)
    can = ROOT.TCanvas("c1","c1", 700, 500)
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
    ROOT.gPad.SetTopMargin(0.12)
    #show_normhist.plotOn( plot,ROOT.RooFit.LineColor(0), ROOT.RooFit.MarkerColor(0))
    if toy:
        show_hist.plotOn( plot, ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson), ROOT.RooFit.Name('hist'))
        SBpdf.plotOn(plot,  ROOT.RooFit.LineColor(4), ROOT.RooFit.LineWidth(2), ROOT.RooFit.Name('bonly_model'), ROOT.RooFit.Normalization(bkg_normSF))
    else: show_hist.plotOn( plot, ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2), ROOT.RooFit.Name('hist'))
    if sideBand: pdf.plotOn(plot,  ROOT.RooFit.LineColor(2), ROOT.RooFit.LineWidth(2), ROOT.RooFit.Range('full'), ROOT.RooFit.Normalization(100000./sum_SB, ROOT.RooAbsReal.Relative), ROOT.RooFit.Name('model'))#ROOT.RooFit.NormRange(fitRange))
    else: pdf.plotOn(plot,  ROOT.RooFit.LineColor(2), ROOT.RooFit.LineWidth(2), ROOT.RooFit.Name('model'))
    plot.GetYaxis().SetTitleSize(0.05)
    plot.GetYaxis().SetTitleOffset(0.8)
    plot.GetYaxis().SetLabelSize(0.05)
    plot.Draw()
    if leftSpace:
        leg = ROOT.TLegend(.15,.7,.55,.8)
    else:
        leg = ROOT.TLegend(.5,.7,.9,.8)
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(50)
    leg.SetTextSize(0.05)
    leg.AddEntry('hist',hist_name,"LP")
    leg.AddEntry('model',pdf.GetName(),"L")
    leg.Draw("SAME")
    CMS_lumi.CMS_lumi(can, 0, 0)
    if leftSpace:
        latex.DrawLatexNDC(.15, .7, note)
    else:
        latex.DrawLatexNDC(.6, .7, note)
    can.cd(2)
    ROOT.gPad.SetPad(0.0, 0.1, 1.0, 0.3)
    ROOT.gPad.SetTopMargin(0)
    ROOT.gPad.SetBottomMargin(0.33)
    ratio.Draw()
    line.Draw("same")
    can.SaveAs(output_dir + title+".pdf")
    x.setBins(int(4*(x.getMax() - x.getMin())))

def multiPlotClass(x, datahist, classList, title="Histogram", output_dir="plots/", sideBand = False, fitRange = '', best_index = 0, CMS = "Preliminary", fullRatio = True, bestLabel = False, bins = 1, ratio_range=[0.9,1.1]):
    ROOT.gStyle.SetOptStat(0)
    CMS_lumi.lumi_sqrtS = "137.61 fb^{-1} (13 TeV) + 62.32 fb^{-1} (13.6 TeV)"
    CMS_lumi.writeExtraText = 1
    CMS_lumi.extraText = "      " + CMS
    CMS_lumi.cmsTextSize = 0.5
    CMS_lumi.lumiTextSize = 0.3
    nbins = int(bins * (x.getMax() - x.getMin()))
    
    #Define Canvas
    can2 = ROOT.TCanvas("c2","c2", 700, 500)
    can2.Divide(1,2)
    can2.cd(1)
    ROOT.gPad.SetPad(0.0, 0.3, 1.0, 0.99)
    ROOT.gPad.SetBottomMargin(0)
    ROOT.gPad.SetTopMargin(0.12)

    #Create plot frame
    plot2 = x.frame()
    x.setBins(nbins)

    #Plot histogram
    show_hist = ROOT.RooDataHist("show_multi_hist", "show_multi_hist", ROOT.RooArgSet(x), datahist)
    if sideBand: show_hist = show_hist.reduce(ROOT.RooFit.CutRange(fitRange))
    show_hist.plotOn(plot2, ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2), ROOT.RooFit.Name('hist'))

    #Plot functions
    color_index = 1
    best_color = 0
    for i,entry in enumerate(classList):
        norm_hist = entry.pdf.generateBinned(x, 100000, True)
        color_index +=1
        if color_index == 10: color_index = 20
        if color_index > 20: color_index += 5
        if i == best_index: best_color = color_index
        if sideBand: sum_SB = norm_hist.sumEntries('x', fitRange)
        else: sum_SB = norm_hist.sumEntries()

        if sideBand: entry.pdf.plotOn(plot2, ROOT.RooFit.LineColor(color_index), ROOT.RooFit.LineWidth(2), ROOT.RooFit.Name(entry.pdf.GetName()), ROOT.RooFit.Range('full'), ROOT.RooFit.Normalization(100000./sum_SB, ROOT.RooAbsReal.Relative))#ROOT.RooFit.NormRange(fitRange))
        else: entry.pdf.plotOn(plot2, ROOT.RooFit.LineColor(color_index), ROOT.RooFit.LineWidth(2), ROOT.RooFit.Name(entry.pdf.GetName()))
    if sideBand:
        
        #Ctreate model histograms for ratio plot
        model_hist = classList[best_index].SBpdf.generateBinned(x, show_hist.sumEntries(), True).createHistogram("model_hist", x, ROOT.RooFit.Binning(nbins))
    else:
        model_hist = classList[best_index].pdf.generateBinned(x, show_hist.sumEntries(), True).createHistogram("model_hist", x, ROOT.RooFit.Binning(nbins))
    plot2.Draw()
    #plot2.GetXaxis().SetTitle("m_{#font[12]{ll}\gamma} (GeV)")
    plot2.GetYaxis().SetTitleOffset(1)
    plot2.GetYaxis().SetTitleSize(0.05)
    plot2.GetYaxis().SetTitleOffset(0.8)
    plot2.GetYaxis().SetLabelSize(0.05)
    plot2.SetTitle("")
    
    if CMS == "Simulation":
        hist_name = "MC Sample"
    else:
        hist_name = "Data"
    leg = ROOT.TLegend(.55,.3,.9,.8)
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.050)
    leg.AddEntry('hist', hist_name, 'LP')
    for ind, entry in enumerate(classList):
        if bestLabel and ind == best_index:
            leg.AddEntry(entry.pdf.GetName(), entry.pdf.GetName()+' (best)',"L")
        else:
            leg.AddEntry(entry.pdf.GetName(), entry.pdf.GetName(),"L")
    leg.Draw("same")

    #Create TH1 histogram for ratio plot
    h_hist = show_hist.createHistogram("h_hist", x, ROOT.RooFit.Binning(65))
    print('xmin = ', x.getMin())

    #Create ratio plot
    ratio = ROOT.TH1D("ratio", "ratio", 65, x.getMin(), x.getMax())
    ratio.Divide(h_hist, model_hist)
    ratio.SetMarkerColor(ROOT.kBlack)
    ratio.SetMarkerStyle(8)
    ratio.SetMarkerSize(1)
    ratio.SetLineColor(ROOT.kBlack)
    ratio.SetTitle("")
    if fullRatio:
        minimum = 99999
        for i in range(int(x.getMax() - x.getMin())):
            bincont = ratio.GetBinContent(i+1)
            if bincont != 0 and bincont < minimum:
                minimum = bincont
        ratio.GetYaxis().SetRangeUser(minimum, ratio.GetBinContent(ratio.GetMaximumBin()))
    else:
        ratio.GetYaxis().SetRangeUser(ratio_range[0], ratio_range[1])
    #Ratio plot seetings
    ratio.GetXaxis().SetTitleOffset(0.9)
    ratio.GetXaxis().SetTitleSize(0.17)
    ratio.GetXaxis().SetTitle("m_{#font[12]{ll}\gamma}/(GeV)")
    ratio.GetXaxis().SetLabelSize(0.17)
    ratio.GetYaxis().SetTitleSize(0.2)
    ratio.GetYaxis().SetLabelSize(0.12)
    ratio.GetYaxis().SetTitleOffset(0.2)
    ratio.GetYaxis().SetTitleSize(0.17)
    ratio.GetYaxis().SetNdivisions(4)
    ratio.GetYaxis().SetTitle("Ratio")

    line = ROOT.TLine( x.getMin(), 1, x.getMax(), 1)
    line.SetLineColor(ROOT.kBlack)
    line.SetLineStyle(7)
    line.SetLineWidth(2)

    CMS_lumi.CMS_lumi(can2, 0, 0)
    #can2.Draw()
    can2.cd(2)
    ROOT.gPad.SetPad(0.0, 0.1, 1.0, 0.3)
    ROOT.gPad.SetTopMargin(0)
    ROOT.gPad.SetBottomMargin(0.33)
    ratio.Draw()
    line.Draw("same")

    can2.SaveAs(output_dir + title + ".pdf")
    x.setBins(260)
    
