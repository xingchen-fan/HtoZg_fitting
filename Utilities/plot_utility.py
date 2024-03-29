import ROOT
from bkg_functions_class import *

def plotClass (x, datahist, pdf, title="Histogram", output_dir="plots/"):
    ROOT.gStyle.SetOptStat(0)
    h_hist = datahist.createHistogram("h_hist", x, ROOT.RooFit.Binning(65))
    model_hist = pdf.generateBinned(x, datahist.sumEntries(), True).createHistogram("model_hist", x, ROOT.RooFit.Binning(65))

    ratio = ROOT.TH1D("ratio", "ratio", 65, x.getMin(), x.getMax())
    ratio.Divide(h_hist, model_hist)
    ratio.SetMarkerColor(ROOT.kRed)
    ratio.SetMarkerStyle(8)
    ratio.SetMarkerSize(1)
    ratio.SetTitle("")
    ratio.GetYaxis().SetRangeUser(0.5, 1.5)
    ratio.GetXaxis().SetTitleOffset(-3)
    ratio.GetXaxis().SetTitleSize(0.15)
    ratio.GetXaxis().SetTitle("m_{llg}")
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

    x.setBins(65)
    show_hist = ROOT.RooDataHist("show_hist", "show_hist", x, datahist)
    ROOT.gStyle.SetOptStat(0)
    can = ROOT.TCanvas("","", 500, 500)
    #pad1 = ROOT.TPad("pad1","",0,0.15,1,1)
    #pad2 = ROOT.TPad("pad2","",0,0,1,0.2)
    #pad1.Draw()
    #pad2.Draw()
    #pad1.cd()
    plot = x.frame()
    plot.SetTitle(title)
    can.Divide(1,2)
    can.cd(1)
    ROOT.gPad.SetPad(0.0, 0.3, 1.0, 1.0)
    ROOT.gPad.SetBottomMargin(0)
    show_hist.plotOn( plot, ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))
    pdf.plotOn(plot,  ROOT.RooFit.LineColor(2), ROOT.RooFit.LineWidth(2))
    plot.Draw()
    can.cd(2)
    ROOT.gPad.SetPad(0.0, 0.1, 1.0, 0.3)
    ROOT.gPad.SetTopMargin(0)
   # ROOT.gPad.SetBottomMargin(0.05)
    ratio.Draw()
    line.Draw("same")
    can.SaveAs(output_dir + title+".pdf")
    x.setBins(260)

    
    
