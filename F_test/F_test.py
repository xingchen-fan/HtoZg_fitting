import ROOT
import os
import sys
import argparse
sys.path.append(os.path.abspath("../Utilities/"))
sys.path.append(os.path.abspath("../CMS_plotter/"))
import CMS_lumi, tdrstyle
from bkg_functions_fit import *
from bkg_functions_class import *
from Xc_Minimizer import *
from plot_utility import *
from sample_reader import *
#ROOT.gInterpreter.AddIncludePath('../Utilities/HZGRooPdfs.h')
ROOT.gSystem.Load('../Utilities/HZGRooPdfs_cxx.so')

parser = argparse.ArgumentParser(description = "F test method (NLL or Chi2)")
parser.add_argument("method")
args = parser.parse_args()
if not(args.method == "Chi2" or args.method == "NLL") :
    print("Please use the correct method.")
    sys.exit(1)


# def singleBernFTest(x, gauss_mu, histogram, cat = "", method = "Chi2", e_type = "Poisson", offset = False):
#     bern2_model = Bern2Minization(x, gauss_mu, histogram, method, e_type, cat, 10, 0.3, 10, 7., 105., \
#                                   printLevel = -1, eps = 1, offSet = offset, strategy = 0)
#     bern2_model[0].Print("v")
#     bern3_model = Bern3Minization(x, gauss_mu, histogram, method, e_type, cat, 10, 0.3, 50, 7., 105.,\
#                                   printLevel = -1, eps = 1, offSet = offset, strategy = 0)
#     bern3_model[0].Print("v")
#     bern4_model = Bern4Minization(x, gauss_mu, histogram, method, e_type, cat, 10, 0.3, 50, 7., 105.,\
#                                   printLevel = -1, eps = 1, offSet = offset, strategy = 0)
#     bern4_model[0].Print("v")
#     bern5_model = Bern5Minization(x, gauss_mu, histogram, method, e_type, cat, 10, 0.3, 50, 7., 105.,\
#                                   printLevel = -1, eps = 1, offSet = offset, strategy = 0)
#     bern5_model[0].Print("v")

    
#     stat = [bern2_model[1], bern3_model[1], bern4_model[1], bern5_model[1]]
#     fs = []
#     if method == "Chi2":
#         for i in range(len(stat) - 1):
#             fs.append(ROOT.Math.fdistribution_cdf_c((stat[i] - stat[i+1])*(260-5-i)/stat[i+1], 1, 260-5-i))
#     elif method == "NLL":
#         for i in range(len(stat) - 1):
#             fs.append(ROOT.Math.chisquared_cdf_c(2*(stat[i] - stat[i+1]), 1))

#     print(method, " = ", stat)
#     print("P-value = ", fs)

def goodness(pdfClass, histogram,  e_type = "Poisson", eps = 0.1, n_bins = 260):
    if e_type == "Poisson": error = ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson)
    elif e_type == "SumW2": error = ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2)
    stat =  ROOT.RooChi2Var("stat_goodness", "goodness test", pdfClass.pdf,  histogram, error)
    Minimizer_Chi2(stat, -1, 100, False, 0)
    r = Minimizer_Chi2(stat, -1, eps, False, 0)
    pdfClass.checkBond()
    print(pdfClass.pdf.GetName(), " goodness = ", ROOT.Math.chisquared_cdf_c(stat.getVal(), n_bins - r.floatParsFinal().getSize()))

    

def singleBernFTest(x, gauss_mu, histogram, cat = "", method = "Chi2", e_type = "Poisson", eps = 0.1, offset = False, strategy = 0, range_ = "", n_bins = 260):
    bern2_model = Bern2Class(x, gauss_mu, cat, 10, 0.3, 10, 7., 105.)
    bern3_model = Bern3Class(x, gauss_mu, cat, 10, 0.3, 50, 7., 105.)
    bern4_model = Bern4Class(x, gauss_mu, cat, 10, 0.3, 50, 7., 105.)
    bern5_model = Bern5Class(x, gauss_mu, cat, 10, 0.3, 50, 7., 105.)
    if e_type == "Poisson": error = ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson)
    elif e_type == "SumW2": error = ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2)
    #cuthistogram = histogram.reduce(ROOT.RooFit.CutRange("left"))
    if method == "Chi2": 
        stat1 = ROOT.RooChi2Var("stat_bern2_" + cat, "stat bern2 " + cat, bern2_model.pdf,  histogram, error)
        stat2 = ROOT.RooChi2Var("stat_bern3_" + cat, "stat bern3 " + cat, bern3_model.pdf,  histogram, error)
        stat3 = ROOT.RooChi2Var("stat_bern4_" + cat, "stat bern4 " + cat, bern4_model.pdf,  histogram, error)
        stat4 = ROOT.RooChi2Var("stat_bern5_" + cat, "stat bern5 " + cat, bern5_model.pdf,  histogram, error)
    elif method == "NLL":
        stat1 = ROOT.RooNLLVar("stat_bern2_" + cat, "stat bern2 " + cat, bern2_model.pdf,  histogram, ROOT.RooFit.Range(range_))
        stat2 = ROOT.RooNLLVar("stat_bern3_" + cat, "stat bern3 " + cat, bern3_model.pdf,  histogram, ROOT.RooFit.Range(range_))
        stat3 = ROOT.RooNLLVar("stat_bern4_" + cat, "stat bern4 " + cat, bern4_model.pdf,  histogram, ROOT.RooFit.Range(range_))
        stat4 = ROOT.RooNLLVar("stat_bern5_" + cat, "stat bern5 " + cat, bern5_model.pdf,  histogram, ROOT.RooFit.Range(range_))

    stats = [stat1, stat2, stat3, stat4]
    if method == "Chi2": 
        for entry in stats:
            print(entry.GetTitle())
            Minimizer_Chi2(entry, -1, 100, False, strategy)
            r = Minimizer_Chi2(entry, -1, eps, offset, strategy)
            r.Print("V")
            # print("Cov q = ", r.covQual(), " status = ", r.status(), end="\n\n")
    elif method == "NLL": 
        for entry in stats:
            print(entry.GetTitle())
            Minimizer_NLL(entry, -1, 100, False, strategy)
            r = Minimizer_NLL(entry, -1, eps, offset, strategy)
            r.Print("V")
    output = [stat1.getVal(), stat2.getVal(), stat3.getVal(), stat4.getVal()]
    #res = bern4_model.pdf.chi2FitTo(cuthistogram, ROOT.RooFit.Range("left"), error, ROOT.RooFit.PrintLevel(-1), ROOT.RooFit.Save(True))
    #res.Print("V")
    fs = []
    bern2_model.checkBond()
    bern3_model.checkBond()
    bern4_model.checkBond()
    bern5_model.checkBond()
    if method == "Chi2":
        for i in range(len(output) - 1):
            fs.append(ROOT.Math.fdistribution_cdf_c((output[i] - output[i+1])*(260-5-i)/output[i+1], 1, 260-5-i))
    elif method == "NLL":
        for i in range(len(output) - 1):
            fs.append(ROOT.Math.chisquared_cdf_c(2*(output[i] - output[i+1]), 1))

    print(method, " = ", output)
    if method == "Chi2": print("goodness = ", [ROOT.Math.chisquared_cdf_c(output[i], n_bins - 4 -i)for i in range(len(output))])
    print("P-value = ", fs)
    plotClass(x, histogram, bern2_model.pdf, title="Bern2 u1", output_dir="plots/", sideBand = True, fitRange = range_)
    plotClass(x, histogram, bern3_model.pdf, title="Bern3 u1", output_dir="plots/", sideBand = True, fitRange = range_)
    plotClass(x, histogram, bern4_model.pdf, title="Bern4 u1", output_dir="plots/", sideBand = True, fitRange = range_)
    plotClass(x, histogram, bern5_model.pdf, title="Bern5 u1", output_dir="plots/", sideBand = True, fitRange = range_)

    CMS_lumi.lumi_sqrtS = "13 TeV"
    CMS_lumi.writeExtraText = 1
    CMS_lumi.extraText = "          Preliminary"
    can2 = ROOT.TCanvas("c2","c2", 500, 500)
    can2.cd()
    plot2 = x.frame()
    x.setBins(65)
    show_hist = ROOT.RooDataHist("show_multi_hist", "show_multi_hist", x, histogram)
    show_hist = show_hist.reduce(ROOT.RooFit.CutRange(range_))
    show_hist.plotOn(plot2,ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))
    bern2_model.pdf.plotOn(plot2, ROOT.RooFit.LineColor(2), ROOT.RooFit.LineWidth(3), ROOT.RooFit.Name("Bern2"),ROOT.RooFit.NormRange(range_))
    bern3_model.pdf.plotOn(plot2, ROOT.RooFit.LineColor(3), ROOT.RooFit.LineWidth(3), ROOT.RooFit.Name("Bern3"),ROOT.RooFit.NormRange(range_))
    bern4_model.pdf.plotOn(plot2, ROOT.RooFit.LineColor(4), ROOT.RooFit.LineWidth(3),ROOT.RooFit.Name("Bern4"),ROOT.RooFit.NormRange(range_))
    bern5_model.pdf.plotOn(plot2, ROOT.RooFit.LineColor(7), ROOT.RooFit.LineWidth(3), ROOT.RooFit.Name("Bern5"), ROOT.RooFit.NormRange(range_))
    plot2.Draw()
    plot2.GetXaxis().SetTitle("m_{#font[12]{ll}\gamma} (GeV)")
    plot2.GetYaxis().SetTitleOffset(1)
    plot2.SetTitle("")
    leg = ROOT.TLegend(.7,.7,.9,.9)
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.035)
    leg.AddEntry("Bern2","Bern 2","L")
    leg.AddEntry("Bern3","Bern 3","L")
    leg.AddEntry("Bern4","Bern 4","L")
    leg.AddEntry("Bern5","Bern 5","L")
    leg.Draw("same")
    CMS_lumi.CMS_lumi(can2, 0, 0)
    #    can2.Draw()
    can2.SaveAs("bern_multi_u1.pdf")
    x.setBins(260)
    
# ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
lowx = 100.

# Define variables
x = ROOT.RooRealVar("x", "mllg", lowx, lowx + 65.)
y = ROOT.RooRealVar("y", "photon pT", 15., 1000.)
w = ROOT.RooRealVar("w", "w", -40., 40.)
bdt = ROOT.RooRealVar("bdt", "bdt", -1, 1)
year = ROOT.RooRealVar("year", "year", 2015, 2019)
lep = ROOT.RooRealVar("lep", "lep", 0, 1) #0 = electron, 1 = muon
ph_eta = ROOT.RooRealVar("ph_eta", "ph_eta", -3, 3)
nlep = ROOT.RooRealVar("nlep", "nlep", 0, 10)
njet = ROOT.RooRealVar("njet", "njet", 0, 10)

mu_gauss = ROOT.RooRealVar("mu_gauss","always 0"       ,0.)
list = [x, y, w, bdt, year, lep, ph_eta, nlep, njet]

# Expected bkg
# expbkg_u1 = u1_bkg_run2.sumEntries()
# expbkg_u2 = u2_bkg_run2.sumEntries()
# expbkg_u3 = u3_bkg_run2.sumEntries()
# expbkg_u4 =  u4_bkg_run2.sumEntries()

# Cornell MC sample dat reader and make RooDataHist
x.setBins(260)
#reader = readDat(list, "/afs/cern.ch/user/f/fanx/public/samples/")
reader = readDat(list, "../../sample/")
reader.numCheck()
reader.dataNumCheck()
# Beijing data sample root reader and make RooDataHist
#x.setBins(260)
#reader = readRoot(x, "~/beijing_sample/data.root")

x.setRange('left', lowx, 120)
x.setRange('right', 130, lowx+65)
x.setRange('full', lowx, lowx+105)

CAT = "u1"

# Goodness of fit test
bern2_model_good = Bern2Class(x, mu_gauss, CAT, 10, 0.3, 10, 7., 105.)
bern3_model_good = Bern3Class(x, mu_gauss, CAT, 10, 0.3, 50, 7., 105.)
bern4_model_good = Bern4Class(x, mu_gauss, CAT, 10, 0.3, 50, 7., 105.)
bern5_model_good = Bern5Class(x, mu_gauss, CAT, 10, 0.3, 50, 7., 105.)
pdf_list = [bern2_model_good, bern3_model_good, bern4_model_good,bern5_model_good]
for entry in pdf_list:
    goodness(entry, reader.data_hist_untagged1_bkg,  e_type = "Poisson", eps = 0.1, n_bins = 260)

singleBernFTest(x, mu_gauss, reader.data_u1, CAT, args.method, "Poisson", eps = 0.1, offset = False, strategy = 0, range_ = "left,right", n_bins = 220)
# singleBernFTest(x, mu_gauss, reader.data_hist_untagged2_bkg, "u2", args.method, "Poisson", False)
# singleBernFTest(x, mu_gauss, reader.data_hist_untagged3_bkg, "u3", args.method, "Poisson", False)
# singleBernFTest(x, mu_gauss, reader.data_hist_untagged4_bkg, "u4", args.method, "Poisson", False)




