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

def goodness(pdfClass, histogram,  e_type = "Poisson", eps = 0.1, n_bins = 260, className = "Default"):
    if e_type == "Poisson": error = ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson)
    elif e_type == "SumW2": error = ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2)
    good_PV = []
    chi2_=[]
    for entry in pdfClass:
        stat =  ROOT.RooChi2Var("stat_goodness", "goodness test", entry.pdf,  histogram, error)
        Minimizer_Chi2(stat, -1, 100, False, 0)
        r = Minimizer_Chi2(stat, -1, eps, False, 0)
        if r.status() != 0: 
                print(entry.pdf.GetName(), " Minimization fails!")
        entry.checkBond()
        chi2_.append(stat.getVal())
        good_PV.append(ROOT.Math.chisquared_cdf_c(stat.getVal(), n_bins - r.floatParsFinal().getSize()))
    print(className, " goodness = ", good_PV)
    print(className, " Chi2 = ", chi2_)

    

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
    plotClass(x, histogram, bern2_model.pdf, title="Bern2 " + cat, output_dir="plots/", sideBand = True, fitRange = range_)
    plotClass(x, histogram, bern3_model.pdf, title="Bern3 " + cat, output_dir="plots/", sideBand = True, fitRange = range_)
    plotClass(x, histogram, bern4_model.pdf, title="Bern4 " + cat, output_dir="plots/", sideBand = True, fitRange = range_)
    plotClass(x, histogram, bern5_model.pdf, title="Bern5 " + cat, output_dir="plots/", sideBand = True, fitRange = range_)

    multiPlotClass(x, histogram, [bern2_model, bern3_model, bern4_model, bern5_model], title="bern multi " + cat, output_dir="plots/", sideBand=True, fitRange= range_)

def singleFTestSidebandNLL(x, pdfList, histogram, cat = '', eps = 0.1, offset = True, strategy = 0, range_= "", calssName = "Default"):
    stats = []
    fitres = []
    for entry in pdfList:
        entry.reset()
        nll_= ROOT.RooNLLVar("stat_" + entry.pdf.GetName(), "stat " + entry.pdf.GetName(), entry.pdf,  histogram, ROOT.RooFit.Range(range_))
        stats.append(nll_)
        Minimizer_NLL(nll_, -1, 100, False, strategy)
        r_=Minimizer_NLL(nll_, -1, eps, offset, strategy)
        r_.Print("V")
        fitres.append(r_)
        entry.checkBond()
        plotClass(x, histogram, entry.pdf, title=entry.pdf.GetName() + "_" + cat, output_dir="plots/", sideBand = True, fitRange = range_)

    multiPlotClass(x, histogram, pdfList, title=calssName+"_multi_" + cat, output_dir="plots/", sideBand=True, fitRange= range_)
    print("NLL = ", [ele.getVal() for ele in stats])
    fs = []
    dif = []
    for i in range(len(stats) - 1):
        fs.append(ROOT.Math.chisquared_cdf_c(2*(stats[i].getVal() - stats[i+1].getVal()), fitres[i+1].floatParsFinal().getSize() - fitres[i].floatParsFinal().getSize()))
        dif.append(fitres[i+1].floatParsFinal().getSize() - fitres[i].floatParsFinal().getSize())
    print("DOF diff = ", dif)
    print("P-value = ", fs)

        

# ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)
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

# Define PDF classes
bern2_model_good = Bern2Class(x, mu_gauss, CAT, 10, 0.3, 10, 7., 105.)
bern3_model_good = Bern3Class(x, mu_gauss, CAT, 10, 0.3, 50, 7., 105.)
bern4_model_good = Bern4Class(x, mu_gauss, CAT, 10, 0.3, 50, 7., 105.)
bern5_model_good = Bern5Class(x, mu_gauss, CAT, 10, 0.3, 50, 7., 105.)
bern_list = [bern2_model_good, bern3_model_good, bern4_model_good,bern5_model_good]

pow1_model = Pow1Class(x, mu_gauss, "bin1")
pow2_model = Pow2Class(x, mu_gauss, "bin1")
pow3_model = Pow3Class(x, mu_gauss, "bin1")
pow_list = [pow1_model, pow2_model]

exp1_model = Exp1Class(x, mu_gauss, "bin1")
exp2_model = Exp2Class(x, mu_gauss, "bin1")
exp3_model = Exp3Class(x, mu_gauss, "bin1")
exp_list = [exp1_model, exp2_model, exp3_model]

lau1_model = Lau1Class(x, mu_gauss, "bin1")
lau2_model = Lau2Class(x, mu_gauss, "bin1")
lau3_model = Lau3Class(x, mu_gauss, "bin1")
lau_list = [lau1_model, lau2_model, lau3_model]

modg_model = ModGausClass(x, "bin1", lowx, lowx+65)

# Goodness of fit test
#goodness(bern_list, reader.data_hist_untagged1_bkg,  e_type = "Poisson", eps = 0.1, n_bins = 260, className="Bern")
goodness(pow_list, reader.data_hist_untagged1_bkg,  e_type = "SumW2", eps = 0.1, n_bins = 260, className="Bern")
# goodness(exp_list, reader.data_hist_untagged1_bkg,  e_type = "Poisson", eps = 0.1, n_bins = 260, className="Bern")
# goodness(lau_list, reader.data_hist_untagged1_bkg,  e_type = "Poisson", eps = 0.1, n_bins = 260, className="Bern")

# F Tset
#singleBernFTest(x, mu_gauss, reader.data_u1, CAT, args.method, "Poisson", eps = 0.1, offset = True, strategy = 0, range_ = "left,right", n_bins = 220)
singleFTestSidebandNLL(x, pow_list, reader.data_u1, cat = CAT, eps = 0.1, offset = True, strategy = 0, range_= "left,right", calssName = "Pow")
# singleFTestSidebandNLL(x, exp_list, reader.data_u1, cat = CAT, eps = 0.1, offset = True, strategy = 0, range_= "left,right", calssName = "Exp")
# singleFTestSidebandNLL(x, lau_list, reader.data_u1, cat = CAT, eps = 0.1, offset = True, strategy = 0, range_= "left,right", calssName = "Lau")





