import ROOT
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
#import argparse
sys.path.append(os.path.abspath("../Utilities/"))
sys.path.append(os.path.abspath("../CMS_plotter/"))
import CMS_lumi, tdrstyle
from bkg_functions_fit import *
from bkg_functions_class import *
from sig_functions_class import *
from Xc_Minimizer import *
from plot_utility import *
from sample_reader import *
from bias_class import *
ROOT.gInterpreter.AddIncludePath('../Utilities/HZGRooPdfs.h')
ROOT.gSystem.Load('../Utilities/HZGRooPdfs_cxx.so')
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)

# Define variables
lowx = 105.
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
"""
N = ROOT.RooRealVar("N", "Extended term", 19000, 10000, 30000);
dummy_N = ROOT.RooRealVar("dummy_N", "dummy N", 10, 0, 1000);
dummy_mu = ROOT.RooRealVar("d_mu", "a dummy mu", 125)
dummy_width = ROOT.RooRealVar("d_sigma", "a dummy width", 0.1)
dummy_sig =  ROOT.RooGaussian("dummy", "dummy for ext", x,dummy_mu, dummy_width)
"""

# Read samples
#reader = readDat(list, "/afs/cern.ch/user/f/fanx/public/samples/")
#reader = readDat(list, "../../sample/")
reader = readRoot(x, "~/beijing_sample/data.root")
#cutHist = reader.data_hist_bin1.reduce(ROOT.RooFit.CutRange('left,right'))
#print ("cut norm = ", cutHist.sumEntries())
N = reader.data_hist_bin1.sumEntries() *10
print("10 times stats = ",N)
x.setBins(260)
x.setRange('left', lowx, 120)
x.setRange('right', 130, lowx+65)
x.setRange('full', lowx, lowx+65)

# Signal model (pre-fit)
MH = ROOT.RooRealVar("MH","MH"       ,124.7, 120., 130.)
dscb_model = DSCB_Class(x, MH, "bin1", 1.78,50, 100, 50, 100,0.845, 2.36)
dscb_model.setConst(True)
N_sig = 10

# Functions to test
bern2_model_seed = Bern2Class(x, mu_gauss, "bin1_seed", 10, 0.3, 10, 7., 105.)
bern3_model_seed = Bern3Class(x, mu_gauss, "bin1_seed", 10, 0.3, 10, 3., 106.)
bern4_model_seed = Bern4Class(x, mu_gauss, "bin1_seed", 10, 0.3, 10, 3., 106.)

bern2_model = Bern2Class(x, mu_gauss, "bin1", 10, 0.3, 10, 7., 105.)
bern3_model = Bern3Class(x, mu_gauss, "bin1", 10, 0.3, 10, 3., 106.)
bern4_model = Bern4Class(x, mu_gauss, "bin1", 10, 0.3, 10, 3., 106.)
bern5_model = Bern5Class(x, mu_gauss, "bin1", 10, 0.3, 10, 3., 106.)
#extmodel = ROOT.RooExtendPdf("extmodel", "Extended model", bern2_model.pdf, N, 'full');
#extmodel = ROOT.RooAddPdf("extmodel","extmodel", ROOT.RooArgList(bern3_model.pdf, dummy_sig), ROOT.RooArgList(N, dummy_N))
#r = bern2_model.pdf.fitTo(reader.data_hist_untagged1_bkg,ROOT.RooFit.Save(True), ROOT.RooFit.PrintLevel(-1), ROOT.RooFit.SumW2Error(True))
#r.Print("v")
#ROOT.RooFit.Range('left,right'),
profile_seed = [bern2_model_seed]#, bern3_model_seed, bern4_model_seed]
profile = [bern2_model, bern3_model, bern4_model, bern5_model]

# Set best-fit values
for entry in profile_seed:
    BiasClass(entry.pdf, reader.data_hist_bin1, True, 'left,right').minimize()
    entry.checkBond()

print("Done seed PDFs")

# Define profile Fit
def profileFit(profile_, sig_model, hist, fix = False, str = 0.):
    min_nll = 0
    ind = 999
    r_sig_ = 0
    r_error_ = 0
    best_=''
    for i, ele in enumerate(profile_):
        #ele.reset()
        if (fix): 
            c1 = ROOT.RooRealVar("c1_"+ele.pdf.GetName(), "c1_"+ele.pdf.GetName(), N)
            c2 = ROOT.RooRealVar("c2_"+ele.pdf.GetName(), "c2_"+ele.pdf.GetName(), str)
        else:
            c1 = ROOT.RooRealVar("c1_"+ele.pdf.GetName(), "c1_"+ele.pdf.GetName(), N, 0, 3.*N)
            c2 = ROOT.RooRealVar("c2_"+ele.pdf.GetName(), "c2_"+ele.pdf.GetName(), 0., -500.*N_sig, 500.*N_sig)
        tot_model = ROOT.RooAddPdf("tot_model_"+ele.pdf.GetName(), "tot_model_"+ele.pdf.GetName(), ROOT.RooArgList(sig_model.pdf, ele.pdf), ROOT.RooArgList(c2, c1))
        bias = BiasClass(tot_model, hist, False)
        bias.minimize()
        if i ==0: 
            ind = 0
            min_nll=bias.corrNLL
            r_sig_ = c2.getVal()
            r_error_ = c2.getError()
            best_= ele.pdf.GetName()
        elif bias.corrNLL< min_nll: 
            ind = i
            min_nll = bias.corrNLL
            r_sig_ = c2.getVal()
            r_error_ = c2.getError()
            best_= ele.pdf.GetName()
    return [ind, min_nll, r_sig_, r_error_, best_]

def scanFitPlot(bkgclass, sig_model, hist, r_sig_, min_nll, scan_size_ = 0.1, N_scan_ = 20):
    # c1 = ROOT.RooRealVar("c1_"+ entry.pdf.GetName(), "c1_"+ entry.pdf.GetName(), N, 0, 3.*N)
    # c2 = ROOT.RooRealVar("c2_"+ entry.pdf.GetName(), "c2_"+ entry.pdf.GetName(), 0., -500.*N_sig, 500.*N_sig)
    # tot_model = ROOT.RooAddPdf("tot_model_"+ entry.pdf.GetName(), "tot_model_"+ entry.pdf.GetName(), ROOT.RooArgList(sig_model.pdf, entry.pdf), ROOT.RooArgList(c2, c1))
    # bias = BiasClass(tot_model, hist, False)
    # bias.minimize()
    # offset = bias.corrNLL
    # pll = bias.nll.createProfile(c2)

    scan_list_ = []
    # for k in range(N_scan):
    #     c2.setVal(abs(r_sig) * (k - N_scan/2) * scan_size)
    #     scan_list_.append(pll.getVal() + offset - min_nll)
    for k in range(N_scan_):
        c1 = ROOT.RooRealVar("c1_"+ bkgclass.pdf.GetName(), "c1_"+ bkgclass.pdf.GetName(), N, 0, 3.*N)
        c2 = ROOT.RooRealVar("c2_"+ bkgclass.pdf.GetName(), "c2_"+ bkgclass.pdf.GetName(), abs(r_sig_) * (k - N_scan_/2) * scan_size_)
        tot_model = ROOT.RooAddPdf("tot_model_"+ bkgclass.pdf.GetName(), "tot_model_"+ bkgclass.pdf.GetName(), ROOT.RooArgList(sig_model.pdf, bkgclass.pdf), ROOT.RooArgList(c2, c1))
        bias = BiasClass(tot_model, hist, False)
        if k == N_scan/2: 
            bias.minimize(skip_hesse = True)
        else: bias.minimize(skip_hesse = False)
        scan_list_.append(bias.corrNLL - min_nll + 0.5) # One fewer DOF
    return scan_list_

def scanFit(profile_, sig_model, hist, r_sig_, scan_size_ = 0.5):
    scan_list_ = []
    scan_all_ = []
    output_all_ = []
    # Left r < 0
    step = 0
    offset_nll = 0.
    scan = True
    while scan:
        choose = []
        for pdf_ in profile_:
            #pdf_.reset()
            c1 = ROOT.RooRealVar("c1_"+ pdf_.pdf.GetName(), "c1_"+ pdf_.pdf.GetName(), N, 0, 2.*N)
            c2 = ROOT.RooRealVar("c2_"+ pdf_.pdf.GetName(), "c2_"+ pdf_.pdf.GetName(), - abs(r_sig_) * step * scan_size_)
            tot_model = ROOT.RooAddPdf("tot_"+ pdf_.pdf.GetName(), "tot_"+ pdf_.pdf.GetName(), ROOT.RooArgList(sig_model.pdf, pdf_.pdf), ROOT.RooArgList(c2, c1))
            bias = BiasClass(tot_model, hist, False)
            if step==0: 
                bias.minimize(skip_hesse = True)
            else: bias.minimize(skip_hesse = True)
            choose.append(bias.corrNLL)
        chose = min(choose)
        if step==0: offset_nll = chose
        scan_list_.insert(0, chose)
        scan_all_.insert(0, choose)
        step += 1
        if chose - offset_nll > 2.: scan = False
        elif step > 40: scan = False
    
    # Right r > 0
    step = 1
    scan = True
    while scan:
        choose = []
        for pdf_ in profile_:
            if step == 1: pdf_.reset()
            c1 = ROOT.RooRealVar("c1_"+ pdf_.pdf.GetName(), "c1_"+ pdf_.pdf.GetName(), N, 0, 2.*N)
            c2 = ROOT.RooRealVar("c2_"+ pdf_.pdf.GetName(), "c2_"+ pdf_.pdf.GetName(), abs(r_sig_) * step * scan_size_)
            tot_model = ROOT.RooAddPdf("tot_"+ pdf_.pdf.GetName(), "tot_"+ pdf_.pdf.GetName(), ROOT.RooArgList(sig_model.pdf, pdf_.pdf), ROOT.RooArgList(c2, c1))
            bias = BiasClass(tot_model, hist, False)
            bias.minimize(skip_hesse = True)
            choose.append(bias.corrNLL)
        chose = min(choose)
        scan_list_.append(chose)
        scan_all_.append(choose)
        step += 1
        if chose - offset_nll > 2.: scan = False
        elif step > 40: scan = False

    for i in range(len(profile_)):
        output_all_.append([x[i] for x in scan_all_])
    
    dNLL_ = [x - min(scan_list_) for x in scan_list_ ]
    return dNLL_, output_all_

# Discrete profiling - Find minimum and (r_down, r_up)
# Scan points of signal_yield * scan_size around 0
N_toy = 10
# N_scan = 40
scan_size = 0.25
for entry in profile_seed:
    r_sig = []
    r_error = []
    best_list = []
    best_error = []
    bad = 0
    pull_list = []
    pull = ROOT.TH1F("pull", "pull", 80, -4, 4)
    can = ROOT.TCanvas("can", "can", 500, 500)
    for j in range(N_toy):
        scan_list = []      
        x.setBins(260)
        hist_toy = entry.pdf.generateBinned(x, ROOT.RooFit.NumEvents(N))
        list = profileFit(profile, dscb_model, hist_toy)
        
        # Method1 ##########################
        # for ele in profile:
        #     ele.reset()
        #     scan_list.append(scanFitPlot(ele, dscb_model, hist_toy, list[2], list[1], scan_size, N_scan))
        # dNLL_offset = []
        # for m in range(N_scan):
        #     dNLL_offset.append(min([scan_list[n][m] for n in range(len(profile))]))
        # dNLL = [x - min(dNLL_offset) for x in dNLL_offset]
        # xs = [scan_size*(x - N_scan/2) for x in range(N_scan)]
        # fig = plt.figure()    
        # plt.plot(xs, scan_list[0])
        # plt.plot(xs, scan_list[1])
        # plt.plot(xs, scan_list[2])
        # plt.plot(xs, dNLL)
        # plt.savefig("plots/NLL_"+entry.pdf.GetName() + str(j) + ".pdf")
        # plt.close(fig)
        #######################################

        # Method2 #############################
        dNLL, output = scanFit(profile, dscb_model, hist_toy, list[2], scan_size)
        xs = [x for x in range(len(dNLL))]
        fig = plt.figure()
        # plt.plot(xs, dNLL)
        plt.plot(xs, output[0])
        plt.plot(xs, output[1])
        plt.plot(xs, output[2])
        plt.plot(xs, output[3])
        plt.savefig("plots/NLL_"+entry.pdf.GetName() + str(j) + ".pdf")
        plt.close(fig)
        #######################################

        left = []
        right = []
        r_error_ = 0
        for i in range(len(dNLL) - 1):
            if dNLL[i] > 0.5 and dNLL[i+1] < 0.5: left.append(i)
            if dNLL[i] < 0.5 and dNLL[i+1] > 0.5: right.append(i)
        if len(left) == 0 or len(right) == 0: 
            r_error_ = -1
        else: r_error_ = (right[len(right) - 1] - left[0])*abs(list[2]) * scan_size
        
        
        
        r_sig.append(list[2])
        best_list.append(list[4])
        best_error.append(list[3])
        r_error.append(r_error_)
        if r_error_ > 0: 
            pull.Fill(list[2]/r_error_)
            pull_list.append(list[2]/r_error_)
        else: bad += 1
        print("Finish toy ", j+1)

    pull.Fit("gaus")
    # pull.Draw("HIST")
    can.SaveAs("plots/Pull_"+entry.pdf.GetName() + "_200.pdf")


    print("r = ", sum(r_sig)/N_toy)
    print("best error = ", sum(best_error)/N_toy)
    print("r error = ", sum(r_error)/N_toy)
    print("bad toys = ", bad)

    # print("r = ", r_sig)
    # print("best error = ", best_error)
    # print("r error = ", r_error)
    # print("bad toys = ", bad)
    # print("best func = ", best_list)
    # print("pull = ", pull_list)
    print("Finish ", entry.pdf.GetName()," toy sample")






#print ("nll = ", bern2Bias.nll.getVal())
#print ("corr nll = ", bern2Bias.corrNLL)
#print ("norm = ", extmodel.getNorm())
#i1 = bern2_model.pdf.createIntegral(x, ROOT.RooFit.Range('left,right'))
#i2 = bern2_model.pdf.createIntegral(x, ROOT.RooFit.Range('full'))
#print ("full int = ", i2.getVal(), " side int = ", i1.getVal())

#plotClass(x, reader.data_hist_bin1, bern2_model.pdf, title = "Bern2", sideBand = True, fitRange = 'left,right')
#plotClass(x, reader.data_hist_bin1, bern3_model.pdf, title = "Bern3", sideBand = True, fitRange = 'left,right')
#plotClass(x, reader.data_hist_bin1, bern4_model.pdf, title = "Bern4", sideBand = True, fitRange = 'left,right')
