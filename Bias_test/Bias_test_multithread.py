import ROOT
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import argparse
import threading

sys.path.append(os.path.abspath("../Utilities/"))
sys.path.append(os.path.abspath("../CMS_plotter/"))
import CMS_lumi, tdrstyle
from bkg_functions_class import *
from sig_functions_class import *
from Xc_Minimizer import *
from plot_utility import *
from sample_reader import *
from bias_class import *
# ROOT.gInterpreter.AddIncludePath('../Utilities/HZGRooPdfs.h')
#ROOT.gSystem.Load('../Utilities/HZGRooPdfs_cxx.so')

ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)

parser = argparse.ArgumentParser(description = "Number of toy samples")
parser.add_argument("N_toy")
args = parser.parse_args()

# Debug Flag
DEBUG = False

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
reader = readDat(list, dir = "../../sample/")
#reader = readRoot(x, "~/beijing_sample/data.root")
#cutHist = reader.data_hist_bin1.reduce(ROOT.RooFit.CutRange('left,right'))
#print ("cut norm = ", cutHist.sumEntries())
N = reader.data_u2.sumEntries()
#print("stats = ",N)
x.setBins(260)
x.setRange('left', lowx, 120)
x.setRange('right', 130, lowx+65)
x.setRange('full', lowx, lowx+65)

# Signal model (pre-fit)
MH = ROOT.RooRealVar("MH","MH"       ,125, 120., 130.)
dscb_model = DSCB_Class(x, MH, "bin1", 1.78,50, 100, 50, 100,0.845, 2.36)
dscb_model.setConst(True)
MH.setConstant(True)
N_sig = 10

# Functions to test
bern2_model_seed = Bern2Class(x, mu_gauss, "bin1_seed", 10, 0.3, 10, 7., 105.)
bern3_model_seed = Bern3Class(x, mu_gauss, "bin1_seed", 10, 0.3, 10, 3., 106.)
bern4_model_seed = Bern4Class(x, mu_gauss, "bin1_seed", 10, 0.3, 10, 3., 106.)
exp1_model_seed = Exp1Class(x, mu_gauss, "bin1_seed")

bern2_model = Bern2Class(x, mu_gauss, "bin1", 10, 0.3, 10, 7., 105.)
bern3_model = Bern3Class(x, mu_gauss, "bin1", 10, 0.3, 10, 3., 106.)
bern4_model = Bern4Class(x, mu_gauss, "bin1", 10, 0.3, 10, 3., 106.)
bern5_model = Bern5Class(x, mu_gauss, "bin1", 10, 0.3, 10, 3., 106.)
pow1_model = Pow1Class(x, mu_gauss, "bin1")
exp1_model = Exp1Class(x, mu_gauss, "bin1")
exp2_model = Exp2Class(x, mu_gauss, "bin1")
lau1_model = Lau1Class(x, mu_gauss, "bin1")
lau2_model = Lau2Class(x, mu_gauss, "bin1")
modg_model = ModGausClass(x, "bin1", 105., 170.)
#extmodel = ROOT.RooExtendPdf("extmodel", "Extended model", bern2_model.pdf, N, 'full');
#extmodel = ROOT.RooAddPdf("extmodel","extmodel", ROOT.RooArgList(bern3_model.pdf, dummy_sig), ROOT.RooArgList(N, dummy_N))
#r = bern2_model.pdf.fitTo(reader.data_hist_untagged1_bkg,ROOT.RooFit.Save(True), ROOT.RooFit.PrintLevel(-1), ROOT.RooFit.SumW2Error(True))
#r.Print("v")
#ROOT.RooFit.Range('left,right'),
profile_seed = [bern2_model_seed]#, bern3_model_seed, bern4_model_seed]
profile = [bern2_model, bern3_model, pow1_model, exp1_model,  modg_model]

# Set best-fit values
for entry in profile_seed:
    BiasClass(entry.pdf, reader.data_u2).minimize() #True, 'left,right'
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
            c2 = ROOT.RooRealVar("c2_"+ele.pdf.GetName(), "c2_"+ele.pdf.GetName(), 0., -50.*N_sig, 50.*N_sig)
        tot_model = ROOT.RooAddPdf("tot_model_"+ele.pdf.GetName(), "tot_model_"+ele.pdf.GetName(), ROOT.RooArgList(sig_model.pdf, ele.pdf), ROOT.RooArgList(c2, c1))
        bias = BiasClass(tot_model, hist, False, extend = True)
        bias.minimize()
        #result = tot_model.fitTo(hist, ROOT.RooFit.PrintLevel(-1),  ROOT.RooFit.Save(True), ROOT.RooFit.Minimizer("Minuit2"), ROOT.RooFit.Strategy(0))
        #result.Print("v")
        if i ==0: 
            ind = 0
            min_nll= bias.corrNLL #result.minNll() + .5*result.floatParsFinal().getSize()
            r_sig_ = c2.getVal()
            r_error_ = c2.getError()
            best_= ele.pdf.GetName()
        elif bias.corrNLL< min_nll: 
            ind = i
            min_nll = bias.corrNLL #result.minNll() + .5*result.floatParsFinal().getSize()
            r_sig_ = c2.getVal()
            r_error_ = c2.getError()
            best_= ele.pdf.GetName()
    return [ind, min_nll, r_sig_, r_error_, best_]

# Method 1 (WRONG)
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

# Method 2
def scanFit(profile_, sig_model, hist, r_sig_, scan_size_ = 0.3):
    scan_list_ = []
    scan_all_ = []
    output_all_ = []
    # Left r <= 0
    step = 0
    offset_nll = 0.
    scan = True
    if DEBUG: print ('Start scanning, N = ', N)
    while scan:
        choose = []
        if DEBUG: print ('Left step ', step)
        scan_sig =  abs(r_sig_) * step * scan_size_
        for pdf_ in profile_:
            #pdf_.reset()
            c1 = ROOT.RooRealVar("c1_"+ pdf_.pdf.GetName(), "c1_"+ pdf_.pdf.GetName(), N+scan_sig, 0., 5.*N)
            c2 = ROOT.RooRealVar("c2_"+ pdf_.pdf.GetName(), "c2_"+ pdf_.pdf.GetName(), -scan_sig )
            tot_model = ROOT.RooAddPdf("tot_"+ pdf_.pdf.GetName(), "tot_"+ pdf_.pdf.GetName(), ROOT.RooArgList(sig_model.pdf, pdf_.pdf), ROOT.RooArgList(c2, c1))
            #tot_model = ROOT.RooRealSumPdf("tot_"+ pdf_.pdf.GetName(), "tot_"+ pdf_.pdf.GetName(), ROOT.RooArgList(sig_model.pdf, pdf_.pdf), ROOT.RooArgList(c2, c1), False)
            bias = BiasClass(tot_model, hist, False, extend = True)
            if step==0: 
                bias.minimize(skip_hesse = True)
            else: bias.minimize(skip_hesse = True)
            choose.append(bias.corrNLL)
            if DEBUG:
                xframe = x.frame(ROOT.RooFit.Title("Left step " + str(step) + " " + pdf_.pdf.GetName()))
                hist.plotOn(xframe)
                tot_model.plotOn(xframe)
                can2 = ROOT.TCanvas("can2", "can2", 500, 500)
                can2.cd()
                xframe.Draw()
                can2.SaveAs("plots/Left_step" + str(step) + "_" + pdf_.pdf.GetName() + ".pdf")
        chose = min(choose)
        if step==0: offset_nll = chose
        scan_list_.insert(0, chose)
        scan_all_.insert(0, choose)
        step += 1
        if chose - offset_nll > 1.: scan = False
        elif step > 30: scan = False
    n_left_ = step - 1
    # Right r > 0
    step = 1
    scan = True
    while scan:
        choose = []
        if DEBUG: print ('Right step ', step)
        pdf_.reset()
        scan_sig =  abs(r_sig_) * step * scan_size_
        for pdf_ in profile_:
            if step == 1: pdf_.reset()
            #pdf_.reset()
            c1 = ROOT.RooRealVar("c1_"+ pdf_.pdf.GetName(), "c1_"+ pdf_.pdf.GetName(), N-scan_sig, 0, 5.*N)
            c2 = ROOT.RooRealVar("c2_"+ pdf_.pdf.GetName(), "c2_"+ pdf_.pdf.GetName(), scan_sig )
            tot_model = ROOT.RooAddPdf("tot_"+ pdf_.pdf.GetName(), "tot_"+ pdf_.pdf.GetName(), ROOT.RooArgList(sig_model.pdf, pdf_.pdf), ROOT.RooArgList(c2, c1))
            #tot_model = ROOT.RooRealSumPdf("tot_"+ pdf_.pdf.GetName(), "tot_"+ pdf_.pdf.GetName(), ROOT.RooArgList(sig_model.pdf, pdf_.pdf), ROOT.RooArgList(c2, c1), False)
            bias = BiasClass(tot_model, hist, False, extend = True)
            bias.minimize(skip_hesse = True)
            choose.append(bias.corrNLL)
            if DEBUG:
                xframe = x.frame(ROOT.RooFit.Title("Right step " + str(step) + " " + pdf_.pdf.GetName()))
                hist.plotOn(xframe)
                tot_model.plotOn(xframe)
                can2 = ROOT.TCanvas("can2", "can2", 500, 500)
                can2.cd()
                xframe.Draw()
                can2.SaveAs("plots/Right_step" + str(step) + "_" + pdf_.pdf.GetName() + ".pdf")
        chose = min(choose)
        scan_list_.append(chose)
        scan_all_.append(choose)
        step += 1
        if chose - offset_nll > 1.: scan = False
        elif step > 30: scan = False
    
    min_nll_ = min(scan_list_)
    for i in range(len(profile_)):
        output_all_.append([inx[i] - min_nll_ for inx in scan_all_])
    
    dNLL_ = [inx - min_nll_ for inx in scan_list_ ]
    return dNLL_, output_all_, n_left_

# Discrete profiling - Find minimum and (r_down, r_up)
# Scan at steps of signal_yield * scan_size around 0
# Stps when N_step >= 30 or deltaNLL > 1

def singleThread(thread_name, r_sig_, best_error_, r_error_list, bad_, pull_hist):
    scan_size = 0.2
    for entry in profile_seed:
        for j in range(int(args.N_toy)):   
            x.setBins(260)
            hist_toy = entry.pdf.generateBinned(x, ROOT.RooFit.NumEvents(N), ROOT.RooFit.Extended(True))
            list = profileFit(profile, dscb_model, hist_toy)

            # Method2 #############################
            dNLL, output, n_left = scanFit(profile, dscb_model, hist_toy, list[2], scan_size)
            if DEBUG:
                if j%100 == 0:
                    xs = [(inx - n_left) * abs(list[2]) * scan_size for inx in range(len(dNLL))]
                    fig = plt.figure()
                    for inx in range(len(profile)): plt.plot(xs, output[inx], marker = 'o')
                    plt.legend([tit.pdf.GetName() for tit in profile])
                    plt.plot(xs, dNLL, '--k')
                    plt.savefig("plots/NLL_"+entry.pdf.GetName() + "_" + str(j) + ".pdf")
                    plt.close(fig)
            #######################################

    # Find r_up and r_down
            left = []
            right = []
            r_error_ = 0
            for i in range(len(dNLL) - 1):
                if dNLL[i] > 0.5 and dNLL[i+1] < 0.5: left.append(i)
                if dNLL[i] < 0.5 and dNLL[i+1] > 0.5: right.append(i)
            if len(left) == 0 or len(right) == 0: 
                r_error_ = -1
    # r_error is approximately (r_up - r_down)/2
            else: r_error_ = (right[len(right) - 1] - left[0])*abs(list[2]) * scan_size /2
                    
            r_sig_.append(list[2])
            #best_list.append(list[4])
            best_error_.append(list[3])
            r_error_list.append(r_error_)
            if r_error_ > 0: 
                pull_hist.Fill(list[2]/r_error_)
                #pull_list.append(list[2]/r_error_)
            else: bad_ += 1
            print("Finish toy ", j+1, " Thread", thread_name)

if __name__ == "__main__":
    r_sig = []
    r_error = []
    best_error = []
    bad = 0
    can = ROOT.TCanvas("can", "can", 500, 500)
    pull = ROOT.TH1F("pull", "pull", 80, -4, 4)

    t1 = threading.Thread(target=singleThread, args=('1', r_sig, best_error, r_error, bad, pull))
    t2 = threading.Thread(target=singleThread, args=('2', r_sig, best_error, r_error, bad, pull))
    t3 = threading.Thread(target=singleThread, args=('3', r_sig, best_error, r_error, bad, pull))
    t4 = threading.Thread(target=singleThread, args=('4', r_sig, best_error, r_error, bad, pull))
    t5 = threading.Thread(target=singleThread, args=('5', r_sig, best_error, r_error, bad, pull))
    t1.start()
    t2.start()
    t3.start()
    t4.start()
    t5.start()

    t1.join()
    t2.join()
    t3.join()
    t4.join()
    t5.join()

    can.cd()
    pull.Fit("gaus")
    pull.GetXaxis().SetTitle("Pull")
    pull.SetTitle(entry.pdf.GetName() + " Pull")
    ROOT.gStyle.SetOptFit(1)
    can.SaveAs("plots/Pull_"+entry.pdf.GetName() + "_1000.pdf")


    print("r = ", sum(r_sig)/int(args.N_toy))
    print("best error = ", sum(best_error)/int(args.N_toy))
    print("r error = ", sum(r_error)/int(args.N_toy))
    print("bad toys = ", bad)

    print("Finish ", entry.pdf.GetName()," toy sample")



