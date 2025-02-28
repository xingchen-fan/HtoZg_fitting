#!/usr/bin/env python3
import ROOT
import os
import sys
import json
#import matplotlib.pyplot as plt
import numpy as np
import argparse
sys.path.append(os.path.abspath("../../Utilities/"))
sys.path.append(os.path.abspath("../../CMS_plotter/"))
import CMS_lumi, tdrstyle
from bkg_functions_class import *
from sig_functions_class import *
from Xc_Minimizer import *
from plot_utility import *
from sample_reader import *
from bias_class import *
# ROOT.gInterpreter.AddIncludePath('../Utilities/HZGRooPdfs.h')
import ctypes
#ROOT.gSystem.Load('../../Utilities/HZGRooPdfs_cxx.so')
ROOT.gInterpreter.Declare("""
RooDataSet readToy(TString filename, int& index){
auto file_ = new TFile(filename, "READ");
TDirectoryFile *dir = (TDirectoryFile *)file_->Get("toys");
RooDataSet *oneToy  = (RooDataSet *)dir->Get(("toy_" + std::to_string(index)).c_str());
return *oneToy;
}
""")

ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)

parser = argparse.ArgumentParser(description = "Number of toy samples")
parser.add_argument('-n', '--Ntoys', help = "N_toy")
parser.add_argument( '-i', '--it', help="It")
parser.add_argument('-f', '--func', help = "Func")
parser.add_argument('-c', '--cat', help = "Category")

args = parser.parse_args()
jfile = open('../../Config/config.json', 'r')
configs = json.load(jfile)
CAT = args.cat
setting = configs[CAT]
lowx = setting["Range"]

# Debug Flag
DEBUG = False
DEBUG_SCAN = False
# Define variables

x = ROOT.RooRealVar("CMS_hzg_mass", "CMS_hzg_mass", lowx, lowx + 65.)
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

# Read samples
ZBreader = readWsp(x, '/afs/cern.ch/user/f/fanx/EOS_space/zebing_sample/HZGamma_data_bkg_workspace_cat2.root', 'data_mass_cat0')
bkg_hist = ZBreader.hist
N = bkg_hist.sumEntries()

#print("stats = ",N)
x.setBins(260)
x.setRange('left', lowx, 120)
x.setRange('right', 130, lowx+65)
x.setRange('full', lowx, lowx+65)

# Signal model (pre-fit)
MH = ROOT.RooRealVar("MH","MH"       ,125, 120., 130.)
dscb_model = DSCB_Class(x, MH, CAT, 1.78,50, 100, 50, 100,0.845, 2.36)
dscb_model.setConst(True)
MH.setConstant(True)
N_sig = 10

# Conodr output file
#output = open('output_'+args.Tag+'.txt', 'w')

# Functions to test

bern2_model = Bern2Class(x, mu_gauss, CAT, setting["bern2"]['p0'], setting["bern2"]['p_init'],setting["bern2"]['bond'],setting["bern2"]['sigma_init'],setting["bern2"]['step_init'])
bern3_model = Bern3Class(x, mu_gauss, CAT, setting["bern3"]['p0'], setting["bern3"]['p_init'],setting["bern3"]['bond'],setting["bern3"]['sigma_init'],setting["bern3"]['step_init'])
bern4_model = Bern4Class(x, mu_gauss, CAT, setting["bern4"]['p0'], setting["bern4"]['p_init'],setting["bern4"]['bond'],setting["bern4"]['sigma_init'],setting["bern4"]['step_init'])
bern5_model = Bern5Class(x, mu_gauss, CAT, setting["bern5"]['p0'], setting["bern5"]['p_init'],setting["bern5"]['bond'],setting["bern5"]['sigma_init'],setting["bern5"]['step_init'])

pow1_model = Pow1Class(x, mu_gauss, CAT, setting["pow1"]['sigma_init'], setting["pow1"]['step_init'], setting["pow1"]['p_init'], setting["pow1"]['p_low'], setting["pow1"]['p_high'], setting["pow1"]['di_gauss'])
pow2_model = Pow2Class(x, mu_gauss, CAT, setting["pow2"]['sigma_init'], setting["pow2"]['step_init'], setting["pow2"]['p1_init'], setting["pow2"]['p1_low'], setting["pow2"]['p1_high'], setting["pow2"]['p2_init'], setting["pow2"]['p2_low'], setting["pow2"]['p2_high'], setting["pow2"]['f1_init'], setting["pow2"]['f2_init'], setting["pow2"]['xmax'], setting["pow2"]['const_f1'], setting["pow2"]['di_gauss'])
pow3_model = Pow3Class(x, mu_gauss, CAT,setting["pow3"]['sigma_init'], setting["pow3"]['step_init'], setting["pow3"]['p1_init'], setting["pow3"]['p1_low'], setting["pow3"]['p1_high'], setting["pow3"]['p2_init'], setting["pow3"]['p2_low'], setting["pow3"]['p2_high'], setting["pow3"]['p3_init'], setting["pow3"]['p3_low'], setting["pow3"]['p3_high'], setting["pow3"]['f1_init'], setting["pow3"]['f2_init'], setting["pow3"]['f3_init'], setting["pow3"]['xmax'], setting["pow3"]['const_f1'], setting["pow3"]['di_gauss'])

exp1_model = Exp1Class(x, mu_gauss, CAT, setting["exp1"]['sigma_init'], setting["exp1"]['step_init'], setting["exp1"]['p_init'], setting["exp1"]['p_low'], setting["exp1"]['p_high'], setting["exp1"]['di_gauss'])
exp2_model = Exp2Class(x, mu_gauss, CAT, setting["exp2"]['sigma_init'], setting["exp2"]['step_init'], setting["exp2"]['p1_init'], setting["exp2"]['p1_low'], setting["exp2"]['p1_high'], setting["exp2"]['p2_init'], setting["exp2"]['p2_low'], setting["exp2"]['p2_high'], setting["exp2"]['f1_init'], setting["exp2"]['f2_init'], setting["exp2"]['xmax'], setting["exp2"]['const_f1'], setting["exp2"]['di_gauss'])
exp3_model = Exp3Class(x, mu_gauss, CAT, setting["exp3"]['sigma_init'], setting["exp3"]['step_init'], setting["exp3"]['p1_init'], setting["exp3"]['p1_low'], setting["exp3"]['p1_high'], setting["exp3"]['p2_init'], setting["exp3"]['p2_low'], setting["exp3"]['p2_high'], setting["exp3"]['p3_init'], setting["exp3"]['p3_low'], setting["exp3"]['p3_high'], setting["exp3"]['f1_init'], setting["exp3"]['f2_init'], setting["exp3"]['f3_init'], setting["exp3"]['xmax'], setting["exp3"]['const_f1'], setting["exp3"]['di_gauss'])

lau2_model = Lau2Class(x, mu_gauss, CAT, setting["lau2"]['sigma_init'], setting["lau2"]['step_init'], setting["lau2"]['p1'], setting["lau2"]['p2'], setting["lau2"]['f_init'], setting["lau2"]['xmax'], setting["lau2"]['const_f1'], setting["lau2"]['di_gauss'])
lau3_model = Lau3Class(x, mu_gauss, CAT, setting["lau3"]['sigma_init'], setting["lau3"]['step_init'], setting["lau3"]['p1'], setting["lau3"]['p2'], setting["lau3"]['p3'], setting["lau3"]['f_init'], setting["lau3"]['xmax'], setting["lau3"]['const_f1'], setting["lau3"]['di_gauss'])
lau4_model = Lau4Class(x, mu_gauss, CAT, setting["lau4"]['sigma_init'], setting["lau4"]['step_init'], setting["lau4"]['p1'], setting["lau4"]['p2'], setting["lau4"]['p3'], setting["lau4"]['p4'], setting["lau4"]['f_init'], setting["lau4"]['xmax'], setting["lau4"]['const_f1'], setting["lau4"]['di_gauss'])

modg_model = ModGausClass(x, CAT, lowx, lowx+65, setting["modg"]['m0'], setting["modg"]['sl'], setting["modg"]['sh'], setting["modg"]['vl'], setting["modg"]['vr'])

profile = []
if "bern2" in setting["FT"]: profile.append(bern2_model)
if "bern3" in setting["FT"]: profile.append(bern3_model)
if "bern4" in setting["FT"]: profile.append(bern4_model)
if "bern5" in setting["FT"]: profile.append(bern5_model)
if "pow1" in setting["FT"]: profile.append(pow1_model)
if "pow2" in setting["FT"]: profile.append(pow2_model)
if "pow3" in setting["FT"]: profile.append(pow3_model)
if "exp1" in setting["FT"]: profile.append(exp1_model)
if "exp2" in setting["FT"]: profile.append(exp2_model)
if "exp3" in setting["FT"]: profile.append(exp3_model)
if "lau2" in setting["FT"]: profile.append(lau2_model)
if "lau3" in setting["FT"]: profile.append(lau3_model)
if "lau4" in setting["FT"]: profile.append(lau4_model)
if "modg" in setting["FT"]: profile.append(modg_model)

print("Done seed PDFs")

# Define profile Fit
def profileFit(profile_, sig_model, hist, fix = False, strength = 0.):
    min_nll = 0
    ind = 999
    r_sig_ = 0
    r_error_ = 0
    best_=''
    for i, ele in enumerate(profile_):
        #ele.reset()
        if (fix): 
            c1 = ROOT.RooRealVar("c1_"+ele.pdf.GetName(), "c1_"+ele.pdf.GetName(), N)
            c2 = ROOT.RooRealVar("c2_"+ele.pdf.GetName(), "c2_"+ele.pdf.GetName(), strength)
        else:
            c1 = ROOT.RooRealVar("c1_"+ele.pdf.GetName(), "c1_"+ele.pdf.GetName(), N, 0, 3.*N)
            c2 = ROOT.RooRealVar("c2_"+ele.pdf.GetName(), "c2_"+ele.pdf.GetName(), 0., -100.*N_sig, 100.*N_sig)
        tot_model = ROOT.RooAddPdf("tot_model_"+ele.pdf.GetName(), "tot_model_"+ele.pdf.GetName(), ROOT.RooArgList(sig_model.pdf, ele.pdf), ROOT.RooArgList(c2, c1))
        bias = BiasClass(tot_model, hist, False, "", extend = True)
        bias.minimize(debug = False)
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
def scanFit(profile_, sig_model, hist, r_scale_, scan_size_ = 0.3):
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
        scan_sig =  abs(r_scale_) * step * scan_size_
        for pdf_ in profile_:
            #pdf_.reset()
            c1 = ROOT.RooRealVar("c1_"+ pdf_.pdf.GetName(), "c1_"+ pdf_.pdf.GetName(), N+scan_sig, 0., 5.*N)
            c2 = ROOT.RooRealVar("c2_"+ pdf_.pdf.GetName(), "c2_"+ pdf_.pdf.GetName(), -scan_sig )
            tot_model = ROOT.RooAddPdf("tot_"+ pdf_.pdf.GetName(), "tot_"+ pdf_.pdf.GetName(), ROOT.RooArgList(sig_model.pdf, pdf_.pdf), ROOT.RooArgList(c2, c1))
            #tot_model = ROOT.RooRealSumPdf("tot_"+ pdf_.pdf.GetName(), "tot_"+ pdf_.pdf.GetName(), ROOT.RooArgList(sig_model.pdf, pdf_.pdf), ROOT.RooArgList(c2, c1), False)
            bias = BiasClass(tot_model, hist, False, "", extend = True)
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
        if chose - offset_nll > 0.75: scan = False
        elif step > 50: scan = False
    n_left_ = step - 1
    # Right r > 0
    step = 1
    scan = True
    while scan:
        choose = []
        if DEBUG: print ('Right step ', step)
        pdf_.reset()
        scan_sig =  abs(r_scale_) * step * scan_size_
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
        if chose - offset_nll > 0.75: scan = False
        elif step > 50: scan = False
    
    min_nll_ = min(scan_list_)
    for i in range(len(profile_)):
        output_all_.append([inx[i] - min_nll_ for inx in scan_all_])
    
    dNLL_ = [inx - min_nll_ for inx in scan_list_ ]
    return dNLL_, output_all_, n_left_

# Discrete profiling - Find minimum and (r_down, r_up)
# Scan at steps of signal_yield * scan_size around 0
# Stps when N_step >= 30 or deltaNLL > 1

scan_size = 0.1
generator = ROOT.TRandom()
r_sig = []
r_error = []
best_list = []
best_error = []
bad = 0
pull_list = []
    
for j in range(int(args.Ntoys)):
    toynum = j+1+ int(args.it)*5
    print (toynum)
    cppj = ctypes.c_int(int(toynum))
    scan_list = []      
    #x.setBins(260)
    #hist_toy = entry.pdf.generateBinned(x, ROOT.RooFit.NumEvents(generator.Poisson(N)))
    file_ = "../../Make_combine_workspaces/higgsCombine."+args.func+"."+ CAT+".GenerateOnly.mH125.123456.root"
    print("file = ", file_)
    hist_toy = ROOT.readToy(file_, cppj)
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
    dNLL, output, n_left = scanFit(profile, dscb_model, hist_toy, list[3], scan_size)
    
    if DEBUG_SCAN :
        xs = [(inx - n_left) * abs(list[3]) * scan_size for inx in range(len(dNLL))]
        #######################Matplot is not installed properly in some CMSSW#####################
        #fig = plt.figure()
        #for inx in range(len(profile)): plt.plot(xs, output[inx], marker = 'o')
        #plt.legend([tit.pdf.GetName() for tit in profile])
        #plt.plot(xs, dNLL, '--k')
        #plt.savefig("plots/NLL_"+entry.pdf.GetName() + "_" + str(j) + ".pdf")
        #plt.close(fig)
        can1 = ROOT.TCanvas('can1', 'can1', 600, 500)
        can1.cd()
        npx = np.array(xs, dtype=float)
        mg = ROOT.TMultiGraph()
        for inx in range(len(profile)):
            graph = ROOT.TGraph(len(xs), npx, np.array(output[inx], dtype=float))
            graph.SetMarkerStyle(20)
            graph.SetMarkerSize(0.5)
            graph.SetMarkerColor(inx + 2)
            graph.SetLineColor(inx+2)
            graph.SetTitle(profile[inx].pdf.GetName())
            mg.Add(graph)
        mg.Draw("ALP")
        can1.BuildLegend()
        can1.Draw()
        can1.SaveAs("plots/NLL_"+args.func + "_" + str(j) + ".pdf")
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
    else: 
        left1NLL = dNLL[left[0]]
        left2NLL = dNLL[left[0] + 1]
        acu_left = left[0] + (left1NLL - 0.5)/(left1NLL - left2NLL)
        right1NLL = dNLL[right[len(right) - 1]]
        right2NLL = dNLL[right[len(right) - 1] + 1]
        acu_right = right[len(right) - 1] + (right2NLL - 0.5)/(right2NLL - right1NLL)
        r_error_ = (acu_right - acu_left)*abs(list[3]) * scan_size /2
                
    r_sig.append(list[2])
    best_list.append(list[4])
    best_error.append(list[3])
    r_error.append(r_error_)
    if r_error_ > 0: 
        #pull.Fill(list[2]/r_error_)
        pull_list.append(list[2]/r_error_)
    else: bad += 1
    print("Finish toy ", j+1)

    
ncover = 0
for inx in range(len(r_sig)):
    if r_error[inx] > 0 and 0 > r_sig[inx] - r_error[inx] and 0 < r_sig[inx] + r_error[inx]:
        ncover += 1
print("r = ", sum(r_sig)/int(args.Ntoys))
print("best error = ", sum(best_error)/int(args.Ntoys))
print("r error = ", sum(r_error)/int(args.Ntoys))
print("bad toys = ", bad)
print("covered = ", ncover)
# print("r = ", r_sig)
# print("best error = ", best_error)
# print("r error = ", r_error)
# print("bad toys = ", bad)
print("best func = ", [best_list[i] for i in range(len(best_list)) if r_error[i] > 0])
print("pull = ", pull_list)
