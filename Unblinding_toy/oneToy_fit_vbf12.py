#!/usr/bin/env python3
import ROOT
import os
import sys
import json
#import matplotlib.pyplot as plt
import numpy as np
import argparse
sys.path.append(os.path.abspath("../Utilities/"))
sys.path.append(os.path.abspath("../CMS_plotter/"))
import CMS_lumi, tdrstyle
from bkg_functions_class import *
from sig_functions_class import *
from Xc_Minimizer import *
from plot_utility import *
from sample_reader import *
from bias_class import *
from profile_class import *
from sig_functions_class import *
import CMS_lumi, tdrstyle
# ROOT.gInterpreter.AddIncludePath('../Utilities/HZGRooPdfs.h')
import ctypes
#ROOT.gSystem.Load('../../Utilities/HZGRooPdfs_cxx.so')
ROOT.gInterpreter.Declare("""
RooDataSet readToy(TString filename, string num){
auto file_ = new TFile(filename, "READ");
TDirectoryFile *dir = (TDirectoryFile *)file_->Get("toys");
RooDataSet *oneToy  = (RooDataSet *)dir->Get(("toy_" + num).c_str());
return *oneToy;
}
""")

ROOT.gInterpreter.Declare("""
RooDataHist readToyWs(TString filename, string CAT, string lep){
auto file_ = new TFile(filename, "READ");
string tag = "";
if (lep != ""){tag = tag + "_" + lep;}
RooWorkspace *ws = (RooWorkspace *)file_->Get(("workspace_bkg"+tag).c_str());
RooDataHist *oneToy  = (RooDataHist *)ws->data(("hist_" + CAT +tag).c_str());
return *oneToy;
}
""")

ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)

parser = argparse.ArgumentParser(description = "S+B Asimov fit")
parser.add_argument('-con', '--config', help = "Configuration", default='../Config/chi2_config_xgboost_nodrop.json')
parser.add_argument('-l', '--lep', help = "Lepton", default='')

args = parser.parse_args()
jfile = open(args.config, 'r')
configs = json.load(jfile)
setting = configs['vbf2']
lowx = setting["Range"]

# Debug Flags
DEBUG = False
DEBUG_SCAN = True
PLOTFIT =False

# Define variables

x = ROOT.RooRealVar("CMS_hzg_mass_vbf12", "CMS_hzg_mass_vbf12", lowx, lowx + 65.)
mu_gauss = ROOT.RooRealVar("mu_gauss","always 0"       ,0.)

x.setBins(260)
x.setRange('left', lowx, 120)
x.setRange('right', 130, lowx+65)
x.setRange('full', lowx, lowx+65)

# Signal model
MH = ROOT.RooRealVar("MH","MH"       ,125, 120., 130.)
sig_model1_el = combineSignalLep(x, MH, 'vbf1', '../Config/config_DSCB.json', 'el')
sig_model1_mu = combineSignalLep(x, MH, 'vbf1', '../Config/config_DSCB.json', 'mu')
sig_model2_el = combineSignalLep(x, MH, 'vbf2', '../Config/config_DSCB.json', 'el')
sig_model2_mu = combineSignalLep(x, MH, 'vbf2', '../Config/config_DSCB.json', 'mu')
c_sig_el = ROOT.RooRealVar('c_sig_el', 'c_sig_el', sig_model1_el.ntot/(sig_model1_el.ntot+sig_model2_el.ntot))
c_sig_mu = ROOT.RooRealVar('c_sig_mu', 'c_sig_mu', sig_model1_mu.ntot/(sig_model1_mu.ntot+sig_model2_mu.ntot))
N_sig = 0

if args.lep == 'el':
    dscb_model = ROOT.RooAddPdf('sig_model12_el', 'sig_model12_el', ROOT.RooArgList(sig_model1_el.pdf, sig_model2_el.pdf), c_sig_el)
    N_sig = sig_model1_el.ntot+sig_model2_el.ntot
elif args.lep == 'mu':
    dscb_model = ROOT.RooAddPdf('sig_model12_mu', 'sig_model12_mu', ROOT.RooArgList(sig_model1_mu.pdf, sig_model2_mu.pdf), c_sig_mu)
    N_sig = sig_model1_mu.ntot+sig_model2_mu.ntot
MH.setConstant(True)
insig=0

# Functions to test

profile_class = profileClass(x, mu_gauss, 'vbf2', args.config)
profile = profile_class.testSelection("FT")
print ("N sig = ", N_sig)
print("Done seed PDFs")

# Define profile Fit
def profileFit(profile_, sig_model, data, fix = False, strength = 0.):
    min_nll = 99999999
    ind = 999
    r_sig_ = -999
    r_error_ = -999
    best_=''
    unstable = []
    hist = ROOT.RooDataHist("hist_hist", "hist_hist", ROOT.RooArgSet(x), data)
    for i, ele in enumerate(profile_):
        #ele.reset()
        if (fix): 
            c1 = ROOT.RooRealVar("c1_"+ele.pdf.GetName(), "c1_"+ele.pdf.GetName(), hist.sumEntries())
            c2 = ROOT.RooRealVar("c2_"+ele.pdf.GetName(), "c2_"+ele.pdf.GetName(), strength)
        else:
            c1 = ROOT.RooRealVar("c1_"+ele.pdf.GetName(), "c1_"+ele.pdf.GetName(), hist.sumEntries(), 0, 3.*hist.sumEntries())
            c2 = ROOT.RooRealVar("c2_"+ele.pdf.GetName(), "c2_"+ele.pdf.GetName(), 0., -100.*N_sig, 100.*N_sig)
        tot_model = ROOT.RooAddPdf("tot_model_"+ele.pdf.GetName(), "tot_model_"+ele.pdf.GetName(), ROOT.RooArgList(sig_model, ele.pdf), ROOT.RooArgList(c2, c1))
        bias = BiasClass(tot_model, hist, False, "", extend = True)
        bias.minimize(debug = False)
        #result = tot_model.fitTo(hist, ROOT.RooFit.PrintLevel(-1),  ROOT.RooFit.Save(True), ROOT.RooFit.Minimizer("Minuit2"), ROOT.RooFit.Strategy(0))
        #result.Print("v")
        if i ==0 and bias.stable == 0:
            ind = 0
            min_nll= bias.corrNLL #result.minNll() + .5*result.floatParsFinal().getSize()
            r_sig_ = c2.getVal()
            r_error_ = c2.getError()
            best_= ele.pdf.GetName()
        elif bias.corrNLL< min_nll and bias.stable == 0: 
            ind = i
            min_nll = bias.corrNLL #result.minNll() + .5*result.floatParsFinal().getSize()
            r_sig_ = c2.getVal()
            r_error_ = c2.getError()
            best_= ele.pdf.GetName()
        elif bias.stable != 0:
            unstable.append(i)
        if bias.stable == 0 and PLOTFIT:
            #hist_hist = ROOT.RooDataHist("hist_hist", "hist_hist", ROOT.RooArgSet(x), hist)
            plotClass(x, hist, tot_model, ele.pdf, title="Toy_hist_"+ele.pdf.GetName(), output_dir="", note="#splitline{NLL = %.2f"%bias.corrNLL + ", Raw NLL = %.2f"%bias.nll.getVal() + "}{Signal strength = %.2f"%(c2.getVal()/N_sig) + ", Error = %.2f"%(c2.getError()/N_sig) + "}", CMS = "Simulation", fullRatio = True, toy = True, bkg_normSF = c1.getVal()/(c1.getVal() + c2.getVal()))
            
    if len(unstable)!=0:
        for bad in unstable:
            del profile_[bad]
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
        c1 = ROOT.RooRealVar("c1_"+ bkgclass.pdf.GetName(), "c1_"+ bkgclass.pdf.GetName(), hist.sumEntries(), 0, 3.*hist.sumEntries())
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
    if DEBUG: print ('Start scanning, N = ', hist.sumEntries())
    while scan:
        choose = []
        if DEBUG: print ('Left step ', step)
        scan_sig =  abs(r_scale_) * step * scan_size_
        for pdf_ in profile_:
            #pdf_.reset()
            c1 = ROOT.RooRealVar("c1_"+ pdf_.pdf.GetName(), "c1_"+ pdf_.pdf.GetName(), hist.sumEntries()+scan_sig, 0., 5.*hist.sumEntries())
            c2 = ROOT.RooRealVar("c2_"+ pdf_.pdf.GetName(), "c2_"+ pdf_.pdf.GetName(), insig*N_sig-scan_sig )
            tot_model = ROOT.RooAddPdf("tot_"+ pdf_.pdf.GetName(), "tot_"+ pdf_.pdf.GetName(), ROOT.RooArgList(sig_model, pdf_.pdf), ROOT.RooArgList(c2, c1))
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
            c1 = ROOT.RooRealVar("c1_"+ pdf_.pdf.GetName(), "c1_"+ pdf_.pdf.GetName(), hist.sumEntries()-scan_sig, 0, 5.*hist.sumEntries())
            c2 = ROOT.RooRealVar("c2_"+ pdf_.pdf.GetName(), "c2_"+ pdf_.pdf.GetName(), insig*N_sig+scan_sig )
            tot_model = ROOT.RooAddPdf("tot_"+ pdf_.pdf.GetName(), "tot_"+ pdf_.pdf.GetName(), ROOT.RooArgList(sig_model, pdf_.pdf), ROOT.RooArgList(c2, c1))
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
    

scan_list = []      
#x.setBins(260)
#hist_toy = entry.pdf.generateBinned(x, ROOT.RooFit.NumEvents(generator.Poisson(N)))
file_ = 'workspaces/workspace_bkg_toy_unblind_vbf12.root'

print("file = ", file_)
hist_toy = ROOT.readToyWs(file_, 'vbf12', args.lep)
list = profileFit(profile, dscb_model, hist_toy)

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
    can1.SaveAs("NLL_vbf12"+ args.lep+".pdf")
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
    pull_list.append((list[2] - insig*N_sig)/r_error_)
else: bad += 1
print("Finish toy")

    
ncover = 0
for inx in range(len(r_sig)):
    if r_error[inx] > 0 and 0 > r_sig[inx] - r_error[inx] and 0 < r_sig[inx] + r_error[inx]:
        ncover += 1
print("bad toys = ", bad)
print("covered = ", ncover)
# print("r = ", r_sig)
# print("best error = ", best_error)
# print("r error = ", r_error)
# print("bad toys = ", bad)
print("best func = ", [best_list[i] for i in range(len(best_list)) if r_error[i] > 0])
print("all r = ", r_sig[0]/N_sig)
print("all r error from fit = ", best_error[0]/N_sig)
print("all r error from scan = ", r_error[0]/N_sig)
print("pull = ", pull_list)
