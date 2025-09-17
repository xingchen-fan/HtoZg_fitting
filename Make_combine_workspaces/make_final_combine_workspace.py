#!/usr/bin/env python3
import math
import ROOT
import os
import json
import sys
import argparse
sys.path.append(os.path.abspath("../Utilities/"))
from bkg_functions_class import *
from sig_functions_class import *
from Xc_Minimizer import *
from plot_utility import *
from sample_reader import *
from profile_class import *
from sig_functions_class import *
from multiprocessing import Process

#ROOT.gInterpreter.AddIncludePath('../Utilities/HZGRooPdfs.h')
#ROOT.gSystem.Load('../Utilities/HZGRooPdfs_cxx.so')
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)

##############################################
# This script prepares the workspaces as input of combine. 
# Combine is needed and custom classes should be installed in combine instead!!!
# Each function in the profile needs to be fitted before entering the workspace.
# Signal model is fitted and fixed.
##############################################

CATEGORIES = ["ggf4", "ggf3", "ggf2", "ggf1", "vbf4", "vbf3", "vbf2", "vbf1"]

#can split according to year, prod mode, whatever, currently just el vs mu
SIGNAL_PROCS = ["Htozg_el", "Htozg_mu", "Htomm"] 

#stored as (name, is_scale)
SYSTEMATICS = []
'''
[("CMS_scale_e", True),
               ("CMS_res_e", False),
               ("CMS_scale_g", True),
               ("CMS_res_g", False),
               ("CMS_scale_m", True)]
'''
rng = ROOT.TRandom3()

def parse_args():
    """Parses arguments

    Returns:
      namespace with parsed arguments
    """
    parser = argparse.ArgumentParser(description = "Make workspace")
    parser.add_argument('-s', '--sig_config', help = 'Signal Configuration')
    parser.add_argument('-b', '--bak_config', help = 'Background Configuration')
    parser.add_argument('-d', '--datacard', help = 'Datacard filename')
    parser.add_argument('-t', '--threads', type=int, default=1)
    args = parser.parse_args()
    return args

def hist_to_toy(hist):
    """Turns a histogram into a set of toys generated from that histogram
    assuming each bin is the mean of a poisson random variable
  
    Args:
      hist: histogram to round (modified in place)
    """
    hist.SetBinContent(0, 0.0)
    hist.SetBinError(0, 0.0)
    hist.SetBinContent(hist.GetNbinsX()+1, 0.0)
    hist.SetBinError(hist.GetNbinsX()+1, 0.0)
    for ibin in range(0,hist.GetNbinsX()+2):
        bin_content = hist.GetBinContent(ibin)
        if bin_content < 0.0:
          bin_content = 0.0
        toy_value = rng.PoissonD(bin_content)
        hist.SetBinContent(ibin, toy_value)
        error = math.sqrt(toy_value)
        hist.SetBinError(ibin, math.sqrt(toy_value))

def round_bins(hist):
    """Rounds bins of a TH1D to the nearest integer, also removing negative 
    yields
  
    Args:
      hist: histogram to round (modified in place)
    """
    for ibin in range(hist.GetNbinsX()+2):
        rounded_value = float(round(hist.GetBinContent(ibin)))
        if (rounded_value < 0.0):
            rounded_value= 0.0
        hist.SetBinContent(ibin, rounded_value)

def write_all_workspaces(datacard_filename, signal_config, background_config):
    """Writes workspaces for each fit category

    Args:
        datacard_filename: name of datacard
        signal_config: name of json with signal configurations
        background_config: name of json with background configuration
    """
    #process each category
    category_processes = []
    for category in CATEGORIES:
        category_processes.append(Process(target=write_workspace, args=
            (category, datacard_filename, signal_config, background_config, False, False)))
        category_processes[-1].start()
    for icat in range(len(CATEGORIES)):
        category_processes[icat].join()
    #merge outputs
    #just deal with this in datacard for now?
    #workspace_filename = datacard_filename[:-4] + ".root"
    #output_file = ROOT.TFile(workspace_filename, "RECREATE")
    #category_filename = datacard_filename[:-4] + f"_{category}.root"
    #append_nuisances_to_datacard(datacard_filename)
    print("Finished processing all categories")

def prepare_sig(x, category, signal_config, datacard_filename):
    workspace_filename_sig = datacard_filename[:-4] + f"_{category}_sig.root"
    output_file_sig = ROOT.TFile(workspace_filename_sig, "RECREATE")
    output_file_sig.cd()
    configs_ = {}
    with open(signal_config, "r") as jfile_:
        configs_ = json.load(jfile_)
    MH = ROOT.RooRealVar("MH","MH"       ,125, 120., 130.)
    MH.setVal(125)
    MH.setConstant(True)
    shape_systematics = []
    for syst_name, is_scale in SYSTEMATICS:
        shape_systematics.append(ROOT.RooRealVar(syst_name, "", 0.0, -5.0,
                                                 5.0))
    for proc in SIGNAL_PROCS:
        proc_name = f"{proc}_cat_{category}"
        proc_settings = configs_[proc_name+"_nominal"]
        dMH = ROOT.RooRealVar(f"dMH_{proc_name}", "", 0.0, -2.0, 2.0)
        sigmaL = ROOT.RooRealVar(f"sigmaL_{proc_name}", "", 0.0, 0.01, 5.0)
        sigmaR = ROOT.RooRealVar(f"sigmaR_{proc_name}", "", 0.0, 0.01, 5.0)
        mu_arglist = ROOT.RooArgList(MH, dMH)
        sigmaL_arglist = ROOT.RooArgList(sigmaL)
        sigmaR_arglist = ROOT.RooArgList(sigmaR)
        mu_formulastr = "(@1"
        sigmaL_formulastr = "@0"
        sigmaR_formulastr = "@0"
        mu_idx = 2
        sigma_idx = 1
        for isyst in range(len(SYSTEMATICS)):
            syst_name = SYSTEMATICS[isyst][0]
            is_scale = SYSTEMATICS[isyst][1]
            syst_var = shape_systematics[isyst]
            alt_settings_up = configs_[proc_name+f"_{syst_name}Up"]
            alt_settings_down = configs_[proc_name+f"_{syst_name}Down"]
            if is_scale:
                variation = ((abs(proc_settings["dMH"]
                    -alt_settings_up["dMH"])+abs(proc_settings["dMH"]
                    -alt_settings_down["dMH"]))/2.0)/abs(proc_settings["dMH"])
                if (variation > 0.01):
                    mu_arglist.add(shape_systematics[isyst])
                    mu_formulastr += f"*(1.0+{variation}*@{mu_idx})"
                    mu_idx += 1
            else:
                variation1 = (((abs(proc_settings["sigmaL"]
                    -alt_settings_up["sigmaL"])+abs(proc_settings["sigmaL"]
                    -alt_settings_down["sigmaL"]))/2.0)
                    /abs(proc_settings["sigmaL"]))
                variation2 = (((abs(proc_settings["sigmaR"]
                    -alt_settings_up["sigmaR"])+abs(proc_settings["sigmaR"]
                    -alt_settings_down["sigmaR"]))/2.0)
                    /abs(proc_settings["sigmaR"]))
                if (variation1 > 0.01):
                    sigmaL_arglist.add(shape_systematics[isyst])
                    sigmaL_formulastr += (f"*(1.0+{variation1}*@{sigma_idx})")
                if (variation2 > 0.01):
                    sigmaR_arglist.add(shape_systematics[isyst])
                    sigmaR_formulastr += (f"*(1.0+{variation2}*@{sigma_idx})")

                sigma_idx += 1
        mu_formulastr += ")+@0"
        mu_final = ROOT.RooFormulaVar(f"mu_comb_{proc_name}","",mu_formulastr,
            mu_arglist)
        sigmaL_final = ROOT.RooFormulaVar(f"sigmaL_comb_{proc_name}","",
            sigmaL_formulastr,sigma_arglist)
        sigmaR_final = ROOT.RooFormulaVar(f"sigmaL_comb_{proc_name}","",
            sigmaR_formulastr,sigmaR_arglist)
        signal_model = DSCB_sys_Class(x, mu_final, sigmaL_final, sigmaR_final, dMH, sigmaL, sigmaR, proc_name, di_sigma=True)
        signal_model.assignValModular(proc_settings)
        ws_signal = ROOT.RooWorkspace("WS_"+proc_name, "WS_"+proc_name)
        getattr(ws_signal,"import")(signal_model.pdf)
        ws_signal.Write()
    output_file_sig.Close()
        
def prepare_bkg(x, category, background_config, datacard_filename, fake_data=False):
    mu_gauss = ROOT.RooRealVar("mu_gauss","always 0"       ,0.)
    rawdata_filename = datacard_filename[:-4] + "_rawdata.root"
    reader = readPico(x, rawdata_filename)
    #NOTE you MUST open output after initializing readPico
    workspace_filename = datacard_filename[:-4] + f"_{category}.root"
    output_file = ROOT.TFile(workspace_filename, "RECREATE")
    output_file.cd()

    #save workspace with data
    hist_name = f"data_obs_cat_{category}"
    hist_data = getattr(reader, hist_name)
    if fake_data:
      normalized_hist = ROOT.TH1D("","",65,lowx,highx)
      hist_data.fillHistogram(normalized_hist, x)
      hist_to_toy(normalized_hist)
      hist_data = ROOT.RooDataHist(hist_name, hist_name, x, normalized_hist)
    ws_name = f"WS_data_obs_cat_{category}"
    ws_data = ROOT.RooWorkspace(ws_name, ws_name)
    getattr(ws_data,"import")(hist_data)
    ws_data.Write()
    nevents_data = hist_data.sumEntries()

    #save workspace with bkg
    profile_ = profileClass(x, mu_gauss, category, background_config)
    profile = profile_.testSelection("FinalTemp")
    stat_list = []
    cuthist = hist_data.reduce(ROOT.RooFit.CutRange('left,right'))
    error = ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson)
    for model in profile:
        stat_list.append(ROOT.RooNLLVar("stat_"+model.pdf.GetName(),
            "stat_"+model.pdf.GetName(), model.SBpdf, cuthist))

    eps = 0.1
    strategy = 0
    stat_vals = []
    for stat in stat_list:
        Minimizer_NLL(stat, -1, 100, False, strategy)
        #Minimizer_Chi2(stat, -1, 1, False, strategy)
        r = Minimizer_NLL(stat, -1, eps, True, strategy)
        print(stat.GetName(), ": Cov q = ", r.covQual(), " status = ",
              r.status())
        stat_vals.append(stat.getVal() + 0.5*r.floatParsFinal().getSize())
    for model in profile:
        model.checkBond()

    best_ = stat_vals.index(min(stat_vals))

    cate = ROOT.RooCategory("pdfindex_"+category,
                            "Index of Pdf which is active for "+category)
    models = ROOT.RooArgList()
    for model in profile:
        models.add(model.pdf)
    bak_name = f"background_cat_{category}"
    pdf_bak_name = "pdf_"+bak_name
    multipdf = ROOT.RooMultiPdf(pdf_bak_name, pdf_bak_name, cate, models)
    # Penalty term
    #multipdf.setCorrectionFactor(0.)
    norm = ROOT.RooRealVar(pdf_bak_name +"_norm",
        "Number of background events", nevents_data, 0, 3*nevents_data)
    w_bkg = ROOT.RooWorkspace("WS_"+bak_name,"WS_"+bak_name)
    getattr(w_bkg, "import")(cate)
    getattr(w_bkg, "import")(norm)
    getattr(w_bkg, "import")(multipdf)
    w_bkg.Write()

    # Plot it
    multiPlotClass(x, hist_data, profile, title="Profile_" +category,
                   output_dir="", sideBand=True, fitRange= 'left,right',
                   best_index = best_)

    output_file.Close()


def write_workspace(category, datacard_filename, signal_config, 
                    background_config, fake_data=False, skip_bkg=False):
    """Writes workspace for given category

    Args:
        category: name of category
        datacard_filename: filename of datacard
        signal_config: config file for signal
        background_config: config file for background
        fake_data: using MC as data
    """
    configs_bkg_ = {}
    with open(background_config, "r") as jfile_:
        configs_bkg_ = json.load(jfile_)
    lowx = configs_bkg_[category]["Range"][0]
    highx = configs_bkg_[category]["Range"][1]

    # Define variables
    # For now use full range, needs to match draw_pico
    mllg_name = f"mllg_cat_{category}"
    x = ROOT.RooRealVar(mllg_name, mllg_name, lowx, highx)
    nbins = int(configs_bkg_[category]["Bins"])
    x.setBinning(ROOT.RooUniformBinning(lowx, highx, nbins),'full')
    #x.setBinning(ROOT.RooUniformBinning(lowx-32.5, highx+32.5, 1000),'cache')
    x.setBins(20000, "cache")
    x.setBins(nbins,'')
    x.setRange('left', lowx, 120)
    x.setRange('right', 130, highx)
    prepare_sig(x, category, signal_config, datacard_filename)
    if not skip_bkg: 
        prepare_bkg(x, category, background_config, datacard_filename, fake_data)
        
def append_nuisances_to_datacard(datacard_filename):
    """Appends nuisance parameters to datacard

    Args:
        datacard_filename: filename of datacard
    """
    with open(datacard_filename, "a") as datacard:
        for category in CATEGORIES:
            datacard.write(f"\npdfindex_{category} discrete")

def main():
    args = parse_args()
    write_all_workspaces(args.datacard, args.sig_config, args.bak_config)

if __name__=="__main__":
    main()
