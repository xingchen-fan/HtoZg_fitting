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

#ROOT.gInterpreter.AddIncludePath('../Utilities/HZGRooPdfs.h')
#ROOT.gSystem.Load('../Utilities/HZGRooPdfs_cxx.so')
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)

##############################################
# This script prepares the workspaces as input of combine. 
# Combine is needed and custom classes should be installed in combine instead!!!
# Each function in the profile needs to be fitted before entering the workspace.
# Signal model is fitted and fixed.
##############################################

#CATEGORIES = ["ggf4", "ggf3", "ggf2", "ggf1", "vbf4", "vbf3", "vbf2", "vbf1",
#              "vh3l", "vhmet", "tthhad", "tthlep"]

#other categories missing rn
CATEGORIES = ["ggf4", "ggf3", "ggf2", "ggf1"]

#can split according to year, prod mode, whatever, currently just el vs mu
SIGNAL_PROCS = ["Htozg_el", "Htozg_mu"] 

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
    workspace_filename = datacard_filename[:-4] + ".root"
    rawdata_filename = datacard_filename[:-4] + "_rawdata.root"
    output_file = ROOT.TFile(workspace_filename, "RECREATE")
    x = ROOT.RooRealVar("CMS_hzg_mass", "CMS_hzg_mass", 95.0, 180.0)
    x.setBins(85)
    reader = readPico(x, rawdata_filename)
    for category in CATEGORIES:
        write_workspace(category, reader, signal_config, background_config,
                        output_file, fake_data=True)
    output_file.Close()
    #append_nuisances_to_datacard(datacard_filename)

def write_workspace(category, reader, signal_config, background_config, 
                    output_file, fake_data=False):
    """Writes workspace for given category

    Args:
        category: name of category
        reader: pico reader
        signal_config: config file for signal
        background_config: config file for background
        output_file: root file to save workspaces to
        fake_data: using MC as data
    """

    output_file.cd()
    configs_bkg_ = {}
    with open(background_config, "r") as jfile_:
        configs_bkg_ = json.load(jfile_)
    lowx = configs_bkg_[category]["Range"]
    highx = lowx+65.0

    # Define variables
    # For now use full range, needs to match draw_pico
    mllg_name = f"mllg_cat_{category}"
    x = ROOT.RooRealVar(mllg_name, mllg_name, lowx, highx)
    x.setBinning(ROOT.RooUniformBinning(lowx, highx, 65),'full')
    x.setBins(65,'')
    x.setBins(1000,'cache')
    x.setRange('left', lowx, 120)
    x.setRange('right', 130, highx)
    mu_gauss = ROOT.RooRealVar("mu_gauss","always 0"       ,0.)
    MH = ROOT.RooRealVar("MH","MH"       ,125, 120., 130.)
    MH.setVal(125)
    MH.setConstant(True)

    #save workspace with data
    #have to now undo the hack we did to get all the histograms to have the
    #same variable 
    hist_name = f"data_obs_cat_{category}"
    data_hist_og = getattr(reader, f"hist_data_obs_cat_{category}")
    normalized_hist = ROOT.TH1D("","",65,lowx,highx)
    data_hist_og.fillHistogram(normalized_hist, 
        ROOT.RooArgList(reader.default_var))
    if fake_data:
        hist_to_toy(normalized_hist)
    hist_data = ROOT.RooDataHist(hist_name, hist_name, x, normalized_hist)
    ws_name = f"WS_data_obs_cat_{category}"
    ws_data = ROOT.RooWorkspace(ws_name, ws_name)
    #getattr(ws_data,"import")(x)
    getattr(ws_data,"import")(hist_data)
    ws_data.Write()
    nevents_data = hist_data.sumEntries()

    #save workspaces with signal
    configs_ = {}
    with open(signal_config, "r") as jfile_:
        configs_ = json.load(jfile_)
  
    for proc in SIGNAL_PROCS:
        proc_name = f"{proc}_cat_{category}"
        proc_settings = configs_[proc_name]
        signal_model = DSCB_Class(x, MH, proc_name, simple_name=True)
        signal_model.assignValModular(proc_settings)
        ws_signal = ROOT.RooWorkspace("WS_"+proc_name, "WS_"+proc_name)
        getattr(ws_signal,"import")(signal_model.pdf)
        ws_signal.Write()

    #save workspace with background
    profile_ = profileClass(x, mu_gauss, category, background_config)
    profile = profile_.testSelection("FT")
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
