#!/usr/bin/env python3
import os
import sys
import ROOT
import argparse
import json
import types
sys.path.append(os.path.abspath("../Utilities/"))
sys.path.append(os.path.abspath("../CMS_plotter/"))
import CMS_lumi, tdrstyle
from sig_functions_class import *
from Xc_Minimizer import *
from plot_utility import *
from sample_reader import *

ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sampletype', help = 'Sample Type', default='')
    parser.add_argument('-i', '--input_file', help = 'Input File', default='')
    parser.add_argument('-c', '--cat', help = 'Category', default='')
    parser.add_argument('-y', '--year', help = 'Year', default='')
    parser.add_argument('-p', '--prod', help = 'Production mode', default='')
    return parser.parse_args()

def get_hist(CAT, prod, lep_flavor):
    """Gets appropriate signal histogram from Rui's samples

    Args:
        CAT: category name
        prod: production mode
        lep_flavor: lepton flavor
    """
    if not CAT in ["ggf1","ggf2","ggf3","ggf4","vbf1","vbf2","vbf3","vbf4"]:
        raise ValueError("Unknown category")
    if not prod in ["ggf","vbf"]:
        raise ValueError("Unknown production mode")
    if not lep_flavor in ["El","Mu"]:
        raise ValueError("Unknown lepton flavor")
    sig_sample = None
    if "ggf" in CAT and args.prod == "ggf":
        sig_sample = readRuiROOTggFSignalggF(x, '/eos/user/r/rzou/SWAN_projects/Classifier/Output_ggF_pinnacles_fix/relpt_peking_run2p3/TreeS/', args.year, 0.81,0.64,0.47)
    elif "ggf" in CAT and args.prod == "vbf":
        sig_sample = readRuiROOTggFSignalVBF(x, '/eos/user/r/rzou/SWAN_projects/Classifier/Output_ggF_pinnacles_fix/relpt_peking_run2p3/TreeS/', args.year, 0.81,0.64,0.47)
    elif "vbf" in CAT and args.prod == "ggf":
        sig_sample = readRuiROOTVBFSignalggF(x, '')
    elif "vbf" in CAT and args.prod == "vbf":
        sig_sample = readRuiROOTVBFSignalVBF(x, '')
    return getattr(sig_sample, f"{CAT}{lep_flavor}")

def save_params_to_json(dscb_model, MH, filename, preset_name):
    """Saves DSCB parameter values to json, overwriting any existing values
    for the given preset_name

    Args:
        dscb_model: DSCB_Class after fit is performed
        MH: RooRealVar of DSCB mu 
        filename: json file to save model to
        preset_name: key in file to access stored parameters
    """
    #generate dictionary with parameters
    fit_params = {}
    #N.B. sigmaL is two-sided sigma when disigma=False
    fit_params["disigma"] = dscb_model.disigma
    for key in ["sigmaL","sigmaR","nL","nR","alphaL","alphaR"]:
        fit_params[key] = getattr(dscb_model,key).getVal()
    fit_params["MH"] = MH.getVal()

    #save to file
    dscb_config = {}
    if os.path.isfile(filename):
        with open(filename, "r") as config_file:
            dscb_config = json.load(config_file)
    dscb_config[preset_name] = fit_params
    with open(filename, "w") as config_file:
      config_file.write(json.dumps(dscb_config,indent=2))

def perform_signal_model_fit(x, MH, hist, name, asymm_gaussian = True, 
                             save_json=''):
    """Performs fit to signal with signal model (DSCB)

    Args:
        x: fitting RooRealVar
        MH: higgs mass variable
        hist: RooDataHist to fit
        name: name of signal model
        asymm_gaussian: whether or not to use symmetric sigma
        save_json: filename to save json parameters to (empty string skips)
    """
    DISIGMA = asymm_gaussian
    sig_model = DSCB_Class(x, MH, name, di_sigma = DISIGMA)
    #print(hist.sumEntries())
    #nll = ROOT.RooNLLVar("nll_"+ name, "nll_"+ name, sig_model.pdf, hist, 
    #    ROOT.RooFit.AsymptoticError(True))
    
    #Minimizer_NLL(nll, -1, 100, False, 0)
    sig_model.pdf.fitTo(hist, ROOT.RooFit.AsymptoticError(True), 
        ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0))
    #sig_model.pdf.chi2FitTo(hist, 
    #    ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2), 
    #    ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0))
    sig_model.setStable()
    res = sig_model.pdf.fitTo(hist, ROOT.RooFit.AsymptoticError(True), 
        ROOT.RooFit.Save(True),ROOT.RooFit.PrintLevel(-1),
        ROOT.RooFit.Strategy(0))
    #res_el = sig_model.pdf.chi2FitTo(hist, 
    #    ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2), 
    #    ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0), 
    #    ROOT.RooFit.Save(True))
    res.Print('v')
    print(name+'_nexp = %.2f'%hist.sumEntries())
    
    if DISIGMA:
        note = '#splitline{#splitline{#sigmaL = ' + '%.2f'%sig_model.sigmaL.getVal() + '#pm%.2f'%sig_model.sigmaL.getError()+', #sigmaR = ' + '%.2f'%sig_model.sigmaR.getVal() + '#pm%.2f'%sig_model.sigmaR.getError() + '}{nL = %.2f'%sig_model.nL.getVal()+', nR = %.2f'%sig_model.nR.getVal()+', #mu = %.2f'%MH.getVal() +  '}}{#alphaL = %.2f'%sig_model.alphaL.getVal() + ', #alphaR = %.2f'%sig_model.alphaR.getVal()+'}'
    
    else:
        note = '#splitline{#splitline{#sigma = ' + '%.2f'%sig_model.sigmaL.getVal() + '#pm%.2f'%sig_model.sigmaL.getError()+'}{nL = %.2f'%sig_model.nL.getVal()+', nR = %.2f'%sig_model.nR.getVal()+', #mu = %.2f'%MH.getVal() +  '}}{#alphaL = %.2f'%sig_model.alphaL.getVal() + ', #alphaR = %.2f'%sig_model.alphaR.getVal()+'}'
    plotClass(x, hist, sig_model.pdf, sig_model.pdf, name, note=note, 
        CMS = "Simulation", fullRatio = True, leftSpace=True, bins = 2)
    if save_json != "":
      save_params_to_json(sig_model, MH, save_json, name)

def do_category_fit(args):
    """Performs fit to signal in appropriate category and draws output plot.
    The stdout from this process can be used with config_parser.py in this
    directory to produce the DSCB config file

    Args:
      argument dictionary as returned by ArgumentParser or similar
    """
    if args.sampletype=="pico":
        x = ROOT.RooRealVar("x", "mllg", 118, 130)
        MH =  ROOT.RooRealVar("MH", "MH", 125, 120, 130)
        x.setBins(int(4* (x.getMax() - x.getMin())))
        pico_reader = readPico(x, args.input_file)
        name = f"{args.prod}_cat_{args.cat}"
        hist = getattr(pico_reader, f"hist_{name}")
        perform_signal_model_fit(x, MH, hist, name, False,
            save_json="../Config/DSCB_config_new.json")
    else:
        CAT = args.cat
        for lepton_flavor in ["El","Mu"]:
            x = ROOT.RooRealVar("x", "mllg", 118, 130)
            MH =  ROOT.RooRealVar("MH", "MH", 125, 120, 130)
            x.setBins(int(4* (x.getMax() - x.getMin())))
            hist = get_hist(CAT, args.prod, lepton_flavor)
            name = CAT+"_"+args.year+"_"+lepton_flavor.lower()+"_"+args.prod
            perform_signal_model_fit(x, MH, hist, name)

def perform_all_category_fits(args):
    """Perform fits to all signal processes and categories 

    Args:
      args: arguments from ArgumentParser; only input_file is needed
    """
    #written for pico, but easily modifiable for Rui samples
    args_proxy = types.SimpleNamespace()
    args_proxy.sampletype = "pico"
    args_proxy.input_file = args.input_file
    for cat in ["ggf4","ggf3","ggf2","ggf1","vbf4","vbf3","vbf2","vbf1","vh3l",
                "vhmet","tthhad","tthlep"]:
        for prod in ["Htozg_el","Htozg_mu"]:
            args_proxy.cat = cat
            args_proxy.prod = prod
            do_category_fit(args_proxy)
    
if __name__=="__main__":
    args = get_args()
    #previous functionality would be equivalent to just calling do_category_fit
    #do_category_fit(args)
    perform_all_category_fits(args)
