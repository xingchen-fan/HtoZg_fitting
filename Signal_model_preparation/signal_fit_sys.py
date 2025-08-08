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
from multiprocessing import Process, Queue
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)

###############################################################################
# Based on Michael's script.
#
# Comment out the multithreading feature. May add back if possible.
#
# Signal models are not seperated by years and production mode in the nominal,
# so such features are removed.
#
# Both pico and Rui's signal root files can be used as input.
##############################################################################

def get_args():
    """Parses arguments

    Returns:
        Namespace with arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sampletype', help = 'Sample Type', default='')
    parser.add_argument('-i', '--input_file', help = 'Input File', default='')
    parser.add_argument('-c', '--cat', help = 'Category', default='')
    parser.add_argument('-v', '--variation', help = 'Variation', 
                        default='nominal')
    parser.add_argument('-con', '--config', help = 'Configuration', default='')
    return parser.parse_args()

def get_hist(x, CAT, lep_flavor):
    """Gets appropriate signal histogram from Rui's samples

    Args:
        CAT: category name
        prod: production mode
        lep_flavor: lepton flavor

    Returns:
       RooDataHist for appropriate category, production mode, and lepton flavor
    """
    if not CAT in ["ggf1","ggf2","ggf3","ggf4","vbf1","vbf2","vbf3","vbf4"]:
        raise ValueError("Unknown category")
    #if not prod in ["ggf","vbf"]:
    #    raise ValueError("Unknown production mode")
    if not lep_flavor in ["El","Mu"]:
        raise ValueError("Unknown lepton flavor")
    sig_sample = None
    if "ggf" in CAT:
        sig_sample = readRuiROOTggFSignal(x, '/eos/user/r/rzou/SWAN_projects/Classifier/Output_ggF_rui_commonparam/', 0.91,0.82,0.61)
    elif "vbf" in CAT:
        sig_sample = readRuiROOTVBFSignal(x, '/eos/user/r/rzou/SWAN_projects/Classifier/Output_VBF_rui_commonparam/', 0.95, 0.91,0.76)
    return getattr(sig_sample, f"{CAT}{lep_flavor}")

def get_param_dict(dscb_model, nexp):
    """Gets dictionary of fit information to save to json

    Args:
        dscb_model: DSCB_Class after fit is performed
        nexp: number of signal events for this process/category

    Returns:
        dictionary with fit information
    """
    fit_params = {}
    #N.B. sigmaL is two-sided sigma when disigma=False
    fit_params["disigma"] = dscb_model.disigma
    fit_params["nexp"] = nexp
    for key in ["dMH", "sigmaL", "sigmaR", "nL", "nR", "alphaL", "alphaR"]:
        fit_params[key] = getattr(dscb_model,key).getVal()
    return fit_params

def save_params_to_json(fit_params, filename, preset_name):
    """Saves DSCB parameter values to json, overwriting any existing values
    for the given preset_name

    Args:
        fit_params: dictionary with fit information
        filename: json file to save model to
        preset_name: key in file to access stored parameters
    """
    dscb_config = {}
    if os.path.isfile(filename):        
        with open(filename, "r") as config_file:
            dscb_config = json.load(config_file)
    dscb_config[preset_name] = fit_params
    with open(filename, "w") as config_file:
      config_file.write(json.dumps(dscb_config,indent=2))

def perform_signal_model_fit(x, MH, hist, name, asymm_gaussian = False):
    """Performs fit to signal with signal model (DSCB)

    Args:
        x: fitting RooRealVar
        MH: higgs mass variable
        hist: RooDataHist to fit
        name: name of signal model
        asymm_gaussian: whether or not to use symmetric sigma

    Returns: 
        dictionary of fit parameters
    """
    dMH = ROOT.RooRealVar("dMH_"+name, "", 0.0, -2.0, 2.0)
    mean_formula = ROOT.RooFormulaVar("mean_formula_"+name, "", "@0+@1", 
        ROOT.RooArgList(MH,dMH))
    sigmaL = ROOT.RooRealVar("sigmaL_"+name, "", 2, 0.01, 5.0)
    sigmaR = ROOT.RooRealVar("sigmaR_"+name, "", 2, 0.01, 5.0)
    sig_model = DSCB_sys_Class(x, mean_formula, sigmaL, sigmaR, dMH, sigmaL, sigmaR, name, di_sigma = False)
    """
    sig_model.pdf.fitTo(hist, ROOT.RooFit.AsymptoticError(True), 
        ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0))
    sig_model.pdf.fitTo(hist, ROOT.RooFit.AsymptoticError(True), 
        ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0))
    """
    if asymm_gaussian:
        sig_model = DSCB_sys_Class(x, mean_formula, sigmaL, sigmaR, dMH, sigmaL, sigmaR, name, alphaL_init = 1.5, alphaR_init = 1.5, nL_init = 4, nR_init = 8, di_sigma = True)
        #sig_model = DSCB_Class(x, MH, name,  sigmaL_init = 1.2, sigmaR_init = 1.2, nL_init = 4, nL_bond = 100, nR_init = 8, nR_bond = 100, alphaL_init = 1.5, alphaR_init = 1.5, di_sigma = True)
        #sig_model.pdf.fitTo(hist, ROOT.RooFit.AsymptoticError(True),
        #ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0))
        nll = ROOT.RooNLLVar("nll_"+ name, "nll_"+ name, sig_model.pdf, hist, ROOT.RooFit.AsymptoticError(True))
        Minimizer_NLL(nll, -1, 100, False, 0)
        Minimizer_NLL(nll, -1, 1, True, 0)
        Minimizer_NLL(nll, -1, 0.001, True, 0)
        #sig_model.setStable()
        res = Minimizer_NLL(nll, -1, 0.001, True, 0)
        #res = sig_model.pdf.fitTo(hist, ROOT.RooFit.AsymptoticError(True),ROOT.RooFit.Save(True),ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0))
        res.Print('v')
    nexp = hist.sumEntries()
    print(name+'_nexp = %.2f'%nexp)
    
    if asymm_gaussian:
        note = '#splitline{#splitline{#sigmaL = ' + '%.2f'%sig_model.sigmaL.getVal() + '#pm%.2f'%sig_model.sigmaL.getError()+', #sigmaR = ' + '%.2f'%sig_model.sigmaR.getVal() + '#pm%.2f'%sig_model.sigmaR.getError() + '}{nL = %.2f'%sig_model.nL.getVal()+', nR = %.2f'%sig_model.nR.getVal()+', #mu = %.2f'%(MH.getVal()+dMH.getVal()) +  '}}{#alphaL = %.2f'%sig_model.alphaL.getVal() + ', #alphaR = %.2f'%sig_model.alphaR.getVal()+'}'
    else:
        note = '#splitline{#splitline{#sigma = ' + '%.2f'%sig_model.sigmaL.getVal() + '#pm%.2f'%sig_model.sigmaL.getError()+'}{nL = %.2f'%sig_model.nL.getVal()+', nR = %.2f'%sig_model.nR.getVal()+', #mu = %.2f'%(MH.getVal()+dMH.getVal()) +  '}}{#alphaL = %.2f'%sig_model.alphaL.getVal() + ', #alphaR = %.2f'%sig_model.alphaR.getVal()+'}'
    plotClass(x, hist, sig_model.pdf, sig_model.pdf, name, note=note, 
        CMS = "Simulation", fullRatio = True, leftSpace=True, bins = 2)
    return get_param_dict(sig_model, nexp)

def do_category_fit(args, queue=None):
    """Performs fit to signal in appropriate category and draws output plot.
    The stdout from this process can be used with config_parser.py in this
    directory to produce the DSCB config file

    Args:
      args: argument dictionary as returned by ArgumentParser or similar
      queue: queue to store multiprocessing output
    """
    if args.sampletype=="pico":
        x = ROOT.RooRealVar("x", "mllg", 118, 130)
        MH = ROOT.RooRealVar("MH", "MH", 125, 120, 130)
        MH.setConstant()
        x.setBins(int(4* (x.getMax() - x.getMin())))
        pico_reader = readPico(x, args.input_file, args.cat, args.prod, args.variation)
        name = f"{args.prod}_cat_{args.cat}_{args.variation}"
        pico_name = f"mcdata_{name}"
        hist = getattr(pico_reader, pico_name)
        output = perform_signal_model_fit(x, MH, hist, name, True)
        if queue != None:
          queue.put(output)
        else:
          return output
    else:
        CAT = args.cat
        x = ROOT.RooRealVar("x", "mllg", 115, 130)
        MH =  ROOT.RooRealVar("MH", "MH", 125, 120, 130)
        MH.setConstant()
        x.setBins(int(4* (x.getMax() - x.getMin())))
        hist = get_hist(x, CAT, args.prod)
        name = CAT+"_"+args.prod.lower()
        output = perform_signal_model_fit(x, MH, hist, name, asymm_gaussian=True)
        return output
    
def perform_all_category_fits(args):
    """Perform fits to all signal processes and categories 
    Kill the multithreading due to lxplus memory restriction
    Args:
      args: arguments from ArgumentParser; only input_file is needed
    """
    #written for pico, but easily modifiable for Rui samples
    args_proxy = types.SimpleNamespace()
    args_proxy.sampletype = args.sampletype
    args_proxy.input_file = args.input_file
    cats = ["ggf4","ggf3","ggf2","ggf1","vbf4","vbf3","vbf2","vbf1","vh3l",
            "vhmet","tthhad","tthlep"]
    if args.cat == 'ggf':
        cats = ["ggf4","ggf3","ggf2","ggf1"]
    elif args.cat == 'vbf':
        cats = ["vbf4","vbf3","vbf2","vbf1"]
    elif args.cat == 'vhtth':
        cats = ["vh3l","vhmet","tthhad","tthlep"]
    else:
        cats = [args.cat]
    cat_fit_processes = []
    cat_fit_names = []
    cat_fit_output = Queue()
    #start all fits with multiprocessing
    if args.sampletype == "pico":
        for cat in cats:
            for prod in ["Htozg_el","Htozg_mu","Htomm"]:
                for syst in ["nominal", "CMS_scale_eUp", "CMS_scale_eDown","CMS_res_eUp", "CMS_res_eDown", "CMS_scale_gUp", "CMS_scale_gDown", "CMS_res_gUp", "CMS_res_gDown","CMS_scale_mUp", "CMS_scale_mDown"]:
                    args_proxy.cat = cat
                    args_proxy.prod = prod
                    args_proxy.variation = syst
                    #cat_fit_processes.append(Process(target=do_category_fit, 
                    #    args=(args_proxy,cat_fit_output,)))
                    #cat_fit_processes[-1].start()
                    print('waiting for '+f"{prod}_cat_{cat}_{syst}")
                    fit_params = do_category_fit(args_proxy)
                    save_params_to_json(fit_params, args.config, f"{prod}_cat_{cat}_{syst}")
                    #cat_fit_names.append(f"{prod}_cat_{cat}_{syst}")
    else:
        for cat in cats:
            for lep in ["El", "Mu"]:
                args_proxy.cat = cat
                args_proxy.prod = lep
                lep_lo = lep.lower()
                print('waiting for '+f"cat_{cat}_{lep_lo}_nominal")
                fit_params = do_category_fit(args_proxy)
                save_params_to_json(fit_params, args.config, f"Htozg_{lep_lo}_cat_{cat}_nominal")
    #save outputs
    """
    for iproc in range(len(cat_fit_names)):
        print('waiting for '+cat_fit_names[iproc])
        cat_fit_processes[iproc].join()
        #fit_params = do_category_fit()
        fit_params = cat_fit_output.get()
        save_params_to_json(fit_params, "../Config/DSCB_config_new.json",
                            cat_fit_names[iproc])
   """     
    print("Successfully finished all fits")
    
if __name__=="__main__":
    args = get_args()
    #previous functionality would be equivalent to just calling do_category_fit
    #do_category_fit(args)
    perform_all_category_fits(args)
