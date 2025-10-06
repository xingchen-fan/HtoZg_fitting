#!/usr/bin/env python3
import os
import sys
import ROOT
import argparse
import json
import logging
import numpy as np
sys.path.append(os.path.abspath("./Utilities/"))
sys.path.append(os.path.abspath("./CMS_plotter/"))
sys.path.append(os.path.abspath("./Signal_model_preparation/"))
import CMS_lumi, tdrstyle
from sig_functions_class import *
from Xc_Minimizer import *
from plot_utility import *

def bootstrapper(modLabel, DISIGMA, binFactor, nEvents, shapeConfigFile):
    iterations = 50
    ROOT.RooRandom.randomGenerator().SetSeed(0)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)

    x = ROOT.RooRealVar("x", "mllg", 115, 130)
    x.setBins(int(binFactor*(x.getMax() - x.getMin())))
    MH =  ROOT.RooRealVar("MH", "MH", 125)
    sigmaEffList = []
    plotPath = "./Signal_model_preparation/plots/"

    initModel = DSCB_Class(x, MH,  modLabel, sigmaL_init = 1.5, sigmaR_init = 1.2, nL_init = 8, nL_bond = 18, nR_init = 12, nR_bond = 18, alphaL_init = 1.0, alphaR_init = 1.5, di_sigma = DISIGMA)
    initModel.assignValN(shapeConfigFile, modLabel)
    initModelPdf = initModel.pdf
    sigma_eff_pair = getEffSigma(x, initModelPdf)
    trueSigmaEff = (sigma_eff_pair[1]-sigma_eff_pair[0])/2.0

    i=0
    while i < iterations:

        sampleTH = initModelPdf.generateBinned(x, nEvents, False).createHistogram("sample_hist%d"%i, x, ROOT.RooFit.Binning(int(binFactor * (x.getMax() - x.getMin()))))
        sampleHist = ROOT.RooDataHist("sample_hist%d"%i, "sample_hist%d"%i, x, sampleTH)
        fit_model = DSCB_Class(x, MH,  "fit_model%d"%i, sigmaL_init = 1.5, sigmaR_init = 1.2, nL_init = 8, nL_bond = 18, nR_init = 12, nR_bond = 18, alphaL_init = 1.0, alphaR_init = 1.5, di_sigma = DISIGMA)
        fit_model.pdf.fitTo(sampleHist, ROOT.RooFit.AsymptoticError(True), ROOT.RooFit.Save(True), ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(2))
        res = fit_model.pdf.fitTo(sampleHist, ROOT.RooFit.AsymptoticError(True), ROOT.RooFit.Save(True),ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(2))

        if(res.status() == 0 or res.status() == 1):
            i+=1
        sigma_eff_pair = getEffSigma(x, fit_model.pdf)
        sigma_eff = (sigma_eff_pair[1]-sigma_eff_pair[0])/2.0
        sigmaEffList.append(sigma_eff)
        #plotClass(x, sampleHist, fit_model.pdf, fit_model.pdf , "disp_%d"%i, plotPath, note="test", CMS = "Simulation", fullRatio = True, leftSpace=True, bins = 2)

    return np.std(sigmaEffList)


if __name__=="__main__":
    DISIGMA = False
    modLabel = "Htozg_ggf1_el_2016APV_ggf_nominal"
    binFactor = 4
    nEvents = 130
    shapeConfigFile = "/afs/cern.ch/work/j/jgrassi/HtoZg_fitting/Config/config_DSCB_flat_jdg.json"
    print('%.2f'%bootstrapper(modLabel, DISIGMA, binFactor, nEvents, shapeConfigFile))

