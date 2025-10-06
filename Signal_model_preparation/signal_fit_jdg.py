#!/usr/bin/env python3
import os
import sys
import ROOT
import argparse
import json
import logging
import gc
sys.path.append(os.path.abspath("./Utilities/"))
sys.path.append(os.path.abspath("./CMS_plotter/"))
sys.path.append(os.path.abspath("./Signal_model_preparation/"))
import CMS_lumi, tdrstyle
from sig_functions_class import *
from Xc_Minimizer import *
from plot_utility import *
from sample_reader import *
from config_parser_jdg import *

def makeCombs(combs):
    genSettings = configs["general_settings"]
    cats          = genSettings["categories"]
    flavs         = genSettings["flavors"]
    years         = genSettings["years"]
    prods         = genSettings["modes"]
    if combs[0]==1:
        cats = ["comb"]
    if combs[1]==1:
        flavs = ["comb"]
    if combs[2]==1:
        years = ["comb"]
    if combs[3]==1:
        prods = ["comb"]
    merges = [(a,b,c,d) for a in cats for b in flavs for c in years for d in prods]
    return merges

def performSignalFits():
    #------------------Input Parsing--------------------#
    #Inputs are taken from the input json file
    DISIGMA = True
    genSettings = configs["general_settings"]
    binFactor     = genSettings["bin_factor"]
    dispBinFactor = genSettings["display_bins"]
    rangeInfo     = genSettings["range_info"]

    ruiInBasePath     = signalSettings["base_path"]
    ggfPath           = signalSettings["ggf_path"]
    vbfPath           = signalSettings["vbf_path"]
    fitCombs          = signalSettings["fit_combs"]
    dispCombs         = signalSettings["disp_combs"]
    shapeConfigFile   = signalSettings["shape_config"]

    bdtSettings = configs["bdt_settings"]
    ggfBins = bdtSettings["ggf_bins"]
    vbfBins = bdtSettings["vbf_bins"]

    fitParams = makeCombs(fitCombs)

    for i in range(len(dispCombs)):
        if fitCombs[i]>dispCombs[i]:
            exit("Error::If histograms are combined for fit, they must be combined for display as well.")

    dispParams = makeCombs(dispCombs)

    plotPath = "./Signal_model_preparation/plots/"
    
    #------------------Loading Histograms for fit--------------------#
    histList = []
    labelList = []

    for fitParam in fitParams:
        x = ROOT.RooRealVar("x", "mllg", 115, 130)
        x.setBins(int(binFactor*(x.getMax() - x.getMin())))
        sigSample = readRuiROOTSignal(x, ruiInBasePath, ggfPath, vbfPath, fitParam, ggfBins, vbfBins, rangeInfo)
        sigHist = sigSample.catflav
        histList.append(sigHist)
        label = '_'.join([s for s in fitParam])
        labelList.append(label)
    
    #------------------Fit Histograms in Turn--------------------#
    for i in range(len(histList)):
        MH =  ROOT.RooRealVar("MH", "MH", 125)
        logging.info("Beginning fit %d"%i)
        sig_model = DSCB_Class(x, MH,  labelList[i], sigmaL_init = 1.5, sigmaR_init = 1.2, nL_init = 8, nL_bond = 18, nR_init = 12, nR_bond = 18, alphaL_init = 1.0, alphaR_init = 1.5, di_sigma = DISIGMA)
        sig_model.pdf.fitTo(histList[i], ROOT.RooFit.AsymptoticError(True), ROOT.RooFit.Save(True), ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0))
        sig_model.setStable()
        res = sig_model.pdf.fitTo(histList[i], ROOT.RooFit.AsymptoticError(True), ROOT.RooFit.Save(True),ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0))
        res.Print('v')
        logging.info('nexp_' + labelList[i] + ' = %.2f'%histList[i].sumEntries())
        sigma_eff_pair = getEffSigma(x, sig_model.pdf)
        sigma_eff = (sigma_eff_pair[1]-sigma_eff_pair[0])/2.0
        if DISIGMA:
            note = '#splitline{#splitline{#splitline{#sigma_{eff} = ' + '%.2f'%sigma_eff + '}{#sigmaL = ' + '%.2f'%sig_model.sigmaL.getVal() + '#pm%.2f'%sig_model.sigmaL.getError()+', #sigmaR = ' + '%.2f'%sig_model.sigmaR.getVal() + '#pm%.2f'%sig_model.sigmaR.getError() + '}}{nL = %.2f'%sig_model.nL.getVal()+', nR = %.2f'%sig_model.nR.getVal()+', #mu = %.2f'%MH.getVal() +  '}}{#alphaL = %.2f'%sig_model.alphaL.getVal() + ', #alphaR = %.2f'%sig_model.alphaR.getVal()+'}'
        else:
            note = '#splitline{#splitline{#sigma = ' + '%.2f'%sig_model.sigmaL.getVal() + '#pm%.2f'%sig_model.sigmaL.getError()+'}{nL = %.2f'%sig_model.nL.getVal()+', nR = %.2f'%sig_model.nR.getVal()+', #mu = %.2f'%MH.getVal() +  '}}{#alphaL = %.2f'%sig_model.alphaL.getVal() + ', #alphaR = %.2f'%sig_model.alphaR.getVal()+'}'
        plotClass(x, histList[i], sig_model.pdf, sig_model.pdf, labelList[i], plotPath, note=note, CMS = "Simulation", fullRatio = True, leftSpace=True, bins = dispBinFactor)

        sampleDetails = labelList[i].split('_')
        logging.info("Ending fit %d"%i)
        logging.info("Parsing output to config json")
        #Storing the most recent output in a temp txt file. . .I am aware this is not a great way of doing this.
        fread = open(logFile, 'r')
        ftext = fread.read()
        beg = ftext.find("Beginning fit %d"%i)
        end = ftext.find("Ending fit %d"%i)
        fread.close()

        fwrite = open('test_config.txt', 'w')
        fwrite.write(ftext[beg:end])
        fwrite.write("Data copied")
        fwrite.close()

        configParser(shapeConfigFile, 'test_config.txt', sampleDetails[0], sampleDetails[1], sampleDetails[2], sampleDetails[3])
        logging.info("Done: parsed the log file info to the config json")
    
    #--------------------Loading histograms for display------------------#
    if combineDisp == False:
        exit()

    dispHistList = []
    dispLabelList = []

    for dispParam in dispParams:
        if combineDisp==False:
            break
        x = ROOT.RooRealVar("x", "mllg", 115, 130)
        x.setBins(int(binFactor*(x.getMax() - x.getMin())))
        sigSample = readRuiROOTSignal(x, ruiInBasePath, ggfPath, vbfPath, dispParam, ggfBins, vbfBins, rangeInfo)
        sigHist = sigSample.catflav
        dispHistList.append(sigHist)
        label = '_'.join([s for s in dispParam])
        dispLabelList.append(label)


    #--------------------Make various combined plots-----------#
    for i in range(len(dispHistList)):
        MH =  ROOT.RooRealVar("MH", "MH", 125)
        logging.info("Making display hist %d"%i)
        pdfList = []
        pdfDict = {}
        coefList = []
        coefDict = {}
        dispDetails = dispLabelList[i].split('_')
        
        for j in range(len(histList)):
            modLabel = "Htozg_" + labelList[j] + "_nominal"
            fitDetails = labelList[j].split('_')
            isSub = True
            for k in range(len(fitDetails)):
                if dispDetails[k]!=fitDetails[k] and dispDetails[k]!='comb':
                    isSub = False
            if isSub == True:
                pdfDict["model%d"%j] = DSCB_Class(x, MH,  modLabel, sigmaL_init = 1.5, sigmaR_init = 1.2, nL_init = 8, nL_bond = 18, nR_init = 12, nR_bond = 18, alphaL_init = 1.0, alphaR_init = 1.5, di_sigma = DISIGMA)
                pdfDict["model%d"%j].assignValN(shapeConfigFile, modLabel)
                pdfList.append(pdfDict["model%d"%j].pdf)
                coefDict["coef%d"%j] = ROOT.RooRealVar('c_'+modLabel,'c_'+modLabel, pdfDict["model%d"%j].nsig/dispHistList[i].sumEntries())
                if j >= 0:
                    coefList.append(coefDict["coef%d"%j])

        add_model = ROOT.RooAddPdf('multi_model_'+dispLabelList[i], 'multi_model_'+dispLabelList[i], ROOT.RooArgList(pdfList), ROOT.RooArgList(coefList))
        add_model.fixCoefNormalization(ROOT.RooArgSet(x)) #what does this do?
        dispNote = '#splitline{MC combined}{Signal model split}'
        plotClass(x, dispHistList[i], add_model, add_model, "disp_" + dispLabelList[i], plotPath, note=dispNote, CMS = "Simulation", fullRatio = True, leftSpace=True, bins = dispBinFactor)

        logging.info("Finished making display hist %d"%i)



if __name__=="__main__":
    print("Beginning signal fitting. . .")

    parser = argparse.ArgumentParser()
    parser.add_argument('-j', '--jsonIn', help = 'Fitting settings in the form of json file name')
    args = parser.parse_args()
    fname = args.jsonIn

    jfile = open(fname)
    configs = json.load(jfile)
    signalSettings = configs["signal_settings"]
    logFile = signalSettings["log_file"]

    clearing = open(logFile, 'w')
    clearing.close()

    logging.basicConfig(level=logging.INFO,stream=sys.stdout)
    ROOT.gSystem.RedirectOutput(logFile)

    combineDisp = True
    performSignalFits()
