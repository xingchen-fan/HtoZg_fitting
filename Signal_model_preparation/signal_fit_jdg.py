#!/usr/bin/env python3
import os
import sys
import ROOT
import argparse
import json
import logging
sys.path.append(os.path.abspath("./Utilities/"))
sys.path.append(os.path.abspath("./CMS_plotter/"))
sys.path.append(os.path.abspath("./Signal_model_preparation/"))
import CMS_lumi, tdrstyle
from sig_functions_class import *
from Xc_Minimizer import *
from plot_utility import *
from sample_reader import *
from config_parser_jdg import *


def performSignalFits():
#------------------Input Parsing--------------------#
#Inputs are taken from the input json file
    DISIGMA = True

    parser = argparse.ArgumentParser()
    parser.add_argument('-j', '--jsonIn', help = 'Signal fitting settings in the form of json file name')
    args = parser.parse_args()
    fname = args.jsonIn

    jfile = open(fname)
    configs = json.load(jfile)

    genSettings = configs["general_settings"]
    cats          = genSettings["categories"]
    flavs         = genSettings["flavors"]
    years         = genSettings["years"]
    prods         = genSettings["modes"]
    binFactor     = genSettings["bin_factor"]
    dispBinFactor = genSettings["display_bins"]
    rangeInfo     = genSettings["range_info"]
    sampParams = [(a,b,c,d) for a in cats for b in flavs for c in years for d in prods]

    signalSettings = configs["signal_settings"]
    ruiInBasePath     = signalSettings["base_path"]
    ggfPath           = signalSettings["ggf_path"]
    vbfPath           = signalSettings["vbf_path"]
    splits            = signalSettings["splits"] #0 split, 1 means combine. [categories, flavors, years, production modes]
    shapeConfigFile   = signalSettings["shape_config"]
    logFile           = signalSettings["log_file"]

    logging.basicConfig(filename=logFile, level=logging.INFO)

    bdtSettings = configs["bdt_settings"]
    ggfBins = bdtSettings["ggf_bins"]
    vbfBins = bdtSettings["vbf_bins"]

    plotPath = "./Signal_model_preparation/plots/"



    #------------------Combining Histograms--------------------#
    histList = []
    labelList = []

    for sampParam in sampParams:
        x = ROOT.RooRealVar("x", "mllg", 115, 130)
        x.setBins(int(binFactor*(x.getMax() - x.getMin())))
        if sampParam[3]=="ggf": sigSample = readRuiROOTSignal(x, ruiInBasePath+ggfPath, sampParam, ggfBins, rangeInfo)
        elif sampParam[3]=="vbf": sigSample = readRuiROOTSignal(x, ruiInBasePath+vbfPath, sampParam, vbfBins, rangeInfo)
        sigHist = sigSample.catflav
        histList.append(sigHist)
        labelList.append(sampParam)


    newHistList = []
    newLabelList = []
    while len(labelList)>1: #this might be some of the saddest code I have ever written. There are definitely better ways
        labelFirst = labelList[0]
        histFirst = histList[0]
        labelList.pop(0)
        histList.pop(0)
        myIter = 0
        newLabel = ""
        for i in range(len(labelFirst)):
            newLabel += ("_comb" if splits[i]==1 else ('_'+labelFirst[i]))
        while len(labelList)>0:
            match = True
            for i in range(len(labelList[0])):
                if (labelList[myIter][i] != labelFirst[i]) and splits[i]==0:
                    match = False

            if match == True:
                histFirst.add(histList[myIter])
                histList.pop(myIter)
                labelList.pop(myIter)
            else:
                myIter+=1
            if myIter >= len(labelList):
                break
        newHistList.append(histFirst)
        newLabelList.append(newLabel)

    #------------------Fit Histograms in Turn--------------------#
    for i in range(len(newHistList)):
        MH =  ROOT.RooRealVar("MH", "MH", 125, 120, 130)

        logging.info("Beginning fit %d"%i)
        sig_model = DSCB_Class(x, MH,  newLabelList[i], sigmaL_init = 1.5, sigmaR_init = 1.2, nL_init = 8, nL_bond = 18, nR_init = 12, nR_bond = 18, alphaL_init = 1.0, alphaR_init = 1.5, di_sigma = DISIGMA)
        sig_model.pdf.fitTo(newHistList[i], ROOT.RooFit.AsymptoticError(True), ROOT.RooFit.Save(True), ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0))
        sig_model.setStable()
        res = sig_model.pdf.fitTo(newHistList[i], ROOT.RooFit.AsymptoticError(True), ROOT.RooFit.Save(True),ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Strategy(0))
        res.Print('v')
        logging.info('nexp_' + newLabelList[i] + ' = %.2f'%newHistList[i].sumEntries())
        if DISIGMA:
            note = '#splitline{#splitline{#sigmaL = ' + '%.2f'%sig_model.sigmaL.getVal() + '#pm%.2f'%sig_model.sigmaL.getError()+', #sigmaR = ' + '%.2f'%sig_model.sigmaR.getVal() + '#pm%.2f'%sig_model.sigmaR.getError() + '}{nL = %.2f'%sig_model.nL.getVal()+', nR = %.2f'%sig_model.nR.getVal()+', #mu = %.2f'%MH.getVal() +  '}}{#alphaL = %.2f'%sig_model.alphaL.getVal() + ', #alphaR = %.2f'%sig_model.alphaR.getVal()+'}'
        else:
            note = '#splitline{#splitline{#sigma = ' + '%.2f'%sig_model.sigmaL.getVal() + '#pm%.2f'%sig_model.sigmaL.getError()+'}{nL = %.2f'%sig_model.nL.getVal()+', nR = %.2f'%sig_model.nR.getVal()+', #mu = %.2f'%MH.getVal() +  '}}{#alphaL = %.2f'%sig_model.alphaL.getVal() + ', #alphaR = %.2f'%sig_model.alphaR.getVal()+'}'
        plotClass(x, newHistList[i], sig_model.pdf, sig_model.pdf, newLabelList[i], plotPath, note=note, CMS = "Simulation", fullRatio = True, leftSpace=True, bins = dispBinFactor)

        sampleDetails = newLabelList[i].split('_')
        sampleDetails.pop(0)
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


if __name__=="__main__":
    print("Beginning signal fitting. . .")
    performSignalFits()
    print("Finished signal fitting. . .")
