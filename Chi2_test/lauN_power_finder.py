#!/usr/bin/env python3
import time
import os
import sys
import ROOT
import json
import argparse
import logging
import itertools
sys.path.append(os.path.abspath("./Utilities/"))
sys.path.append(os.path.abspath("./CMS_plotter/"))
import CMS_lumi, tdrstyle
#from bkg_functions_fit import *
from bkg_functions_class import *
from Xc_Minimizer import *
from plot_utility import *
from sample_reader import *
#ROOT.gInterpreter.AddIncludePath('../Utilities/HZGRooPdfs.h')
#ROOT.gSystem.Load('../Utilities/HZGRooPdfs_cxx.so')

def makeCombs(combs):
    genSettings = configs["general_settings"]
    cats          = genSettings["categories"]
    flavs         = genSettings["flavors"]
    years         = genSettings["years"]
    if combs[0]==1:
        cats = ["comb"]
    if combs[1]==1:
        flavs = ["comb"]
    if combs[2]==1:
        years = ["comb"]
    merges = [(a,b,c) for a in cats for b in flavs for c in years]
    return merges

def lauNPowerFinder():
    #----------------------Input parsing----------------------#
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)

    isTest = 0
    nTerms = 2
    #temporary, consider moving to __init__

    if(nTerms > 5):
        exit("Error::Laurent max terms should be no greater than 5")


    genSettings = configs["general_settings"]
    cats      = genSettings["categories"]
    flavs     = genSettings["flavors"]
    years     = genSettings["years"]
    prods     = genSettings["modes"]
    binFactor = genSettings["bin_factor"]
    sideband  = genSettings["blinding_window"]
    rangeInfo = genSettings["range_info"]

    bkgSettings = configs["bkg_settings"]
    splits        = bkgSettings["splits"] #1 split, 0 means combine. [categories, flavors, years]
    ruiInBasePath = bkgSettings["base_path"]
    ggfPath       = bkgSettings["ggf_path"]
    vbfPath       = bkgSettings["vbf_path"]

    fitParams = makeCombs(splits)

    chi2Settings  = configs["lau_power_finder_settings"]
    logFile       = chi2Settings["log_file"]
    nTerms        = chi2Settings["lau_terms"]
    highPower     = chi2Settings["high_power"]
    lowPower      = chi2Settings["low_power"]
    startRange    = chi2Settings["start_range"]
    widthRange    = chi2Settings["width"]
    nllBool       = chi2Settings["nll"]

    logging.basicConfig(filename=logFile, level=logging.INFO)

    bdtSettings = configs["bdt_settings"]
    ggfBins = bdtSettings["ggf_bins"]
    vbfBins = bdtSettings["vbf_bins"]


    if highPower <= lowPower:
        raise Exception('Higher power must be larger than lower power!')
    tic = time.perf_counter()

    if isTest == 1:
        fitParams = [fitParams[0]]

    nBins = widthRange*binFactor

    #-----------------------Loop to create histograms-------------------#
    logging.info("Creating histograms")

    histList = []
    labelList = []

    for fitParam in fitParams:
        x = ROOT.RooRealVar("x", "mllg", startRange, startRange+widthRange)
        x.setBins(nBins)
        x.setRange('left', startRange, sideband[0])
        x.setRange('right', sideband[1], startRange+widthRange)
        x.setRange('full', startRange, startRange+widthRange)
        sbSample = readRuiROOTData(x, ruiInBasePath, ggfPath, vbfPath, fitParam, ggfBins, vbfBins, rangeInfo)
        sbHist = sbSample.rdHist
        histList.append(sbHist)
        label = '_'.join([s for s in fitParam])
        labelList.append(label)

    toc = time.perf_counter()
    totalTime = ((toc - tic)/60)
    logging.info("Finished creating histograms in %.2f minutes"%totalTime)

    list_ = []
    powers = [-1,-1,-1,-1,-1]

    #------------------Fit Histograms in Turn--------------------#
    if isTest:
        tic = time.perf_counter()
        i = 0
        counter = 0
        sampleDetails = labelList[i].split('_')
        lowPower = 3
        highPower = 6
        nTerms = 2
        for power1 in range(int(lowPower), int(highPower)):
            for power2 in range(power1+1, int(highPower)):
                logging.info("Beginning test fit on "+labelList[i]+" with powers" + str(power1) + " and " + str(power2))
                powers[0]=-1*power1
                powers[1]=-1*power2
                mu_gauss = ROOT.RooRealVar("mu_gauss","always 0"       ,0.)
                logging.info(str(sampleDetails))
                lauN_model = LauNClass(x, mu_gauss, sampleDetails[0], sigma_init = 4., step_init = 105, powers=powers, terms = 2, f_init = [0.1,0.1,0.1,0.1,0.1], xmax = startRange+widthRange, const_f1 = False, fix_sigma = False)
                #lauN_model = Lau2Class(x, mu_gauss, sampleDetails[0], sigma_init = 4., sigma2_init = 5., step_init = 105, p1 = powers[0], p2 = powers[1], f1_init = 0.1, f2_init = 0.1, xmax = startRange+widthRange, const_f1 = False, di_gauss = False, fix_sigma = False, gc_init = 1)
                cuthistogram = histList[i].reduce(ROOT.RooFit.CutRange('left,right'))
                if nllBool == True:
                    chi2 = lauN_model.SBpdf.createNLL(cuthistogram)
                else:
                    chi2 = lauN_model.SBpdf.createChi2(cuthistogram)
                Minimizer_Chi2(chi2, -1, 100, False, 0)
                res=Minimizer_Chi2(chi2, 1, 0.1, True, 0)
                lauN_model.checkBond()
                res.Print('v')

                #create chi2 for probability checking, but not for fitting
                chi2_ = lauN_model.SBpdf.createChi2(cuthistogram)
                chi2_val = chi2_.getVal()
                chi2_pV= ROOT.Math.chisquared_cdf_c(chi2_val, nBins - res.floatParsFinal().getSize())
                print ('Chi2 = ', chi2_val)
                print ('P-value = ', chi2_pV)
                plotClass(x, histList[i], lauN_model.pdf, lauN_model.SBpdf, labelList[i]+"pow"+str(powers[0])+"-"+str(powers[1]), "./Chi2_test/plots/", True, 'left,right')
                if(counter == 0 or chi2_val < minChi2):
                    minChi2 = chi2_val
                    powers = [powers[0],powers[1]]
                counter+=1

 

    for i in range(len(histList)):
        sampleDetails = labelList[i].split('_')
        if isTest:
            break

        logging.info("Beginning fit on "+labelList[i]+" with " + str(nTerms) + " terms")
        chi2_list = []
        pv_list = []

        values = range(int(lowPower), int(highPower)+1)
        powerCombs = list(itertools.combinations(values, nTerms))
        logging.info(str(powerCombs))

        for powers in powerCombs:#loop over possible values of that power term.
            powers_ = list(powers)
            powers = [i*-1 for i in powers_]
            logging.info("Testing " + str(powers) + " as Laurent powers")
            mu_gauss = ROOT.RooRealVar("mu_gauss","always 0"       ,0.)
            if(len(powers)<5):
                filler = [-1]*(5-len(powers))
                powers.extend(filler)
            lauN_model = LauNClass(x, mu_gauss, sampleDetails[0], sigma_init = 4., step_init = 105, powers=powers, terms = nTerms, f_init = [0.1,0.1,0.1,0.1,0.1], xmax = startRange+widthRange, const_f1 = False, fix_sigma = False)
            cuthistogram = histList[i].reduce(ROOT.RooFit.CutRange('left,right'))
            if nllBool == True:
                chi2 = lauN_model.SBpdf.createNLL(cuthistogram)
            else:
                chi2 = lauN_model.SBpdf.createChi2(cuthistogram)
            Minimizer_Chi2(chi2, -1, 100, False, 0)
            res=Minimizer_Chi2(chi2, -1, 0.1, True, 0)
            lauN_model.checkBond()
            res.Print('v')

            #create chi2 for probability checking, but not for fitting
            chi2_ = lauN_model.SBpdf.createChi2(cuthistogram)
            chi2_dof = chi2_.getVal()/(nBins - res.floatParsFinal().getSize() - 1.0)
            chi2_pV = ROOT.Math.chisquared_cdf_c(chi2_.getVal(), nBins - res.floatParsFinal().getSize())

            chi2_list.append(chi2_dof)
            pv_list.append(chi2_pV)
        logging.info("The list of powers and chi squared values are: "+ str(powerCombs) + " and " + str(chi2_list))
        logging.info("P values: " + str(pv_list))
        min_chi2 = min(chi2_list)
        index_ = chi2_list.index(min_chi2)
        logging.info("Best set of powers was " + str(powerCombs[index_]) + " with chi2 " + str(min_chi2))
        logging.info("Completed laurent power optimization for " + labelList[i] + " with " + str(nTerms) + " terms.")
        toc = time.perf_counter()
        totalTime = ((toc - tic)/60)
        logging.info('Time spent: %.2f'%totalTime + 'min')

    toc = time.perf_counter()
    totalTime = ((toc - tic)/60)
    logging.info('Time spent: %.2f'%totalTime + 'min')


if __name__=="__main__":
    print("Beginning Laurent power finding. . .")

    parser = argparse.ArgumentParser()
    parser.add_argument('-j', '--jsonIn', help = 'Fitting settings in the form of json file name')
    args = parser.parse_args()
    fname = args.jsonIn

    jfile = open(fname)
    configs = json.load(jfile)
    chi2Settings = configs["lau_power_finder_settings"]
    logFile = chi2Settings["log_file"]

    clearing = open(logFile, 'w')
    clearing.close()

    logging.basicConfig(level=logging.INFO,stream=sys.stdout)
    ROOT.gSystem.RedirectOutput(logFile)

    lauNPowerFinder()


