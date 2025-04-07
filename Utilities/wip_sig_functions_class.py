import ROOT
import json
import sys

class DSCB_Class:
    def __init__(self, x, MH, cat = "", sigmaL_init = 1.2, sigmaR_init = 1.2, nL_init = 4, nL_bond = 100, nR_init = 4, nR_bond = 100, alphaL_init = 0.5, alphaR_init = 0.5, di_sigma = False):
        self.disigma = di_sigma
        self.sigmaL = ROOT.RooRealVar("sigmaL_"+cat,"sigmaL_"+cat       , sigmaL_init, 0.01, 5.)
        self.sigmaR = ROOT.RooRealVar("sigmaR_"+cat,"sigmaR_"+cat       , sigmaR_init, 0.01, 5.)
        self.nL =  ROOT.RooRealVar("nL_"+cat, "nL_"+cat, nL_init, 0.01, nL_bond)
        self.nR =  ROOT.RooRealVar("nR_"+cat, "nR_"+cat, nR_init, 0.01, nR_bond)
        self.alphaL = ROOT.RooRealVar("alphaL_"+cat, "alphaL_"+cat, alphaL_init, 0.01, 5.)
        self.alphaR = ROOT.RooRealVar("alphaR_"+cat, "alphaR_"+cat, alphaR_init, 0.01, 5.)
        if self.disigma:
            self.pdf = ROOT.RooCrystalBall.RooCrystalBall("combine_sig_"+cat, "model_DS_"+cat, x, MH, self.sigmaL, self.sigmaR, self.alphaL, self.nL, self.alphaR, self.nR)
        else:
            self.pdf = ROOT.RooCrystalBall.RooCrystalBall("combine_sig_"+cat, "model_DS_"+cat, x, MH, self.sigmaL, self.alphaL, self.nL, self.alphaR, self.nR)

    def setStable(self):
        if self.nL.getVal() > 30 or self.nL.getError() > 50:
            self.nL.setVal(20)
            self.nL.setError(0)
            self.nL.setConstant(True)
        if self.nR.getVal() > 30 or self.nR.getError() > 50:
            self.nR.setVal(20)
            self.nR.setError(0)
            self.nR.setConstant(True)
        if self.nL.getVal() < 0.1:
            self.nL.setVal(0)
            self.nL.setError(0)
            self.nL.setConstant(True)
        if self.nR.getVal() < 0.1:
            self.nR.setVal(0)
            self.nR.setError(0)
            self.nR.setConstant(True)
    
    def setConst(self, constant = True):
        self.sigmaL.setConstant(constant)
        self.sigmaR.setConstant(constant)
        self.nL.setConstant(constant)
        self.nR.setConstant(constant)
        self.alphaL.setConstant(constant)
        self.alphaR.setConstant(constant)

    def assignVal(self, config='', year="", cat="", lep="", prod="", debug = False):
        jfile_ = open(config, 'r')
        configs_ = json.load(jfile_)
        setting = configs_[cat]
        self.sigmaL.setVal(setting[year][lep]["sigmaL "+prod])
        self.sigmaR.setVal(setting[year][lep]["sigmaR "+prod]) 
        self.nL.setVal(setting[year][lep]["nL "+prod])
        self.nR.setVal(setting[year][lep]["nR "+prod])
        self.alphaL.setVal(setting[year][lep]["alphaL "+prod])
        self.alphaR.setVal(setting[year][lep]["alphaR "+prod])
        self.disigma = setting[year][lep]["disigma "+prod]
        if debug:
            print("sigmaL = ",  self.sigmaL.getVal())
            print("sigmaR = ",  self.sigmaR.getVal())
            print("nL = ",  self.nL.getVal())
            print("nR = ",  self.nR.getVal())
            print("aL = ",  self.alphaL.getVal())
            print("aR = ",  self.alphaR.getVal())
        self.setConst(True)
       
