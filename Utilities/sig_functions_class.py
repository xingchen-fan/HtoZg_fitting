import ROOT
import json
import sys

class DSCB_Class:
    def __init__(self, x, MH, cat = "", sigmaL_init = 1.2, sigmaR_init = 1.2, nL_init = 4, nL_bond = 100, nR_init = 4, nR_bond = 100, alphaL_init = 0.5, alphaR_init = 0.5, di_sigma = False, MH_init = 125):
        self.disigma = di_sigma
        self.sigmaL = ROOT.RooRealVar("sigmaL_"+cat,"sigmaL_"+cat       , sigmaL_init, 0.01, 5.)
        self.sigmaR = ROOT.RooRealVar("sigmaR_"+cat,"sigmaR_"+cat       , sigmaR_init, 0.01, 5.)
        self.nL =  ROOT.RooRealVar("nL_"+cat, "nL_"+cat, nL_init, 0.01, nL_bond)
        self.nR =  ROOT.RooRealVar("nR_"+cat, "nR_"+cat, nR_init, 0.01, nR_bond)
        self.alphaL = ROOT.RooRealVar("alphaL_"+cat, "alphaL_"+cat, alphaL_init, 0.01, 5.)
        self.alphaR = ROOT.RooRealVar("alphaR_"+cat, "alphaR_"+cat, alphaR_init, 0.01, 5.)
        self.dMH = ROOT.RooRealVar("dMH_"+cat, "dMH_"+cat, MH_init-MH.getVal(), -3, 3)
        self.mean = ROOT.RooFormulaVar("DSCB_mean_"+cat, "@0+@1", ROOT.RooArgList(MH, self.dMH))
        if self.disigma:
            self.pdf = ROOT.RooCrystalBall.RooCrystalBall("model_DS_"+cat, "model_DS_"+cat, x, self.mean, self.sigmaL, self.sigmaR, self.alphaL, self.nL, self.alphaR, self.nR)
        else:
            self.pdf = ROOT.RooCrystalBall.RooCrystalBall("model_DS_"+cat, "model_DS_"+cat, x, self.mean, self.sigmaL, self.alphaL, self.nL, self.alphaR, self.nR)

    def setStable(self):
        if self.nL.getVal() > 30 or self.nL.getError() > 20:
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
        self.dMH.setConstant(constant)

    def assignVal(self, MH, config='', year="", cat="", lep="", prod="", debug = False):
        jfile_ = open(config, 'r')
        configs_ = json.load(jfile_)
        setting = configs_[cat]
        self.sigmaL.setVal(setting[year][lep]["sigmaL "+prod])
        self.sigmaR.setVal(setting[year][lep]["sigmaR "+prod]) 
        self.nL.setVal(setting[year][lep]["nL "+prod])
        self.nR.setVal(setting[year][lep]["nR "+prod])
        self.alphaL.setVal(setting[year][lep]["alphaL "+prod])
        self.alphaR.setVal(setting[year][lep]["alphaR "+prod])
        self.dMH.setVal(setting[year][lep]["MH "+prod] - MH.getVal())
        if debug:
            print("sigmaL = ",  self.sigmaL.getVal())
            print("sigmaR = ",  self.sigmaR.getVal())
            print("nL = ",  self.nL.getVal())
            print("nR = ",  self.nR.getVal())
            print("aL = ",  self.alphaL.getVal())
            print("aR = ",  self.alphaR.getVal())
        self.setConst(True)
        
                

class combineSignal:
    def __init__(self, x, MH, cat='', config=''):
        jfile_ = open(config, 'r')
        configs_ = json.load(jfile_)
        self.setting =  configs_[cat]
        self.ntot = 0
        for era in ["2016", "2016APV", "2017", "2018", "2022", "2022EE", "2023", "2023BPix"]:
            for lep in ["el", "mu"]:
                for prod in ["ggf", "vbf"]:
                    self.ntot += self.setting[era][lep]["nexp "+prod]
        self.c1_ggf = ROOT.RooRealVar('c2016_'+cat+'_el_ggf','c2016_'+cat+'_el_ggf', self.setting["2016"]["el"]["nexp ggf"]/self.ntot)
        self.c2_ggf = ROOT.RooRealVar('c2016_'+cat+'_mu_ggf','c2016_'+cat+'_mu_ggf', self.setting["2016"]["mu"]["nexp ggf"]/self.ntot)
        self.c3_ggf = ROOT.RooRealVar('c2016APV_'+cat+'_el_ggf','c2016APV_'+cat+'_el_ggf', self.setting["2016APV"]["el"]["nexp ggf"]/self.ntot)
        self.c4_ggf = ROOT.RooRealVar('c2016APV_'+cat+'_mu_ggf','c2016APV_'+cat+'_mu_ggf', self.setting["2016APV"]["mu"]["nexp ggf"]/self.ntot)
        self.c5_ggf = ROOT.RooRealVar('c2017_'+cat+'_el_ggf','c2017_'+cat+'_el_ggf', self.setting["2017"]["el"]["nexp ggf"]/self.ntot)
        self.c6_ggf = ROOT.RooRealVar('c2017_'+cat+'_mu_ggf','c2017_'+cat+'_mu_ggf', self.setting["2017"]["mu"]["nexp ggf"]/self.ntot)
        self.c7_ggf = ROOT.RooRealVar('c2018_'+cat+'_el_ggf','c2018_'+cat+'_el_ggf', self.setting["2018"]["el"]["nexp ggf"]/self.ntot)
        self.c8_ggf = ROOT.RooRealVar('c2018_'+cat+'_mu_ggf','c2018_'+cat+'_mu_ggf', self.setting["2018"]["mu"]["nexp ggf"]/self.ntot)
        self.c9_ggf = ROOT.RooRealVar('c2022_'+cat+'_el_ggf','c2022_'+cat+'_el_ggf', self.setting["2022"]["el"]["nexp ggf"]/self.ntot)
        self.c10_ggf = ROOT.RooRealVar('c2022_'+cat+'_mu_ggf','c2022_'+cat+'_mu_ggf', self.setting["2022"]["mu"]["nexp ggf"]/self.ntot)
        self.c11_ggf = ROOT.RooRealVar('c2022EE_'+cat+'_el_ggf','c2022EE_'+cat+'_el_ggf', self.setting["2022EE"]["el"]["nexp ggf"]/self.ntot)
        self.c12_ggf = ROOT.RooRealVar('c2022EE_'+cat+'_mu_ggf','c2022EE_'+cat+'_mu_ggf', self.setting["2022EE"]["mu"]["nexp ggf"]/self.ntot)
        self.c13_ggf = ROOT.RooRealVar('c2023_'+cat+'_el_ggf','c2023_'+cat+'_el_ggf', self.setting["2023"]["el"]["nexp ggf"]/self.ntot)
        self.c14_ggf = ROOT.RooRealVar('c2023_'+cat+'_mu_ggf','c2023_'+cat+'_mu_ggf', self.setting["2023"]["mu"]["nexp ggf"]/self.ntot)
        self.c15_ggf = ROOT.RooRealVar('c2023BPix_'+cat+'_el_ggf','c2023BPix_'+cat+'_el_ggf', self.setting["2023BPix"]["el"]["nexp ggf"]/self.ntot)
        self.c16_ggf = ROOT.RooRealVar('c2023BPix_'+cat+'_mu_ggf','c2023BPix_'+cat+'_mu_ggf', self.setting["2023BPix"]["mu"]["nexp ggf"]/self.ntot)

        self.c1_vbf = ROOT.RooRealVar('c2016_'+cat+'_el_vbf','c2016_'+cat+'_el_vbf', self.setting["2016"]["el"]["nexp vbf"]/self.ntot)
        self.c2_vbf = ROOT.RooRealVar('c2016_'+cat+'_mu_vbf','c2016_'+cat+'_mu_vbf', self.setting["2016"]["mu"]["nexp vbf"]/self.ntot)
        self.c3_vbf = ROOT.RooRealVar('c2016APV_'+cat+'_el_vbf','c2016APV_'+cat+'_el_vbf', self.setting["2016APV"]["el"]["nexp vbf"]/self.ntot)
        self.c4_vbf = ROOT.RooRealVar('c2016APV_'+cat+'_mu_vbf','c2016APV_'+cat+'_mu_vbf', self.setting["2016APV"]["mu"]["nexp vbf"]/self.ntot)
        self.c5_vbf = ROOT.RooRealVar('c2017_'+cat+'_el_vbf','c2017_'+cat+'_el_vbf', self.setting["2017"]["el"]["nexp vbf"]/self.ntot)
        self.c6_vbf = ROOT.RooRealVar('c2017_'+cat+'_mu_vbf','c2017_'+cat+'_mu_vbf', self.setting["2017"]["mu"]["nexp vbf"]/self.ntot)
        self.c7_vbf = ROOT.RooRealVar('c2018_'+cat+'_el_vbf','c2018_'+cat+'_el_vbf', self.setting["2018"]["el"]["nexp vbf"]/self.ntot)
        self.c8_vbf = ROOT.RooRealVar('c2018_'+cat+'_mu_vbf','c2018_'+cat+'_mu_vbf', self.setting["2018"]["mu"]["nexp vbf"]/self.ntot)
        self.c9_vbf = ROOT.RooRealVar('c2022_'+cat+'_el_vbf','c2022_'+cat+'_el_vbf', self.setting["2022"]["el"]["nexp vbf"]/self.ntot)
        self.c10_vbf = ROOT.RooRealVar('c2022_'+cat+'_mu_vbf','c2022_'+cat+'_mu_vbf', self.setting["2022"]["mu"]["nexp vbf"]/self.ntot)
        self.c11_vbf = ROOT.RooRealVar('c2022EE_'+cat+'_el_vbf','c2022EE_'+cat+'_el_vbf', self.setting["2022EE"]["el"]["nexp vbf"]/self.ntot)
        self.c12_vbf = ROOT.RooRealVar('c2022EE_'+cat+'_mu_vbf','c2022EE_'+cat+'_mu_vbf', self.setting["2022EE"]["mu"]["nexp vbf"]/self.ntot)
        self.c13_vbf = ROOT.RooRealVar('c2023_'+cat+'_el_vbf','c2023_'+cat+'_el_vbf', self.setting["2023"]["el"]["nexp vbf"]/self.ntot)
        self.c14_vbf = ROOT.RooRealVar('c2023_'+cat+'_mu_vbf','c2023_'+cat+'_mu_vbf', self.setting["2023"]["mu"]["nexp vbf"]/self.ntot)
        self.c15_vbf = ROOT.RooRealVar('c2023BPix_'+cat+'_el_vbf','c2023BPix_'+cat+'_el_vbf', self.setting["2023BPix"]["el"]["nexp vbf"]/self.ntot)
        #self.c16_vbf = ROOT.RooRealVar('c2023BPix_'+cat+'_mu_vbf','c2023BPix_'+cat+'_mu_vbf', self.setting["2023BPix"]["mu"]["nexp vbf"]/self.ntot)


        self.model_2016_el_ggf = DSCB_Class(x, MH, cat+"_2016_el_ggf", di_sigma = self.setting["2016"]["el"]["disigma ggf"])
        self.model_2016APV_el_ggf = DSCB_Class(x, MH, cat+"_2016APV_el_ggf", di_sigma = self.setting["2016APV"]["el"]["disigma ggf"])
        self.model_2017_el_ggf = DSCB_Class(x, MH, cat+"_2017_el_ggf", di_sigma = self.setting["2017"]["el"]["disigma ggf"])
        self.model_2018_el_ggf = DSCB_Class(x, MH, cat+"_2018_el_ggf", di_sigma = self.setting["2018"]["el"]["disigma ggf"])
        self.model_2022_el_ggf = DSCB_Class(x, MH, cat+"_2022_el_ggf", di_sigma = self.setting["2022"]["el"]["disigma ggf"])
        self.model_2022EE_el_ggf = DSCB_Class(x, MH, cat+"_2022EE_el_ggf", di_sigma = self.setting["2022EE"]["el"]["disigma ggf"])
        self.model_2023_el_ggf = DSCB_Class(x, MH, cat+"_2023_el_ggf", di_sigma = self.setting["2023"]["el"]["disigma ggf"])
        self.model_2023BPix_el_ggf = DSCB_Class(x, MH, cat+"_2023BPix_el_ggf", di_sigma = self.setting["2023BPix"]["el"]["disigma ggf"])
        self.model_2016_mu_ggf = DSCB_Class(x, MH, cat+"_2016_mu_ggf", di_sigma = self.setting["2016"]["mu"]["disigma ggf"])
        self.model_2016APV_mu_ggf = DSCB_Class(x, MH, cat+"_2016APV_mu_ggf", di_sigma = self.setting["2016APV"]["mu"]["disigma ggf"])
        self.model_2017_mu_ggf = DSCB_Class(x, MH, cat+"_2017_mu_ggf", di_sigma = self.setting["2017"]["mu"]["disigma ggf"])
        self.model_2018_mu_ggf = DSCB_Class(x, MH, cat+"_2018_mu_ggf", di_sigma = self.setting["2018"]["mu"]["disigma ggf"])
        self.model_2022_mu_ggf = DSCB_Class(x, MH, cat+"_2022_mu_ggf", di_sigma = self.setting["2022"]["mu"]["disigma ggf"])
        self.model_2022EE_mu_ggf = DSCB_Class(x, MH, cat+"_2022EE_mu_ggf", di_sigma = self.setting["2022EE"]["mu"]["disigma ggf"])
        self.model_2023_mu_ggf = DSCB_Class(x, MH, cat+"_2023_mu_ggf", di_sigma = self.setting["2023"]["mu"]["disigma ggf"])
        self.model_2023BPix_mu_ggf = DSCB_Class(x, MH, cat+"_2023BPix_mu_ggf", di_sigma = self.setting["2023BPix"]["mu"]["disigma ggf"])

        self.model_2016_el_vbf = DSCB_Class(x, MH, cat+"_2016_el_vbf", di_sigma = self.setting["2016"]["el"]["disigma vbf"])
        self.model_2016APV_el_vbf = DSCB_Class(x, MH, cat+"_2016APV_el_vbf", di_sigma = self.setting["2016APV"]["el"]["disigma vbf"])
        self.model_2017_el_vbf = DSCB_Class(x, MH, cat+"_2017_el_vbf", di_sigma = self.setting["2017"]["el"]["disigma vbf"])
        self.model_2018_el_vbf = DSCB_Class(x, MH, cat+"_2018_el_vbf", di_sigma = self.setting["2018"]["el"]["disigma vbf"])
        self.model_2022_el_vbf = DSCB_Class(x, MH, cat+"_2022_el_vbf", di_sigma = self.setting["2022"]["el"]["disigma vbf"])
        self.model_2022EE_el_vbf = DSCB_Class(x, MH, cat+"_2022EE_el_vbf", di_sigma = self.setting["2022EE"]["el"]["disigma vbf"])
        self.model_2023_el_vbf = DSCB_Class(x, MH, cat+"_2023_el_vbf", di_sigma = self.setting["2023"]["el"]["disigma vbf"])
        self.model_2023BPix_el_vbf = DSCB_Class(x, MH, cat+"_2023BPix_el_vbf", di_sigma = self.setting["2023BPix"]["el"]["disigma vbf"])
        self.model_2016_mu_vbf = DSCB_Class(x, MH, cat+"_2016_mu_vbf", di_sigma = self.setting["2016"]["mu"]["disigma vbf"])
        self.model_2016APV_mu_vbf = DSCB_Class(x, MH, cat+"_2016APV_mu_vbf", di_sigma = self.setting["2016APV"]["mu"]["disigma vbf"])
        self.model_2017_mu_vbf = DSCB_Class(x, MH, cat+"_2017_mu_vbf", di_sigma = self.setting["2017"]["mu"]["disigma vbf"])
        self.model_2018_mu_vbf = DSCB_Class(x, MH, cat+"_2018_mu_vbf", di_sigma = self.setting["2018"]["mu"]["disigma vbf"])
        self.model_2022_mu_vbf = DSCB_Class(x, MH, cat+"_2022_mu_vbf", di_sigma = self.setting["2022"]["mu"]["disigma vbf"])
        self.model_2022EE_mu_vbf = DSCB_Class(x, MH, cat+"_2022EE_mu_vbf", di_sigma = self.setting["2022EE"]["mu"]["disigma vbf"])
        self.model_2023_mu_vbf = DSCB_Class(x, MH, cat+"_2023_mu_vbf", di_sigma = self.setting["2023"]["mu"]["disigma vbf"])
        self.model_2023BPix_mu_vbf = DSCB_Class(x, MH, cat+"_2023BPix_mu_vbf", di_sigma = self.setting["2023BPix"]["mu"]["disigma vbf"])
        
        
        #model_el = DSCB_Class(x, MH, "2016_el")
        #model_mu = DSCB_Class(x, MH, "2016_mu")
        #model_el.assignVal(config, "2016", cat, "el")
        #model_mu.assignVal(config, "2016", cat, "mu")

        self.model_2016_el_ggf.assignVal(MH, config, "2016", cat, "el", "ggf")
        self.model_2016APV_el_ggf.assignVal(MH, config, "2016APV", cat, "el", "ggf")
        self.model_2017_el_ggf.assignVal(MH, config, "2017", cat, "el", "ggf")
        self.model_2018_el_ggf.assignVal(MH, config, "2018", cat, "el", "ggf")
        self.model_2022_el_ggf.assignVal(MH, config, "2022", cat, "el", "ggf")
        self.model_2022EE_el_ggf.assignVal(MH, config, "2022EE", cat, "el", "ggf")
        self.model_2023_el_ggf.assignVal(MH, config, "2023", cat, "el", "ggf")
        self.model_2023BPix_el_ggf.assignVal(MH, config, "2023BPix", cat, "el", "ggf")
        self.model_2016_mu_ggf.assignVal(MH, config, "2016", cat, "mu", "ggf")
        self.model_2016APV_mu_ggf.assignVal(MH, config, "2016APV", cat, "mu", "ggf")
        self.model_2017_mu_ggf.assignVal(MH, config, "2017", cat, "mu", "ggf")
        self.model_2018_mu_ggf.assignVal(MH, config, "2018", cat, "mu", "ggf")
        self.model_2022_mu_ggf.assignVal(MH, config, "2022", cat, "mu", "ggf")
        self.model_2022EE_mu_ggf.assignVal(MH, config, "2022EE", cat, "mu", "ggf")
        self.model_2023_mu_ggf.assignVal(MH, config, "2023", cat, "mu", "ggf")
        self.model_2023BPix_mu_ggf.assignVal(MH, config, "2023BPix", cat, "mu", "ggf")

        self.model_2016_el_vbf.assignVal(MH, config, "2016", cat, "el", "vbf")
        self.model_2016APV_el_vbf.assignVal(MH, config, "2016APV", cat, "el", "vbf")
        self.model_2017_el_vbf.assignVal(MH, config, "2017", cat, "el", "vbf")
        self.model_2018_el_vbf.assignVal(MH, config, "2018", cat, "el", "vbf")
        self.model_2022_el_vbf.assignVal(MH, config, "2022", cat, "el", "vbf")
        self.model_2022EE_el_vbf.assignVal(MH, config, "2022EE", cat, "el", "vbf")
        self.model_2023_el_vbf.assignVal(MH, config, "2023", cat, "el", "vbf")
        self.model_2023BPix_el_vbf.assignVal(MH, config, "2023BPix", cat, "el", "vbf")
        self.model_2016_mu_vbf.assignVal(MH, config, "2016", cat, "mu", "vbf")
        self.model_2016APV_mu_vbf.assignVal(MH, config, "2016APV", cat, "mu", "vbf")
        self.model_2017_mu_vbf.assignVal(MH, config, "2017", cat, "mu", "vbf")
        self.model_2018_mu_vbf.assignVal(MH, config, "2018", cat, "mu", "vbf")
        self.model_2022_mu_vbf.assignVal(MH, config, "2022", cat, "mu", "vbf")
        self.model_2022EE_mu_vbf.assignVal(MH, config, "2022EE", cat, "mu", "vbf")
        self.model_2023_mu_vbf.assignVal(MH, config, "2023", cat, "mu", "vbf")
        self.model_2023BPix_mu_vbf.assignVal(MH, config, "2023BPix", cat, "mu", "vbf")


        #self.pdf = ROOT.RooAddPdf('combine_sig_'+cat, 'combine_sig_'+cat, pdf_list, ROOT.RooArgList(self.c1, self.c2, self.c3, self.c4, self.c5, self.c6, self.c7, self.c8, self.c9, self.c10, self.c11, self.c12, self.c13, self.c14, self.c15, self.c16), False)
        self.pdf = ROOT.RooAddPdf("combine_sig_"+cat, "combine_sig_"+cat, ROOT.RooArgList(self.model_2016_el_ggf.pdf, self.model_2016_mu_ggf.pdf, self.model_2016APV_el_ggf.pdf, self.model_2016APV_mu_ggf.pdf, self.model_2017_el_ggf.pdf, self.model_2017_mu_ggf.pdf, self.model_2018_el_ggf.pdf, self.model_2018_mu_ggf.pdf, self.model_2022_el_ggf.pdf, self.model_2022_mu_ggf.pdf, self.model_2022EE_el_ggf.pdf, self.model_2022EE_mu_ggf.pdf, self.model_2023_el_ggf.pdf, self.model_2023_mu_ggf.pdf, self.model_2023BPix_el_ggf.pdf, self.model_2023BPix_mu_ggf.pdf, self.model_2016_el_vbf.pdf, self.model_2016_mu_vbf.pdf, self.model_2016APV_el_vbf.pdf, self.model_2016APV_mu_vbf.pdf, self.model_2017_el_vbf.pdf, self.model_2017_mu_vbf.pdf, self.model_2018_el_vbf.pdf, self.model_2018_mu_vbf.pdf, self.model_2022_el_vbf.pdf, self.model_2022_mu_vbf.pdf, self.model_2022EE_el_vbf.pdf, self.model_2022EE_mu_vbf.pdf, self.model_2023_el_vbf.pdf, self.model_2023_mu_vbf.pdf, self.model_2023BPix_el_vbf.pdf, self.model_2023BPix_mu_vbf.pdf), ROOT.RooArgList(self.c1_ggf, self.c2_ggf, self.c3_ggf, self.c4_ggf, self.c5_ggf, self.c6_ggf, self.c7_ggf, self.c8_ggf, self.c9_ggf, self.c10_ggf, self.c11_ggf, self.c12_ggf, self.c13_ggf, self.c14_ggf, self.c15_ggf, self.c16_ggf, self.c1_vbf, self.c2_vbf, self.c3_vbf, self.c4_vbf, self.c5_vbf, self.c6_vbf, self.c7_vbf, self.c8_vbf, self.c9_vbf, self.c10_vbf, self.c11_vbf, self.c12_vbf, self.c13_vbf, self.c14_vbf, self.c15_vbf))
        
