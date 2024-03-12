import ROOT

class DSCB_Class:
    def __init__(self, x, MH, cat = "", sigma_init = 1.2, nL_init = 4, nL_bond = 50, nR_init = 4, nR_bond = 50, alphaL_init = 0.5, alphaR_init = 0.5):
        self.sigma = ROOT.RooRealVar("sigma_"+cat,"sigma_"+cat       , sigma_init, 0.01, 5.)
        self.nL =  ROOT.RooRealVar("nL_"+cat, "nL_"+cat, nL_init, 0.01, nL_bond)
        self.nR =  ROOT.RooRealVar("nR_"+cat, "nR_"+cat, nR_init, 0.01, nR_bond)
        self.alphaL = ROOT.RooRealVar("alphaL_"+cat, "alphaL_"+cat, alphaL_init, 0.01, 5.)
        self.alphaR = ROOT.RooRealVar("alphaR_"+cat, "alphaR_"+cat, alphaR_init, 0.01, 5.)

        self.pdf = ROOT.RooCrystalBall.RooCrystalBall("sig_model_DS_"+cat, "sig_model_DS_"+cat, x, MH, self.sigma, self.alphaL, self.nL, self.alphaR, self.nR)

        def setConst(self, constant = True):
            self.sigma.setConstant(constant)
            self.nL.setConstant(constant)
            self.nR.setConstant(constant)
            self.alphaL.setConstant(constant)
            self.alphaR.setConstant(constant)


