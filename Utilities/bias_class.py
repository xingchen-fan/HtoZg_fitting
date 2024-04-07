import ROOT
from Xc_Minimizer import *

class BiasClass:
    def __init__(self, pdf, hist, sideBand = False, fitRange = ''):
        if sideBand:
            self.nll = ROOT.RooNLLVar("nll_"+pdf.GetName(), "nll_"+pdf.GetName(), pdf, hist, ROOT.RooFit.Range(fitRange))
        else:
            self.nll = ROOT.RooNLLVar("nll_"+pdf.GetName(), "nll_"+pdf.GetName(), pdf, hist)
        self.corrNLL = 0.
    def minimize(self, printLevel = -1, eps = 1, offSet = True):
        Minimizer_NLL(self.nll, -1, 100, False, 0)
        r = Minimizer_NLL(self.nll, printLevel, eps, offSet, 0)
        print("Cov q = ", r.covQual(), " status = ", r.status(), end="\n")
        #r.Print("V")
        self.corrNLL = self.nll.getVal() + 0.5*r.floatParsFinal().getSize()

        
