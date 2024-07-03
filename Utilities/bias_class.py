import ROOT
from Xc_Minimizer import *

class BiasClass:
    def __init__(self, pdf, hist, sideBand = False, fitRange = "", extend = False):
        if sideBand:
            self.nll = ROOT.RooNLLVar("nll_"+pdf.GetName(), "nll_"+pdf.GetName(), pdf, hist, ROOT.RooFit.Range(fitRange))
        else:
            self.nll = ROOT.RooNLLVar("nll_"+pdf.GetName(), "nll_"+pdf.GetName(), pdf, hist, ROOT.RooFit.Extended(extend))
        
        self.corrNLL = 0.
    def minimize(self, printLevel = -1, eps = 0.1, offSet = True, skip_hesse = False):
        Minimizer_NLL(self.nll, -1, 100, False, 0, skip_hesse)
        r = Minimizer_NLL(self.nll, printLevel, eps, offSet, 0, skip_hesse)
        #print("Cov q = ", r.covQual(), " status = ", r.status(), end="\n")
        if r.status() != 0: 
            print(self.nll.GetName(), " Minimization fails!")
            r.Print("V")
        #r.Print("V")
        self.corrNLL = self.nll.getVal() + 0.5*r.floatParsFinal().getSize()

        
