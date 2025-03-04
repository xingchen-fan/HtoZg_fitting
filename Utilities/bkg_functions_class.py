import ROOT

def bondComp(par, tol):
    if par.getValV()+tol > par.getMax() or par.getVal()-tol < par.getMin():
        print (par.GetName(), " = ", par.getValV(), " hits the bondary (", par.getMin(), ", ", par.getMax(), ")")
        return True
    else: return False
    
#### Bernstein
class Bern2Class:
    def __init__(self, x, gauss_mu, cat="", p0=10, p_init=0.3, bond=20, sigma_init=7., step_init=105.):
        self.init_list = [p0, p_init, sigma_init, step_init]
        self.p0 = ROOT.RooRealVar("b2p0_" + cat, "b2p1_" + cat, p0)
        self.p1 = ROOT.RooRealVar("b2p1_" + cat, "b2p1_" + cat, p_init,-bond, bond)
        self.p2 = ROOT.RooRealVar("b2p2_" + cat, "b2p2_" + cat, p_init,-bond, bond)
        self.sigma = ROOT.RooRealVar("sigma_bern2_" + cat,"sigma_bern2_" + cat,sigma_init,  0.1, 15.)
        self.stepval = ROOT.RooRealVar("stepval_bern2_" + cat, "stepval_bern2_" + cat, step_init, 100., 120.)
        self.pdf = ROOT.RooGaussStepBernstein("bern2_" +cat + "_model", "Bernstein2 (X) gauss " + cat, x,  gauss_mu, self.sigma, self.stepval, ROOT.RooArgList(self.p0,self.p1,self.p2))
        self.SBpdf = ROOT.RooGenericPdf("SBbern2_" +cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
        self.name = "bern2_"+ cat
    def checkBond(self):
        tol = 0.001
        par_list = [self.p1, self.p2, self.sigma, self.stepval]
        if any(bondComp(par, tol) for par in par_list):
            print ("The pdf ", self.pdf.GetName(), " needs refit.")
    def reset(self):
         self.p0.setVal(self.init_list[0])
         self.p1.setVal(self.init_list[1])
         self.p2.setVal(self.init_list[1])
         self.sigma.setVal(self.init_list[2])
         self.stepval.setVal(self.init_list[3])

class Bern3Class:
    def __init__(self, x, gauss_mu, cat="", p0=10, p_init=0.3, bond=20, sigma_init=7., step_init=105.):
        self.init_list = [p0, p_init, sigma_init, step_init]
        self.p0 = ROOT.RooRealVar("b3p0_" + cat, "b3p1_" + cat, p0)
        self.p1 = ROOT.RooRealVar("b3p1_" + cat, "b3p1_" + cat, p_init,-bond, bond)
        self.p2 = ROOT.RooRealVar("b3p2_" + cat, "b3p2_" + cat, p_init,-bond, bond)
        self.p3 = ROOT.RooRealVar("b3p3_" + cat, "b3p3_" + cat, p_init,-bond, bond)
        self.sigma = ROOT.RooRealVar("sigma_bern3_" + cat,"sigma_bern3_" + cat, sigma_init,  0.1, 15.)
        self.stepval = ROOT.RooRealVar("stepval_bern3_" + cat, "stepval_bern3_" + cat, step_init, 100., 120.)
        self.pdf = ROOT.RooGaussStepBernstein("bern3_" +cat + "_model", "Bernstein3 (X) gauss " + cat, x,  gauss_mu, self.sigma, self.stepval, ROOT.RooArgList(self.p0,self.p1,self.p2,self.p3))
        self.SBpdf = ROOT.RooGenericPdf("SBbern3_" +cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
        self.name = "bern3_"+ cat
    def checkBond(self):
        tol = 0.001
        par_list = [self.p1, self.p2, self.p3, self.sigma, self.stepval]
        if any(bondComp(par, tol) for par in par_list):
            print ("The pdf ", self.pdf.GetName(), " needs refit.")
    def reset(self):
         self.p0.setVal(self.init_list[0])
         self.p1.setVal(self.init_list[1])
         self.p2.setVal(self.init_list[1])
         self.p3.setVal(self.init_list[1])
         self.sigma.setVal(self.init_list[2])
         self.stepval.setVal(self.init_list[3])

class Bern4Class:
    def __init__(self, x, gauss_mu, cat="", p0=10, p_init=0.3, bond=20, sigma_init=7., step_init=105.):
        self.init_list = [p0, p_init, sigma_init, step_init]
        self.p0 = ROOT.RooRealVar("b4p0_" + cat, "b4p1_" + cat, p0)
        self.p1 = ROOT.RooRealVar("b4p1_" + cat, "b4p1_" + cat, p_init, -bond, bond)
        self.p2 = ROOT.RooRealVar("b4p2_" + cat, "b4p2_" + cat, p_init, -bond, bond)
        self.p3 = ROOT.RooRealVar("b4p3_" + cat, "b4p3_" + cat, p_init, -bond, bond)
        self.p4 = ROOT.RooRealVar("b4p4_" + cat, "b4p4_" + cat, p_init, -bond, bond)
        self.sigma = ROOT.RooRealVar("sigma_bern4_" + cat,"sigma_bern4_" + cat, sigma_init,  0.1, 15.)
        self.stepval = ROOT.RooRealVar("stepval_bern4_" + cat, "stepval_bern4_" + cat, step_init, 100., 120.)
        self.pdf = ROOT.RooGaussStepBernstein("bern4_" +cat + "_model", "Bernstein4 (X) gauss " + cat, x,  gauss_mu, self.sigma, self.stepval, ROOT.RooArgList(self.p0,self.p1,self.p2, self.p3, self.p4))
        self.SBpdf = ROOT.RooGenericPdf("SBbern4_" +cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
        self.name = "bern4_"+ cat
    def checkBond(self):
        tol = 0.001
        par_list = [self.p1, self.p2, self.p3, self.p4, self.sigma, self.stepval]
        if any(bondComp(par, tol) for par in par_list):
            print ("The pdf ", self.pdf.GetName(), " needs refit.")
    def reset(self):
         self.p0.setVal(self.init_list[0])
         self.p1.setVal(self.init_list[1])
         self.p2.setVal(self.init_list[1])
         self.p3.setVal(self.init_list[1])
         self.p4.setVal(self.init_list[1])
         self.sigma.setVal(self.init_list[2])
         self.stepval.setVal(self.init_list[3])

class Bern5Class:
    def __init__(self, x, gauss_mu, cat="", p0=1., p_init=0.3, bond=20, sigma_init=7., step_init=105.):
        self.init_list = [p0, p_init, sigma_init, step_init]
        self.p0 = ROOT.RooRealVar("b5p0_" + cat, "b5p1_" + cat, p0)
        self.p1 = ROOT.RooRealVar("b5p1_" + cat, "b5p1_" + cat, p_init,-bond, bond)
        self.p2 = ROOT.RooRealVar("b5p2_" + cat, "b5p2_" + cat, p_init,-bond, bond)
        self.p3 = ROOT.RooRealVar("b5p3_" + cat, "b5p3_" + cat, p_init,-bond, bond)
        self.p4 = ROOT.RooRealVar("b5p4_" + cat, "b5p4_" + cat, p_init,-bond, bond)
        self.p5 = ROOT.RooRealVar("b5p5_" + cat, "b5p5_" + cat, p_init,-bond, bond)
        self.sigma = ROOT.RooRealVar("sigma_bern5_" + cat,"sigma_bern5_" + cat,sigma_init,  0.1, 15.)
        self.stepval = ROOT.RooRealVar("stepval_bern5_" + cat, "stepval_bern5_" + cat, step_init, 100., 120.)
        self.pdf = ROOT.RooGaussStepBernstein("bern5_" +cat + "_model", "Bernstein5 (X) gauss " + cat, x,  gauss_mu, self.sigma, self.stepval, ROOT.RooArgList(self.p0,self.p1,self.p2, self.p3, self.p4, self.p5))
        self.SBpdf = ROOT.RooGenericPdf("SBbern5_" +cat + "_model", "SBbern5_" +cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
        self.name = "bern5_"+ cat
    def checkBond(self):
        tol = 0.001
        par_list = [self.p1, self.p2, self.p3, self.p4, self.p5, self.sigma, self.stepval]
        if any(bondComp(par, tol) for par in par_list):
            print ("The pdf ", self.pdf.GetName(), " needs refit.")
    def reset(self):
         self.p0.setVal(self.init_list[0])
         self.p1.setVal(self.init_list[1])
         self.p2.setVal(self.init_list[1])
         self.p3.setVal(self.init_list[1])
         self.p4.setVal(self.init_list[1])
         self.p5.setVal(self.init_list[1])
         self.sigma.setVal(self.init_list[2])
         self.stepval.setVal(self.init_list[3])


#### Pow series
class Pow1Class:
    def __init__(self, x, gauss_mu, cat="", sigma_init = 7., step_init = 108., p_init = -6., p_low = -10., p_high = 0., di_gauss = False):
        self.init_list = [sigma_init, step_init, p_init]
        self.t = ROOT.RooRealVar("pow1t_" + cat, "t pow1" + cat, step_init, 100., 118.)
        self.p = ROOT.RooRealVar("pow1p_" + cat, "p1 pow1" + cat, p_init, p_low, p_high)
        self.sigma = ROOT.RooRealVar("sigma_pow1_" + cat,"sigma_pow1_"+cat, sigma_init,  1., 15.)
        self.sigma2 = ROOT.RooRealVar("sigma2_pow1_" + cat,"sigma2_pow1_"+cat, sigma_init * 2,  0.1, 15.)
        self.gauss = ROOT.RooGaussian("gaussxpow1_"+cat, "gaussian PDF pow1 " + cat, x, gauss_mu, self.sigma)
        self.gauss2 = ROOT.RooGaussian("gauss2xpow1_"+cat, "gaussian2 PDF pow1 " + cat, x, gauss_mu, self.sigma2)
        self.gc = ROOT.RooRealVar("gc_pow1_" + cat, "gc_pow1_" + cat, 1, 0, 10)
        self.addG = ROOT.RooAddPdf("addG_pow1_" + cat, "addG_pow1_" + cat, self.gauss2, self.gauss, self.gc)
        self.transX = ROOT.RooFormulaVar("transX_"+cat, "@0 - 100", ROOT.RooArgList(x)) 
        self.transt = ROOT.RooFormulaVar("transt_"+cat, "@0 - 100", ROOT.RooArgList(self.t))
        self.step = ROOT.RooGenericPdf("step_pow1_" + cat, "step_pow1_" + cat,\
                                        "( ((@0-@1)*153.85 <0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85) ) )*(@0^@2)", ROOT.RooArgList(x,self.t,self.p))
        x.setBins(20000, "cache")
        if di_gauss == True:
            self.pdf = ROOT.RooFFTConvPdf("pow1_" + cat + "_model", "step pow1 (X) gauss" + cat, x, self.step, self.addG)
        else:
            self.pdf = ROOT.RooFFTConvPdf("pow1_" + cat + "_model", "step pow1 (X) gauss" + cat, x, self.step, self.gauss)
        self.pdf.setBufferFraction(0.5)
        self.SBpdf = ROOT.RooGenericPdf("SBpow1_" +cat + "_model", "SBpow1_" +cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
        self.name = "pow1_"+ cat
    def checkBond(self):
        tol = 0.001
        par_list = [self.t, self.p, self.sigma]
        if any(bondComp(par, tol) for par in par_list):
            print ("The pdf ", self.pdf.GetName(), " needs refit.")
    def reset(self):
        self.t.setVal(self.init_list[1])
        self.p.setVal(self.init_list[2])
        self.sigma.setVal(self.init_list[0])

# class Pow2Class:
#     def __init__(self, x, gauss_mu, cat="", sigma_init = 7., step_init = 105., p1_init = -5., p2_init = -9., f_init = 0.5):
#         self.init_list = [sigma_init, step_init, p1_init, p2_init, f_init]
#         self.t = ROOT.RooRealVar("pow2t_" + cat, "t pow2" + cat, step_init, 95., 115.)
#         self.p1 = ROOT.RooRealVar("pow2p1_" + cat, "p1 pow2" + cat, p1_init, -10., 0.)
#         self.p2 = ROOT.RooRealVar("pow2p2_" + cat, "p2 pow2" + cat, p2_init, -100., -5.)
#         self.f = ROOT.RooRealVar("pow2f_" + cat, "f pow2" + cat, f_init, 0., 1.)
#         self.sigma = ROOT.RooRealVar("sigma_pow2_" + cat,"sigma_pow2_"+cat, sigma_init,  1., 15.)
#         self.gauss = ROOT.RooGaussian("gaussxpow2_"+cat, "gaussian PDF pow2 " + cat, x, gauss_mu, self.sigma)
#         # self.step = ROOT.RooGenericPdf("step_pow2_" + cat, "step_pow2" + cat,"( ((@0-@1)*153.85 <0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85/4) ) )*(@0^@2 + @4 * @0^@3)", ROOT.RooArgList(x,self.t,self.p1, self.p2, self.f))
#         self.step = ROOT.RooGenericPdf("step_pow2_" + cat, "step_pow2" + cat,"( ((@0-@1)*153.85 <0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85/4) ) )*((1-@4)/((200^(1+@2) - @1^(1+@2))/(1+@2)) * @0^@2 + @4 /((200^(1+@3) - @1^(1+@3))/(1+@3)) * @0^@3)", ROOT.RooArgList(x,self.t,self.p1, self.p2, self.f))
#         # self.step = ROOT.RooGenericPdf("step_pow2_" + cat, "step_pow2" + cat,"( ((@0-@1)*153.85 <0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85/4) ) )*((1 - @4)*@0^@2 + @4 /((165^(1+@3) - @1^(1+@3))/(1+@3)) * @0^@3)", ROOT.RooArgList(x,self.t,self.p1, self.p2, self.f))
#         # self.step = ROOT.RooGenericPdf("step_pow2_" + cat, "step_pow2" + cat,"( ((@0-@1)*153.85 <0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85/4) ) )*(((@0-67.5)/130)^@2 + @4 *((@0-67.5)/130)^@3)", ROOT.RooArgList(x,self.t,self.p1, self.p2, self.f))
#         x.setBins(20000, "cache")
#         self.pdf = ROOT.RooFFTConvPdf("pow2_" + cat + "_model", "step pow2 (X) gauss " + cat, x, self.step, self.gauss)
#         self.pdf.setBufferFraction(0.5)
#         self.name = "pow2_"+ cat
#     def checkBond(self):
#         tol = 0.0001
#         par_list = [self.t, self.p1, self.p2, self.f, self.sigma]
#         if any(bondComp(par, tol) for par in par_list):
#             print ("The pdf ", self.pdf.GetName(), " needs refit.")
#     def reset(self):
#         self.t.setVal(self.init_list[1])
#         self.p1.setVal(self.init_list[2])
#         self.p2.setVal(self.init_list[3])
#         self.f.setVal(self.init_list[4])
#         self.sigma.setVal(self.init_list[0])

# class Pow2Class:
#     def __init__(self, x, gauss_mu, cat="", sigma_init = 7., step_init = 105., p1_init = -5., p2_init = -9., f_init = 0.9):
#         self.init_list = [sigma_init, step_init, p1_init, p2_init, f_init]
#         self.t = ROOT.RooRealVar("pow2t_" + cat, "t pow2" + cat, step_init, 95., 115.OA)
#         self.p1 = ROOT.RooRealVar("pow2p1_" + cat, "p1 pow2" + cat, p1_init, -10., 1.)
#         self.p2 = ROOT.RooRealVar("pow2p2_" + cat, "p2 pow2" + cat, p2_init, -10, -5.)
#         self.f = ROOT.RooRealVar("pow2f_" + cat, "f pow2" + cat, f_init, 0., 1.)
#         self.sigma = ROOT.RooRealVar("sigma_pow2_" + cat,"sigma_pow2_"+cat, sigma_init,  1., 15.)
#         self.gauss = ROOT.RooGaussian("gaussxpow2_"+cat, "gaussian PDF pow2 " + cat, x, gauss_mu, self.sigma)
#         # self.step = ROOT.RooGenericPdf("step_pow2_" + cat, "step_pow2" + cat,"( ((@0-@1)*153.85 <0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85/4) ) )*(@0^@2 + @4 * @0^@3)", ROOT.RooArgList(x,self.t,self.p1, self.p2, self.f))
#         # self.step1 = ROOT.RooGenericPdf("step1_pow2_" + cat, "step1_pow2_" + cat,"( ((@0-@1)*153.85 <0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85/4) ) )*(@0^@2)", ROOT.RooArgList(x,self.t,self.p1))
#         # self.step2 = ROOT.RooGenericPdf("step2_pow2_" + cat, "step2_pow2_" + cat,"( ((@0-@1)*153.85 <0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85/4) ) )*(@0^@2)", ROOT.RooArgList(x,self.t,self.p2))
#         self.step1 = ROOT.RooGenericPdf("step1_pow2_" + cat, "step1_pow2_" + cat,"(@0>@1)*(@0^@2)", ROOT.RooArgList(x,self.t,self.p1))
#         self.step2 = ROOT.RooGenericPdf("step2_pow2_" + cat, "step2_pow2_" + cat,"(@0>@1)*(@0^@2)", ROOT.RooArgList(x,self.t,self.p2))
#         self.step = ROOT.RooAddPdf("step_pow2_"+ cat, "step_pow2_" + cat, self.step1, self.step2, self.f)
#         # self.step = ROOT.RooGenericPdf("step_pow2_" + cat, "step_pow2" + cat,"( ((@0-@1)*153.85 <0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85/4) ) )*((1 - @4)*@0^@2 + @4 /((165^(1+@3) - @1^(1+@3))/(1+@3)) * @0^@3)", ROOT.RooArgList(x,self.t,self.p1, self.p2, self.f))
#         # self.step = ROOT.RooGenericPdf("step_pow2_" + cat, "step_pow2" + cat,"( ((@0-@1)*153.85 <0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85/4) ) )*(((@0-67.5)/130)^@2 + @4 *((@0-67.5)/130)^@3)", ROOT.RooArgList(x,self.t,self.p1, self.p2, self.f))
#         x.setBins(20000, "cache")
#         self.pdf = ROOT.RooFFTConvPdf("pow2_" + cat + "_model", "step pow2 (X) gauss " + cat, x, self.step, self.gauss)
#         self.pdf.setBufferFraction(0.5)
#         self.name = "pow2_"+ cat
#     def checkBond(self):
#         tol = 0.0001
#         par_list = [self.t, self.p1, self.p2, self.f, self.sigma]
#         if any(bondComp(par, tol) for par in par_list):
#             print ("The pdf ", self.pdf.GetName(), " needs refit.")
#     def reset(self):
#         self.t.setVal(self.init_list[1])
#         self.p1.setVal(self.init_list[2])
#         self.p2.setVal(self.init_list[3])
#         self.f.setVal(self.init_list[4])
#         self.sigma.setVal(self.init_list[0])

class Pow2Class:
    def __init__(self, x, gauss_mu, cat="", sigma_init = 7., step_init = 107., p1_init = -5.8, p1_low = -8., p1_high = 0., p2_init = -15, p2_low = -20., p2_high = -5., f1_init = 0.9, f2_init = 0.9, xmax = 165., const_f1 = False, di_gauss = False):
        self.init_list = [sigma_init, step_init, p1_init, p2_init, f1_init, f2_init]
        self.xmax = ROOT.RooRealVar("xmax_pow2_" + cat, "xmax_pow2_" + cat, xmax)
        self.t = ROOT.RooRealVar("pow2t_" + cat, "t pow2" + cat, step_init, 95., 118.)
        self.p1 = ROOT.RooRealVar("pow2p1_" + cat, "p1 pow2" + cat, p1_init, p1_low, p1_high)
        self.p2 = ROOT.RooRealVar("pow2p2_" + cat, "p2 pow2" + cat, p2_init, p2_low, p2_high)
        self.f1 = ROOT.RooRealVar("pow2f1_" + cat, "f1 pow2" + cat, f1_init, 0., 10.)
        self.f2 = ROOT.RooRealVar("pow2f2_" + cat, "f2 pow2" + cat, f2_init, 0., 10.)
        self.sigma = ROOT.RooRealVar("sigma_pow2_" + cat,"sigma_pow2_"+cat, sigma_init,  0.1, 15.)
        self.sigma2 = ROOT.RooRealVar("sigma2_pow2_" + cat,"sigma2_pow2_"+cat, sigma_init * 2,  0.1, 15.)
        self.offset = ROOT.RooRealVar("os_" + cat, "os_" + cat, 1e-10)
        self.gauss = ROOT.RooGaussian("gaussxpow2_"+cat, "gaussian PDF pow2 " + cat, x, gauss_mu, self.sigma)
        self.gauss2 = ROOT.RooGaussian("gauss2xpow2_"+cat, "gaussian2 PDF pow2 " + cat, x, gauss_mu, self.sigma2)
        self.gc = ROOT.RooRealVar("gc_pow2_" + cat, "gc_pow2_" + cat, 1, 0, 10)
        self.addG = ROOT.RooAddPdf("addG_pow2_" + cat, "addG_pow2_" + cat, self.gauss2, self.gauss, self.gc)
        #self.step1 = ROOT.NormPow("step1_pow2_" + cat, "step1_pow2_" + cat, x,self.t,self.p1, 165)
        #self.step2 = ROOT.NormPow("step2_pow2_" + cat, "step2_pow2_" + cat,x,self.t,self.p2, 165)
        #self.step = ROOT.RooAddPdf("step_pow2_"+ cat, "step_pow2_" + cat, self.step1, self.step2, self.f1)
        #self.transX = ROOT.RooFormulaVar("transX_"+cat, "@0 - 50", ROOT.RooArgList(x))
        #self.transt = ROOT.RooFormulaVar("transt_"+cat, "@0 - 50", ROOT.RooArgList(self.t))
        if const_f1 == True:
            self.f1.setConstant(True)
        self.norm1 = ROOT.RooFormulaVar("norm1_pow2_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.t, self.p1, self.xmax))
        self.norm2 = ROOT.RooFormulaVar("norm2_pow2_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.t, self.p2, self.xmax))
        self.step = ROOT.RooGenericPdf("step_pow2_" + cat, "step_pow2" + cat,"( ((@0-@1)*153.85 <0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85) ) )*(@5/@6 * @0^@2 + @4 / @7 * @0^@3)", ROOT.RooArgList(x, self.t,self.p1, self.p2, self.f1, self.f2, self.norm1, self.norm2))
        x.setBins(20000, "cache")
        if di_gauss == True:
            self.pdf = ROOT.RooFFTConvPdf("pow2_" + cat + "_model", "step pow2 (X) gauss" + cat, x, self.step, self.addG)
        else:
            self.pdf = ROOT.RooFFTConvPdf("pow2_" + cat + "_model", "step pow2 (X) gauss" + cat, x, self.step, self.gauss)
        self.pdf.setBufferFraction(0.5)
        self.SBpdf = ROOT.RooGenericPdf("SBpow2_" +cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
        self.name = "pow2_"+ cat
    def checkBond(self):
        tol = 0.0001
        par_list = [self.t, self.p1, self.p2, self.f1, self.sigma, self.f2]
        if any(bondComp(par, tol) for par in par_list):
            print ("The pdf ", self.pdf.GetName(), " needs refit.")
    def reset(self):
        self.t.setVal(self.init_list[1])
        self.p1.setVal(self.init_list[2])
        self.p2.setVal(self.init_list[3])
        self.f1.setVal(self.init_list[4])
        self.f2.setVal(self.init_list[5])
        self.sigma.setVal(self.init_list[0])

# class Pow3Class:
#     def __init__(self, x, gauss_mu, cat="", sigma_init = 7., step_init = 105., p1_init = -5., p2_init = -7., p3_init = -9., f_init = 0.3):
#         self.init_list = [sigma_init, step_init, p1_init, p2_init, p3_init, f_init]
#         self.t = ROOT.RooRealVar("pow3t_" + cat, "t pow3" + cat, step_init, 95., 115.)
#         self.p1 = ROOT.RooRealVar("pow3p1_" + cat, "p1 pow3" + cat, p1_init, -10., 0.)
#         self.p2 = ROOT.RooRealVar("pow3p2_" + cat, "p2 pow3" + cat, p2_init, -30., -5.)
#         self.p3 = ROOT.RooRealVar("pow3p3_" + cat, "p3 pow3" + cat, p3_init, -100., -7.)
#         self.f1 = ROOT.RooRealVar("pow3f1_" + cat, "f1 pow3" + cat, f_init, 0., 1.)
#         self.f2 = ROOT.RooRealVar("pow3f2_" + cat, "f2 pow3" + cat, f_init, 0., 1.)
#         self.sigma = ROOT.RooRealVar("sigma_pow3_" + cat,"sigma_pow3_"+cat, sigma_init,  1., 15.)
#         self.gauss = ROOT.RooGaussian("gaussxpow3_"+cat, "gaussian PDF pow3 " + cat, x, gauss_mu, self.sigma)
#         # self.step = ROOT.RooGenericPdf("step_pow3_" + cat, "step_pow3" + cat, "( ((@0-@1)*153.85 <0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85/4) ) )*( @0^@2 + @5 * @0^@3 + @6 * @0^@4)", ROOT.RooArgList(x,self.t,self.p1, self.p2, self.p3, self.f1, self.f2))
#         self.step = ROOT.RooGenericPdf("step_pow3_" + cat, "step_pow3" + cat, "( ((@0-@1)*153.85 <0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85/4) ) )*( (1-@5-@6)/((200^(1+@2) - @1^(1+@2))/(1+@2)) *@0^@2 + @5/ ((200^(1+@3) - @1^(1+@3))/(1+@3)) * @0^@3 + @6/((200^(1+@4) - @1^(1+@4))/(1+@4)) * @0^@4)", ROOT.RooArgList(x,self.t,self.p1, self.p2, self.p3, self.f1, self.f2))
#         # self.step = ROOT.RooGenericPdf("step_pow3_" + cat, "step_pow3" + cat, "( ((@0-@1)*153.85 <0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85/4) ) )*( ((@0-67.5)/130)^@2 + @5 * ((@0-67.5)/130)^@3 + @6 * ((@0-67.5)/130)^@4)", ROOT.RooArgList(x,self.t,self.p1, self.p2, self.p3, self.f1, self.f2))
#         x.setBins(20000, "cache")
#         self.pdf = ROOT.RooFFTConvPdf("pow3_" + cat + "_model", "step pow3 (X) gauss " + cat, x, self.step, self.gauss)
#         self.pdf.setBufferFraction(0.5)
#         self.name = "pow3_"+ cat
#     def checkBond(self):
#         tol = 0.0001
#         par_list = [self.t, self.p1, self.p2, self.p3, self.f1, self.f2, self.sigma]
#         if any(bondComp(par, tol) for par in par_list):
#             print ("The pdf ", self.pdf.GetName(), " needs refit.")
#     def reset(self):
#         self.t.setVal(self.init_list[1])
#         self.p1.setVal(self.init_list[2])
#         self.p2.setVal(self.init_list[3])
#         self.p3.setVal(self.init_list[4])
#         self.f1.setVal(self.init_list[5])
#         self.f2.setVal(self.init_list[5])
#         self.sigma.setVal(self.init_list[0])

class Pow3Class:
    def __init__(self, x, gauss_mu, cat="", sigma_init = 7., step_init = 108., p1_init = -5.7, p1_low = -8., p1_high = 0., p2_init = -8., p2_low = -20., p2_high = -5., p3_init = -10., p3_low = -40., p3_high = -5., f1_init = 0.9, f2_init = 0.2, f3_init = 0.3, xmax = 165., const_f1 = False, di_gauss = False):
        self.init_list = [sigma_init, step_init, p1_init, p2_init, p3_init, f1_init, f2_init,  f3_init]
        self.xmax = ROOT.RooRealVar("xmax_pow3_" + cat, "xmax_pow3_" + cat, xmax)
        self.t = ROOT.RooRealVar("pow3t_" + cat, "t pow3" + cat, step_init, 100., 118.)
        self.p1 = ROOT.RooRealVar("pow3p1_" + cat, "p1 pow3" + cat, p1_init, p1_low, p1_high)
        self.p2 = ROOT.RooRealVar("pow3p2_" + cat, "p2 pow3" + cat, p2_init, p2_low, p2_high)
        self.p3 = ROOT.RooRealVar("pow3p3_" + cat, "p3 pow3" + cat, p3_init, p3_low, p3_high)
        self.f1 = ROOT.RooRealVar("pow3f1_" + cat, "f1 pow3" + cat, f1_init, 0., 10.)
        self.f2 = ROOT.RooRealVar("pow3f2_" + cat, "f2 pow3" + cat, f2_init, 0., 10.)
        self.f3 = ROOT.RooRealVar("pow3f3_" + cat, "f3 pow3" + cat, f3_init, 0., 10.)
        self.offset = ROOT.RooRealVar("os_" + cat, "os_" + cat, 1e-10)
        if const_f1 == True:
            self.f1.setConstant(True)
        self.sigma = ROOT.RooRealVar("sigma_pow3_" + cat,"sigma_pow3_"+cat, sigma_init,  1., 15.)
        self.sigma2 = ROOT.RooRealVar("sigma2_pow3_" + cat,"sigma2_pow3_"+cat, sigma_init * 2,  0.1, 15.)
        self.gauss = ROOT.RooGaussian("gaussxpow3_"+cat, "gaussian PDF pow3 " + cat, x, gauss_mu, self.sigma)
        self.gauss2 = ROOT.RooGaussian("gauss2xpow3_"+cat, "gaussian2 PDF pow3 " + cat, x, gauss_mu, self.sigma2)
        self.gc = ROOT.RooRealVar("gc_pow3_" + cat, "gc_pow3_" + cat, 1, 0, 10)
        self.addG = ROOT.RooAddPdf("addG_pow3_" + cat, "addG_pow3_" + cat, self.gauss2, self.gauss, self.gc)
        #self.step1 = ROOT.RooGenericPdf("step1_pow3_" + cat, "step1_pow3" + cat, "(@0>@1)*(@0^@2)", ROOT.RooArgList(x,self.t,self.p1))
        #self.step2 = ROOT.RooGenericPdf("step2_pow3_" + cat, "step2_pow3" + cat, "(@0>@1)*(@0^@2)", ROOT.RooArgList(x,self.t,self.p2))
        #self.step3 = ROOT.RooGenericPdf("step3_pow3_" + cat, "step3_pow3" + cat, "(@0>@1)*(@0^@2)", ROOT.RooArgList(x,self.t,self.p3))
        #self.step = ROOT.RooAddPdf("step_pow3_"+ cat, "step_pow3_" + cat, ROOT.RooArgList(self.step1, self.step2, self.step3), ROOT.RooArgList(self.f1, self.f2))
        #self.transX = ROOT.RooFormulaVar("transX_"+cat, "@0 - 50", ROOT.RooArgList(x))
        #self.transt = ROOT.RooFormulaVar("transt_"+cat, "@0 - 50", ROOT.RooArgList(self.t))
        self.norm1 = ROOT.RooFormulaVar("norm1_pow3_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.t, self.p1, self.xmax))
        self.norm2 = ROOT.RooFormulaVar("norm2_pow3_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.t, self.p2, self.xmax))
        self.norm3 = ROOT.RooFormulaVar("norm3_pow3_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.t, self.p3, self.xmax))
        self.step = ROOT.RooGenericPdf("step_pow3_" + cat, "step_pow3" + cat, "( ((@0-@1)*307.69 <0.0) ? 0.0 : (((@0-@1)*307.69 >1.0) ? 1.0 : ((@0-@1)*307.69) ) )*( @5/ @8 *@0^@2 + @6/ @9 * @0^@3 + @7/ @10 * @0^@4)", ROOT.RooArgList(x,self.t,self.p1, self.p2, self.p3, self.f1, self.f2, self.f3, self.norm1, self.norm2, self.norm3))
        x.setBins(40000, "cache")
        if di_gauss == True:
            self.pdf = ROOT.RooFFTConvPdf("pow3_" + cat + "_model", "step pow3 (X) gauss " + cat, x, self.step, self.addG)
        else:
            self.pdf = ROOT.RooFFTConvPdf("pow3_" + cat + "_model", "step pow3 (X) gauss " + cat, x, self.step, self.gauss)
        self.pdf.setBufferFraction(0.5)
        self.SBpdf = ROOT.RooGenericPdf("SBpow3_" +cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
        self.name = "pow3_"+ cat
    def checkBond(self):
        tol = 0.0001
        par_list = [self.t, self.p1, self.p2, self.p3, self.f1, self.f2, self.f3, self.sigma]
        if any(bondComp(par, tol) for par in par_list):
            print ("The pdf ", self.pdf.GetName(), " needs refit.")
    def reset(self):
        self.t.setVal(self.init_list[1])
        self.p1.setVal(self.init_list[2])
        self.p2.setVal(self.init_list[3])
        self.p3.setVal(self.init_list[4])
        self.f1.setVal(self.init_list[5])
        self.f2.setVal(self.init_list[6])
        self.f3.setVal(self.init_list[7])
        self.sigma.setVal(self.init_list[0])

#### Exp series
class Exp1Class:
    def __init__(self, x, gauss_mu, cat="", sigma_init = 7., step_init = 108., p_init = -0.02, p_low = -2, p_high = 0., di_gauss = False):
        self.init_list = [sigma_init, step_init, p_init]
        self.sigma = ROOT.RooRealVar("sigma_exp1_" + cat,"sigma_exp1" + cat       ,sigma_init,  1., 15.)
        self.sigma2 = ROOT.RooRealVar("sigma2_exp1_" + cat,"sigma2_exp1_"+cat, sigma_init * 2,  0.1, 15.)
        self.t = ROOT.RooRealVar("exp_t1_" + cat, "t exp1 " + cat, step_init, 95., 118.)
        self.p1 = ROOT.RooRealVar("exp_p1_" + cat, "p1 exp1 " +cat, p_init, p_low, p_high)
        self.step = ROOT.RooGenericPdf("step_exp1_" + cat, "step_exp1 " + cat, "( ((@0-@1)*153.85<0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85) ) )*(TMath::Exp(@0*@2))", ROOT.RooArgList(x,self.t,self.p1))
        self.gauss = ROOT.RooGaussian("gaussxexp1_" + cat, "gaussian PDF exp1 " + cat, x, gauss_mu, self.sigma)
        self.gauss2 = ROOT.RooGaussian("gauss2xexp1_"+cat, "gaussian2 PDF exp1 " + cat, x, gauss_mu, self.sigma2)
        self.gc = ROOT.RooRealVar("gc_exp1_" + cat, "gc_exp1_" + cat, 1, 0, 10)
        self.addG = ROOT.RooAddPdf("addG_exp1_" + cat, "addG_exp1_" + cat, self.gauss2, self.gauss, self.gc)
        x.setBins(20000, "cache")
        if di_gauss == True:
            self.pdf = ROOT.RooFFTConvPdf("exp1_"+cat + "_model", "step exp1 (X) gauss" + cat, x, self.step, self.addG)
        else:
            self.pdf = ROOT.RooFFTConvPdf("exp1_"+cat + "_model", "step exp1 (X) gauss" + cat, x, self.step, self.gauss)
        self.pdf.setBufferFraction(0.5)
        self.SBpdf = ROOT.RooGenericPdf("SBexp1_" +cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
        self.name = "exp1_"+ cat
    def checkBond(self):
        tol = 0.001
        par_list = [self.t, self.p1, self.sigma]
        if any(bondComp(par, tol) for par in par_list):
            print ("The pdf ", self.pdf.GetName(), " needs refit.")
    def reset(self):
        self.t.setVal(self.init_list[1])
        self.p1.setVal(self.init_list[2])
        self.sigma.setVal(self.init_list[0])

class Exp2Class:
    def __init__(self, x, gauss_mu, cat="", sigma_init = 7., step_init = 105., p1_init = -0.02, p1_low = -1, p1_high = 0., p2_init = -0.02, p2_low = -1, p2_high = 0., f1_init = 0.5, f2_init = 0.5, xmax = 165., const_f1 = False, di_gauss = False):
        self.init_list = [sigma_init, step_init, p1_init, p2_init, f1_init, f2_init]
        self.xmax = ROOT.RooRealVar("xmax_exp2_" + cat, "xmax_exp2_" + cat, xmax)
        self.sigma = ROOT.RooRealVar("sigma_exp2_" + cat,"sigma_exp2" + cat       ,sigma_init,  1., 15.)
        self.sigma2 = ROOT.RooRealVar("sigma2_exp2_" + cat,"sigma2_exp2_"+cat, sigma_init * 2,  0.1, 15.)
        self.offset = ROOT.RooRealVar("os_" + cat, "os_" + cat, 1e-10)
        self.t = ROOT.RooRealVar("exp2_t1_" + cat, "t exp2 " + cat, step_init, 90., 115.)
        self.p1 = ROOT.RooRealVar("exp2_p1_" + cat, "p1 exp2 " +cat, p1_init, p1_low, p1_high)
        self.p2 = ROOT.RooRealVar("exp2_p2_" + cat, "p2 exp2 " +cat, p2_init, p2_low, p2_high)
        self.f1 = ROOT.RooRealVar("exp2_f1_" + cat, "f1 exp2 "+ cat, f1_init, 0., 10.)
        self.f2 = ROOT.RooRealVar("exp2_f2_" + cat, "f2 exp2 "+ cat, f2_init, 0., 10.)
        if const_f1 == True:
            self.f1.setConstant(True)
        self.norm1 = ROOT.RooFormulaVar("norm1_exp2_"+cat, "(TMath::Exp(@2*@1) - TMath::Exp(@0*@1))/@1", ROOT.RooArgList(self.t, self.p1, self.xmax))
        self.norm2 = ROOT.RooFormulaVar("norm2_exp2_"+cat, "(TMath::Exp(@2*@1) - TMath::Exp(@0*@1))/@1", ROOT.RooArgList(self.t, self.p2, self.xmax))
        self.step = ROOT.RooGenericPdf("step_exp2_" + cat, "step_exp2_" + cat, "( ((@0-@1)*153.85<0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85) ) )*(@5 / @6 * TMath::Exp(@0*@2) + @4 / @7 * TMath::Exp(@0*@3))", ROOT.RooArgList(x,self.t,self.p1,self.p2,self.f1, self.f2, self.norm1, self.norm2))
        self.gauss = ROOT.RooGaussian("gaussxexp2_" + cat, "gaussian PDF exp2 " + cat, x, gauss_mu, self.sigma)
        self.gauss2 = ROOT.RooGaussian("gauss2xexp2_"+cat, "gaussian2 PDF exp2 " + cat, x, gauss_mu, self.sigma2)
        self.gc = ROOT.RooRealVar("gc_exp2_" + cat, "gc_exp2_" + cat, 1, 0, 10)
        self.addG = ROOT.RooAddPdf("addG_exp2_" + cat, "addG_exp2_" + cat, self.gauss2, self.gauss, self.gc)
        x.setBins(20000, "cache")
        if di_gauss == True:
            self.pdf = ROOT.RooFFTConvPdf("exp2_"+cat + "_model", "step exp2 (X) gauss" + cat, x, self.step, self.addG)
        else:
            self.pdf = ROOT.RooFFTConvPdf("exp2_"+cat + "_model", "step exp2 (X) gauss" + cat, x, self.step, self.gauss)
        self.pdf.setBufferFraction(0.5)
        self.SBpdf = ROOT.RooGenericPdf("SBexp2_" +cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
        self.name = "exp2_"+ cat
    def checkBond(self):
        tol = 0.001
        par_list = [self.t, self.p1, self.p2, self.f1, self.f2, self.sigma]
        if any(bondComp(par, tol) for par in par_list):
            print ("The pdf ", self.pdf.GetName(), " needs refit.")
    def reset(self):
        self.t.setVal(self.init_list[1])
        self.p1.setVal(self.init_list[2])
        self.p2.setVal(self.init_list[3])
        self.f1.setVal(self.init_list[4])
        self.f2.setVal(self.init_list[5])
        self.sigma.setVal(self.init_list[0])

class Exp3Class:
    def __init__(self, x, gauss_mu, cat="", sigma_init = 7., step_init = 108., p1_init = -0.02, p1_low = -2, p1_high = 0., p2_init = -0.02, p2_low = -2, p2_high = 0., p3_init = -0.02, p3_low = -2, p3_high = 0., f1_init = 0.33, f2_init = 0.5, f3_init = 0.5, xmax = 165., const_f1 = False, di_gauss = False):
        self.init_list = [sigma_init, step_init, p1_init, p2_init, p3_init, f1_init, f2_init, f3_init]
        self.sigma = ROOT.RooRealVar("sigma_exp3_" + cat,"sigma_exp3" + cat       ,sigma_init,  1., 15.)
        self.sigma2 = ROOT.RooRealVar("sigma2_exp3_" + cat,"sigma2_exp3_"+cat, sigma_init * 2,  0.1, 15.)
        self.xmax = ROOT.RooRealVar("xmax_exp3_" + cat, "xmax_exp3_" + cat, xmax)
        self.t = ROOT.RooRealVar("exp3_t1_" + cat, "t exp3 " + cat, step_init, 100., 115.)
        self.p1 = ROOT.RooRealVar("exp3_p1_" + cat, "p1 exp3 " +cat, p1_init, -1., 0)
        self.p2 = ROOT.RooRealVar("exp3_p2_" + cat, "p2 exp3 " +cat, p2_init, -1., 0)
        self.p3 = ROOT.RooRealVar("exp3_p3_" + cat, "p3 exp3 " +cat, p3_init, -2., 0)
        self.f1 = ROOT.RooRealVar("exp3_f1_" + cat, "f1 exp3 "+ cat, f1_init, 0, 10.)
        self.f2 = ROOT.RooRealVar("exp3_f2_" + cat, "f2 exp3 "+ cat, f2_init, 0, 10.)
        self.f3 = ROOT.RooRealVar("exp3_f3_" + cat, "f3 exp3 "+ cat, f3_init, 0., 10.)
        if const_f1 == True:
            self.f1.setConstant(True)
        self.norm1 = ROOT.RooFormulaVar("norm1_exp3_"+cat, "(TMath::Exp(@2*@1) - TMath::Exp(@0*@1))/@1", ROOT.RooArgList(self.t, self.p1, self.xmax))
        self.norm2 = ROOT.RooFormulaVar("norm2_exp3_"+cat, "(TMath::Exp(@2*@1) - TMath::Exp(@0*@1))/@1", ROOT.RooArgList(self.t, self.p2, self.xmax))
        self.norm3 = ROOT.RooFormulaVar("norm3_exp3_"+cat, "(TMath::Exp(@2*@1) - TMath::Exp(@0*@1))/@1", ROOT.RooArgList(self.t, self.p3, self.xmax))
        self.step = ROOT.RooGenericPdf("step_exp3_" + cat, "step_exp3_" + cat, "( ((@0-@1)*153.85<0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85) ) )*( @5 / @8 * TMath::Exp(@0*@2) + @6 / @9 * TMath::Exp(@0*@3) + @7 / @10  * TMath::Exp(@0*@4))", ROOT.RooArgList(x,self.t,self.p1,self.p2, self.p3, self.f1, self.f2, self.f3, self.norm1, self.norm2, self.norm3))
        self.gauss = ROOT.RooGaussian("gaussxexp3_" + cat, "gaussian PDF exp3 " + cat, x, gauss_mu, self.sigma)
        self.gauss2 = ROOT.RooGaussian("gauss2xexp3_"+cat, "gaussian2 PDF exp3 " + cat, x, gauss_mu, self.sigma2)
        self.gc = ROOT.RooRealVar("gc_exp3_" + cat, "gc_exp3_" + cat, 1, 0, 10)
        self.addG = ROOT.RooAddPdf("addG_exp3_" + cat, "addG_exp3_" + cat, self.gauss2, self.gauss, self.gc)
        x.setBins(20000, "cache")
        if di_gauss == True:
            self.pdf = ROOT.RooFFTConvPdf("exp3_"+cat + "_model", "step exp3 (X) gauss" + cat, x, self.step, self.addG)
        else:
            self.pdf = ROOT.RooFFTConvPdf("exp3_"+cat + "_model", "step exp3 (X) gauss" + cat, x, self.step, self.gauss)
        self.pdf.setBufferFraction(0.5)
        self.SBpdf = ROOT.RooGenericPdf("SBexp3_" +cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
        self.name = "exp3_"+ cat
    def checkBond(self):
        tol = 0.0001
        par_list = [self.t, self.p1, self.p2, self.p3, self.f1, self.f2, self.f3, self.sigma]
        if any(bondComp(par, tol) for par in par_list):
            print ("The pdf ", self.pdf.GetName(), " needs refit.")
    def reset(self):
        self.t.setVal(self.init_list[1])
        self.p1.setVal(self.init_list[2])
        self.p2.setVal(self.init_list[3])
        self.p3.setVal(self.init_list[4])
        self.f1.setVal(self.init_list[5])
        self.f2.setVal(self.init_list[6])
        self.f3.setVal(self.init_list[7])
        self.sigma.setVal(self.init_list[0])

#### Laurent series
class Lau2Class:
    def __init__(self, x, gauss_mu, cat="", sigma_init = 7., step_init = 108., p1 = -7, p2 = -6, f_init = 0.1, xmax = 165., const_f1 = False, di_gauss = False):
        self.init_list = [sigma_init, step_init, p1, p2, f_init]
        self.sigma = ROOT.RooRealVar("sigma_lau2_" + cat,"sigma_lau2_" + cat       ,sigma_init,  1., 15.)
        self.sigma2 = ROOT.RooRealVar("sigma2_lau2_" + cat,"sigma2_lau2_"+cat, sigma_init * 2,  0.1, 15.)
        self.offset = ROOT.RooRealVar("os_" + cat, "os_" + cat, 1e-10)
        self.xmax = ROOT.RooRealVar("xmax_lau2_" + cat, "xmax_lau2_" + cat, xmax)
        self.t = ROOT.RooRealVar("lau2_t_" + cat, "t lau2 " + cat, step_init, 95., 118.)
        self.f1 = ROOT.RooRealVar("lau2_f1_" + cat, "f1 lau2 " + cat, f_init, 0, 10.)
        self.f2 = ROOT.RooRealVar("lau2_f2_" + cat, "f2 lau2 " + cat, f_init, 0., 10.)
        self.p1 = ROOT.RooRealVar("lau2_p1_" + cat, "p1 lau2 " + cat, p1)
        self.p2 = ROOT.RooRealVar("lau2_p2_" + cat, "p2 lau2 " + cat, p2)
        self.norm1 = ROOT.RooFormulaVar("norm1_lau2_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.t, self.p1, self.xmax))
        self.norm2 = ROOT.RooFormulaVar("norm2_lau2_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.t, self.p2, self.xmax))
        self.step = ROOT.RooGenericPdf("step_lau2_" + cat, "step_lau2_" + cat, "( ((@0-@1)*153.85<0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85) ) )*(@2/@6*@0^(@4) + @3/@7*@0^(@5))", \
                                        ROOT.RooArgList(x,self.t,self.f1, self.f2, self.p1, self.p2, self.norm1, self.norm2))
        self.gauss = ROOT.RooGaussian("gaussxlau2_" + cat, "gaussian PDF lau2 " + cat, x, gauss_mu, self.sigma)
        self.gauss2 = ROOT.RooGaussian("gauss2xlau2_"+cat, "gaussian2 PDF lau2 " + cat, x, gauss_mu, self.sigma2)
        self.gc = ROOT.RooRealVar("gc_lau2_" + cat, "gc_lau2_" + cat, 1, 0, 10)
        self.addG = ROOT.RooAddPdf("addG_lau2_" + cat, "addG_lau2_" + cat, self.gauss2, self.gauss, self.gc)
        if const_f1 == True:
            self.f1.setConstant(True)
        x.setBins(20000, "cache")
        if di_gauss == True:
            self.pdf = ROOT.RooFFTConvPdf("lau2_" +cat+ "_model", "step lau2 (X) gauss " + cat, x, self.step, self.addG)
        else:
            self.pdf = ROOT.RooFFTConvPdf("lau2_" +cat+ "_model", "step lau2 (X) gauss " + cat, x, self.step, self.gauss)
        self.pdf.setBufferFraction(0.5)
        self.SBpdf = ROOT.RooGenericPdf("SBlau2_" +cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
        self.name = "lau2_"+ cat
    def checkBond(self):
        tol = 0.0001
        par_list = [self.t, self.p1, self.p2, self.f1, self.f2, self.sigma]
        if any(bondComp(par, tol) for par in par_list):
            print ("The pdf ", self.pdf.GetName(), " needs refit.")
    def reset(self):
        self.t.setVal(self.init_list[1])
        self.p1.setVal(self.init_list[2])
        self.p2.setVal(self.init_list[3])
        self.f1.setVal(self.init_list[4])
        self.f2.setVal(self.init_list[4])
        self.sigma.setVal(self.init_list[0])

class Lau3Class:
    def __init__(self, x, gauss_mu, cat="", sigma_init = 7., step_init = 108., p1 = -8, p2 = -7, p3 = -6, f_init = 0.1, xmax = 165., const_f1 = False, di_gauss = False):
        self.init_list = [sigma_init, step_init, p1, p2, p3, f_init]
        self.sigma = ROOT.RooRealVar("sigma_lau3_" + cat,"sigma_lau3_" + cat       ,sigma_init,  1., 15.)
        self.sigma2 = ROOT.RooRealVar("sigma2_lau3_" + cat,"sigma2_lau3_"+cat, sigma_init * 2,  0.1, 15.)
        self.t = ROOT.RooRealVar("lau3_t_" + cat, "t lau3 " + cat, step_init, 95., 118.)
        self.xmin = ROOT.RooRealVar("lau3_xmin_" + cat, "xmin lau3 " + cat,  100)
        self.xmax = ROOT.RooRealVar("xmax_lau3_" + cat, "xmax_lau3_" + cat, xmax)
        self.offset = ROOT.RooRealVar("os_" + cat, "os_" + cat, 1e-5, 0., 0.01)
        self.f1 = ROOT.RooRealVar("lau3_f1_" + cat, "f1 lau3 " + cat, f_init, 0., 10.)
        self.f2 = ROOT.RooRealVar("lau3_f2_" + cat, "f2 lau3 " + cat, f_init, 0., 10.)
        self.f3 = ROOT.RooRealVar("lau3_f3_" + cat, "f3 lau3 " + cat, f_init, 0., 10.)
        self.p1 = ROOT.RooRealVar("lau3_p1_" + cat, "p1 lau3 " + cat, p1)
        self.p2 = ROOT.RooRealVar("lau3_p2_" + cat, "p2 lau3 " + cat, p2)
        self.p3 = ROOT.RooRealVar("lau3_p3_" + cat, "p3 lau3 " + cat, p3)
        self.norm1 = ROOT.RooFormulaVar("norm1_lau3_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.xmin, self.p1, self.xmax))
        self.norm2 = ROOT.RooFormulaVar("norm2_lau3_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.xmin, self.p2, self.xmax))
        self.norm3 = ROOT.RooFormulaVar("norm3_lau3_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.xmin, self.p3, self.xmax))
        self.step = ROOT.RooGenericPdf("step_lau3_" + cat, "step_lau3_" + cat, "( ((@0-@1)*153.85<0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85) ) )*(@2/@8*@0^(@4) + @3/@9*@0^(@5) + @7/@10 * @0^(@6))", \
                                    ROOT.RooArgList(x,self.t, self.f1, self.f2, self.p1, self.p2, self.p3, self.f3, self.norm1, self.norm2, self.norm3))
        self.gauss = ROOT.RooGaussian("gaussxlau3_" + cat, "gaussian PDF lau3 " + cat, x, gauss_mu, self.sigma)
        self.gauss2 = ROOT.RooGaussian("gauss2xlau3_"+cat, "gaussian2 PDF lau3 " + cat, x, gauss_mu, self.sigma2)
        self.gc = ROOT.RooRealVar("gc_lau3_" + cat, "gc_lau3_" + cat, 1, 0, 10)
        self.addG = ROOT.RooAddPdf("addG_lau3_" + cat, "addG_lau3_" + cat, self.gauss2, self.gauss, self.gc)
        if const_f1 == True:
            self.f1.setConstant(True)
        x.setBins(20000, "cache")
        if di_gauss == True:
            self.pdf = ROOT.RooFFTConvPdf("lau3_" +cat+ "_model", "step lau3 (X) gauss " + cat, x, self.step, self.addG)
        else:
            self.pdf = ROOT.RooFFTConvPdf("lau3_" +cat+ "_model", "step lau3 (X) gauss " + cat, x, self.step, self.gauss)
        self.pdf.setBufferFraction(0.5)
        self.SBpdf = ROOT.RooGenericPdf("SBlau3_" +cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
        self.name = "lau3_"+ cat
    def checkBond(self):
        tol = 0.0001
        par_list = [self.t, self.p1, self.p2, self.p3, self.f1, self.f2, self.sigma]
        if any(bondComp(par, tol) for par in par_list):
            print ("The pdf ", self.pdf.GetName(), " needs refit.")
    def reset(self):
        self.t.setVal(self.init_list[1])
        self.p1.setVal(self.init_list[2])
        self.p2.setVal(self.init_list[3])
        self.p3.setVal(self.init_list[4])
        self.f1.setVal(self.init_list[5])
        self.f2.setVal(self.init_list[5])
        self.f3.setVal(self.init_list[5])
        self.sigma.setVal(self.init_list[0])
    
class Lau4Class:
    def __init__(self, x, gauss_mu, cat="", sigma_init = 7., step_init = 108., p1 = -8, p2 = -7, p3 = -6, p4 = -5, f_init = 0.1, xmax = 165., const_f1 = False, di_gauss = False):
        self.init_list = [sigma_init, step_init, p1, p2, p3, p4, f_init]
        self.sigma = ROOT.RooRealVar("sigma_lau4_" + cat,"sigma_lau4_" + cat       ,sigma_init,  1., 15.)
        self.sigma2 = ROOT.RooRealVar("sigma2_lau4_" + cat,"sigma2_lau4_"+cat, sigma_init * 2,  0.1, 15.)
        self.t = ROOT.RooRealVar("lau4_t_" + cat, "t lau4 " + cat, step_init, 95., 118.)
        self.xmax = ROOT.RooRealVar("xmax_lau4_" + cat, "xmax_lau4_" + cat, xmax)
        self.offset = ROOT.RooRealVar("os_" + cat, "os_" + cat, 1e-5, 0., 0.01)
        self.f1 = ROOT.RooRealVar("lau4_f1_" + cat, "f1 lau4 " + cat, f_init, 0., 10.)
        self.f2 = ROOT.RooRealVar("lau4_f2_" + cat, "f2 lau4 " + cat, f_init, 0., 10.)
        self.f3 = ROOT.RooRealVar("lau4_f3_" + cat, "f3 lau4 " + cat, f_init, 0., 10.)
        self.f4 = ROOT.RooRealVar("lau4_f4_" + cat, "f4 lau4 " + cat, f_init, 0., 10.)
        self.p1 = ROOT.RooRealVar("lau4_p1_" + cat, "p1 lau4 " + cat, p1)
        self.p2 = ROOT.RooRealVar("lau4_p2_" + cat, "p2 lau4 " + cat, p2)
        self.p3 = ROOT.RooRealVar("lau4_p3_" + cat, "p3 lau4 " + cat, p3)
        self.p4 = ROOT.RooRealVar("lau4_p4_" + cat, "p4 lau4 " + cat, p4)
        self.norm1 = ROOT.RooFormulaVar("norm1_lau4_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.t, self.p1, self.xmax))
        self.norm2 = ROOT.RooFormulaVar("norm2_lau4_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.t, self.p2, self.xmax))
        self.norm3 = ROOT.RooFormulaVar("norm3_lau4_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.t, self.p3, self.xmax))
        self.norm4 = ROOT.RooFormulaVar("norm4_lau4_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.t, self.p4, self.xmax))

        self.step = ROOT.RooGenericPdf("step_lau4_" + cat, "step_lau4_" + cat, "( ((@0-@1)*153.85<0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85) ) )*(@2/@10*@0^(@5) + @3/@11*@0^(@6) + @4/@12*@0^(@7) + @9/@13*@0^(@8))", \
                                      ROOT.RooArgList(x,self.t, self.f1, self.f2, self.f3, self.p1, self.p2, self.p3, self.p4, self.f4, self.norm1, self.norm2, self.norm3, self.norm4))
        self.gauss = ROOT.RooGaussian("gaussxlau4_" + cat, "gaussian PDF lau4 " + cat, x, gauss_mu, self.sigma)
        self.gauss2 = ROOT.RooGaussian("gauss2xlau4_"+cat, "gaussian2 PDF lau4 " + cat, x, gauss_mu, self.sigma2)
        self.gc = ROOT.RooRealVar("gc_lau4_" + cat, "gc_lau4_" + cat, 1, 0, 10)
        self.addG = ROOT.RooAddPdf("addG_lau4_" + cat, "addG_lau4_" + cat, self.gauss2, self.gauss, self.gc)
        if const_f1 == True:
            self.f1.setConstant(True)
        x.setBins(20000, "cache")
        if di_gauss == True:
            self.pdf = ROOT.RooFFTConvPdf("lau4_" +cat+ "_model", "step lau4 (X) gauss " + cat, x, self.step, self.addG)
        else:
            self.pdf = ROOT.RooFFTConvPdf("lau4_" +cat+ "_model", "step lau4 (X) gauss " + cat, x, self.step, self.gauss)
        self.pdf.setBufferFraction(0.5)
        self.SBpdf = ROOT.RooGenericPdf("SBlau4_" +cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
        self.name = "lau4_"+ cat
    def checkBond(self):
        tol = 0.0001
        par_list = [self.t, self.p1, self.p2, self.p3, self.p4, self.f1, self.f2, self.f3, self.sigma]
        if any(bondComp(par, tol) for par in par_list):
            print ("The pdf ", self.pdf.GetName(), " needs refit.")

    def reset(self):
        self.t.setVal(self.init_list[1])
        self.p1.setVal(self.init_list[2])
        self.p2.setVal(self.init_list[3])
        self.p3.setVal(self.init_list[4])
        self.p4.setVal(self.init_list[5])
        self.f1.setVal(self.init_list[6])
        self.f2.setVal(self.init_list[6])
        self.f2.setVal(self.init_list[6])
        self.f3.setVal(self.init_list[6])
        self.f4.setVal(self.init_list[6])
        self.sigma.setVal(self.init_list[0])

class Lau5Class:
    def __init__(self, x, gauss_mu, cat="", sigma_init = 7., step_init = 108., p1 = -8, p2 = -7, p3 = -6, p4 = -5, p5 = -4, f_init = 0.1, xmax = 165., const_f1 = False):
        self.init_list = [sigma_init, step_init, p1, p2, p3, p4, p5, f_init]
        self.sigma = ROOT.RooRealVar("sigma_lau5_" + cat,"sigma_lau5_" + cat       ,sigma_init,  1., 15.)
        self.xmax = ROOT.RooRealVar("xmax_lau5_" + cat, "xmax_lau5_" + cat, xmax)
        self.t = ROOT.RooRealVar("lau5_t_" + cat, "t lau5 " + cat, step_init, 95., 118.)
        self.offset = ROOT.RooRealVar("os_" + cat, "os_" + cat, 1e-5, 0., 0.01)
        self.f1 = ROOT.RooRealVar("lau5_f1_" + cat, "f1 lau5 " + cat, f_init, 0., 1.)
        self.f2 = ROOT.RooRealVar("lau5_f2_" + cat, "f2 lau5 " + cat, f_init, 0., 1.)
        self.f3 = ROOT.RooRealVar("lau5_f3_" + cat, "f3 lau5 " + cat, f_init, 0., 1.)
        self.f4 = ROOT.RooRealVar("lau5_f4_" + cat, "f4 lau5 " + cat, f_init, 0., 1.)
        self.f5 = ROOT.RooRealVar("lau5_f5_" + cat, "f5 lau5 " + cat, f_init, 0., 1.)
        self.p1 = ROOT.RooRealVar("lau5_p1_" + cat, "p1 lau5 " + cat, p1)
        self.p2 = ROOT.RooRealVar("lau5_p2_" + cat, "p2 lau5 " + cat, p2)
        self.p3 = ROOT.RooRealVar("lau5_p3_" + cat, "p3 lau5 " + cat, p3)
        self.p4 = ROOT.RooRealVar("lau5_p4_" + cat, "p4 lau5 " + cat, p4)
        self.p5 = ROOT.RooRealVar("lau5_p5_" + cat, "p5 lau5 " + cat, p5)
        self.norm1 = ROOT.RooFormulaVar("norm1_lau5_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.t, self.p1, self.xmax))
        self.norm2 = ROOT.RooFormulaVar("norm2_lau5_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.t, self.p2, self.xmax))
        self.norm3 = ROOT.RooFormulaVar("norm3_lau5_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.t, self.p3, self.xmax))
        self.norm4 = ROOT.RooFormulaVar("norm4_lau5_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.t, self.p4, self.xmax))
        self.norm5 = ROOT.RooFormulaVar("norm5_lau5_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.t, self.p5, self.xmax))
        
        self.step = ROOT.RooGenericPdf("step_lau5_" + cat, "step_lau5_" + cat, "( ((@0-@1)*153.85<0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85) ) )*(@2/@12*@0^@7 + @3/@13*@0^@8 + @4/@14*@0^@9 + @5/@15*@0^@10 + @6/@16*@0^@11)", \
                                      ROOT.RooArgList(x,self.t, self.f1, self.f2, self.f3, self.f4, self.f5, self.p1, self.p2, self.p3, self.p4, self.p5, self.norm1, self.norm2, self.norm3, self.norm4, self.norm5))
        self.gauss = ROOT.RooGaussian("gaussxlau5_" + cat, "gaussian PDF lau5 " + cat, x, gauss_mu, self.sigma)
        if const_f1 == True:
            self.f1.setConstant(True)
        x.setBins(20000, "cache")
        self.pdf = ROOT.RooFFTConvPdf("lau5_" +cat+ "_model", "step lau5 (X) gauss " + cat, x, self.step, self.gauss)
        self.pdf.setBufferFraction(0.5)
        self.SBpdf = ROOT.RooGenericPdf("SBlau5_" +cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
        self.name = "lau5_"+ cat
    def checkBond(self):
        tol = 0.0001
        par_list = [self.t, self.p1, self.p2, self.p3, self.p4, self.p5, self.f1, self.f2, self.f3, self.f4, self.f5, self.sigma]
        if any(bondComp(par, tol) for par in par_list):
            print ("The pdf ", self.pdf.GetName(), " needs refit.")

    def reset(self):
        self.t.setVal(self.init_list[1])
        self.p1.setVal(self.init_list[2])
        self.p2.setVal(self.init_list[3])
        self.p3.setVal(self.init_list[4])
        self.p4.setVal(self.init_list[5])
        self.p5.setVal(self.init_list[6])
        self.f1.setVal(self.init_list[7])
        self.f2.setVal(self.init_list[7])
        self.f2.setVal(self.init_list[7])
        self.f3.setVal(self.init_list[7])
        self.f4.setVal(self.init_list[7])
        self.f5.setVal(self.init_list[7])
        self.sigma.setVal(self.init_list[0])


class ModGausClass:
    def __init__(self, x, cat="", lowx = 105., highx = 170., m0 = 105, sl = 1, sh = 40, vl = 2, vr = 1):
        self.m0 = ROOT.RooRealVar("m0_"+cat,"mass peak value [GeV]" , m0,100,130)
        self.vl = ROOT.RooRealVar("nuL_" +cat,"low-end power"       , vl,  0,  15)
        self.vr = ROOT.RooRealVar("nuRange_"+cat,"power range"       , vr, -5,  5)
        self.s0 = ROOT.RooRealVar("NOUSE_sigma0_"+cat,"peak width"       , 0)#,  1., 10.)
        self.sl = ROOT.RooRealVar("sigmaL_"+cat,"low-end width"    , sl, 1., 40)
        self.sh = ROOT.RooRealVar("sigmaH_"+cat,"high-end width"   , sh, 1.,60)
        self.pdf = ROOT.ModGaus("modg_"+cat+"_model","modg_"+cat+"_model", x, self.m0, self.vl, self.vr, self.s0, self.sl, self.sh, lowx, highx)
        self.SBpdf = ROOT.RooGenericPdf("SBmodg_" +cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
        self.name = "modg_"+ cat
    def checkBond(self):
        tol = 0.001
        par_list = [self.m0, self.vl, self.vr, self.s0, self.sl, self.sh]
        if any(bondComp(par, tol) for par in par_list):
            print ("The pdf ", self.pdf.GetName(), " needs refit.")
    def reset(self):
        pass


# Core Functions
# Remember to change the mllg var name to "CMS_hzg_mass"
class CoreBern2Class:
    def __init__(self, x, cat="", p0=10, p_init=0.3, bond=20, fileName="ZGCoreShape_01jet_NAFCorr", shapeName=""):
        self.init_list = [p0, p_init]
        self.p0 = ROOT.RooRealVar("cb2p0_" + cat, "cb2p1_" + cat, p0)
        self.p1 = ROOT.RooRealVar("cb2p1_" + cat, "cb2p1_" + cat, p_init,-bond, bond)
        self.p2 = ROOT.RooRealVar("cb2p2_" + cat, "cb2p2_" + cat, p_init,-bond, bond)
        file_= ROOT.TFile("/afs/cern.ch/user/f/fanx/CMSSW_12_6_0_patch1/src/HtoZg_fitting/Core_func_shape/"+fileName+".root", "READ")
        w = file_.Get("w")
        self.ZGMCShape = w.pdf(shapeName)
        w.var(x.GetName()).setRange(x.getMin(), x.getMax())
        self.bern = ROOT.RooBernstein("bern2_core_" + cat, "bern2_core_" + cat, x, ROOT.RooArgList(self.p0, self.p1, self.p2))
        self.pdf = ROOT.RooEffProd("corebern2_" +cat + "_model", "Bernstein2 X Core " + cat, self.ZGMCShape, self.bern)
        self.SBpdf = ROOT.RooGenericPdf("SBcorebern2_" +cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
        self.name = "core_bern2_"+ cat
    def checkBond(self):
        tol = 0.001
        par_list = [self.p1, self.p2]
        if any(bondComp(par, tol) for par in par_list):
            print ("The pdf ", self.pdf.GetName(), " needs refit.")
    def reset(self):
         self.p0.setVal(self.init_list[0])
         self.p1.setVal(self.init_list[1])
         self.p2.setVal(self.init_list[1])

class CoreBern3Class:
    def __init__(self, x, cat="", p0=10, p_init=0.3, bond=20, fileName="ZGCoreShape_01jet_NAFCorr", shapeName=""):
        self.init_list = [p0, p_init]
        self.p0 = ROOT.RooRealVar("cb3p0_" + cat, "cb3p1_" + cat, p0)
        self.p1 = ROOT.RooRealVar("cb3p1_" + cat, "cb3p1_" + cat, p_init,-bond, bond)
        self.p2 = ROOT.RooRealVar("cb3p2_" + cat, "cb3p2_" + cat, p_init,-bond, bond)
        self.p3 = ROOT.RooRealVar("cb3p3_" + cat, "cb3p3_" + cat, p_init,-bond, bond)
        file_= ROOT.TFile("/afs/cern.ch/user/f/fanx/CMSSW_12_6_0_patch1/src/HtoZg_fitting/Core_func_shape/"+fileName+".root", "READ")
        w = file_.Get("w")
        self.ZGMCShape = w.pdf(shapeName)
        w.var(x.GetName()).setRange(x.getMin(), x.getMax())
        self.bern = ROOT.RooBernstein("bern3_core_" + cat, "bern3_core_" + cat, x, ROOT.RooArgList(self.p0, self.p1, self.p2, self.p3))
        self.pdf = ROOT.RooEffProd("corebern3_" +cat + "_model", "Bernstein3 X Core " + cat, self.ZGMCShape, self.bern)
        self.SBpdf = ROOT.RooGenericPdf("SBcorebern3_" +cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
        self.name = "core_bern3_"+ cat
    def checkBond(self):
        tol = 0.001
        par_list = [self.p1, self.p2, self.p3]
        if any(bondComp(par, tol) for par in par_list):
            print ("The pdf ", self.pdf.GetName(), " needs refit.")
    def reset(self):
         self.p0.setVal(self.init_list[0])
         self.p1.setVal(self.init_list[1])
         self.p2.setVal(self.init_list[1])
         self.p3.setVal(self.init_list[1])

class CoreBern4Class:
    def __init__(self, x, cat="", p0=10, p_init=0.3, bond=20, fileName="ZGCoreShape_01jet_NAFCorr", shapeName=""):
        self.init_list = [p0, p_init]
        self.p0 = ROOT.RooRealVar("cb4p0_" + cat, "cb4p1_" + cat, p0)
        self.p1 = ROOT.RooRealVar("cb4p1_" + cat, "cb4p1_" + cat, p_init,-bond, bond)
        self.p2 = ROOT.RooRealVar("cb4p2_" + cat, "cb4p2_" + cat, p_init,-bond, bond)
        self.p3 = ROOT.RooRealVar("cb4p3_" + cat, "cb4p3_" + cat, p_init,-bond, bond)
        self.p4 = ROOT.RooRealVar("cb4p4_" + cat, "cb4p4_" + cat, p_init,-bond, bond)
        file_= ROOT.TFile("/afs/cern.ch/user/f/fanx/CMSSW_12_6_0_patch1/src/HtoZg_fitting/Core_func_shape/"+fileName+".root", "READ")
        w = file_.Get("w")
        self.ZGMCShape = w.pdf(shapeName)
        w.var(x.GetName()).setRange(x.getMin(), x.getMax())
        self.bern = ROOT.RooBernstein("bern4_core_" + cat, "bern4_core_" + cat, x, ROOT.RooArgList(self.p0, self.p1, self.p2, self.p3, self.p4))
        self.pdf = ROOT.RooEffProd("corebern4_" +cat + "_model", "Bernstein4 X Core " + cat, self.ZGMCShape, self.bern)
        self.SBpdf = ROOT.RooGenericPdf("SBcorebern4_" +cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
        self.name = "core_bern4_"+ cat
    def checkBond(self):
        tol = 0.001
        par_list = [self.p1, self.p2, self.p3, self.p4]
        if any(bondComp(par, tol) for par in par_list):
            print ("The pdf ", self.pdf.GetName(), " needs refit.")
    def reset(self):
         self.p0.setVal(self.init_list[0])
         self.p1.setVal(self.init_list[1])
         self.p2.setVal(self.init_list[1])
         self.p3.setVal(self.init_list[1])
         self.p4.setVal(self.init_list[1])

class CoreBern5Class:
    def __init__(self, x, cat="", p0=10, p_init=0.3, bond=20, fileName="ZGCoreShape_01jet_NAFCorr", shapeName=""):
        self.init_list = [p0, p_init]
        self.p0 = ROOT.RooRealVar("cb5p0_" + cat, "cb5p1_" + cat, p0)
        self.p1 = ROOT.RooRealVar("cb5p1_" + cat, "cb5p1_" + cat, p_init,-bond, bond)
        self.p2 = ROOT.RooRealVar("cb5p2_" + cat, "cb5p2_" + cat, p_init,-bond, bond)
        self.p3 = ROOT.RooRealVar("cb5p3_" + cat, "cb5p3_" + cat, p_init,-bond, bond)
        self.p4 = ROOT.RooRealVar("cb5p4_" + cat, "cb5p4_" + cat, p_init,-bond, bond)
        self.p5 = ROOT.RooRealVar("cb5p5_" + cat, "cb5p5_" + cat, p_init,-bond, bond)

        file_= ROOT.TFile("/afs/cern.ch/user/f/fanx/CMSSW_12_6_0_patch1/src/HtoZg_fitting/Core_func_shape/"+fileName+".root", "READ")
        w = file_.Get("w")
        self.ZGMCShape = w.pdf(shapeName)
        w.var(x.GetName()).setRange(x.getMin(), x.getMax())
        self.bern = ROOT.RooBernstein("bern5_core_" + cat, "bern5_core_" + cat, x, ROOT.RooArgList(self.p0, self.p1, self.p2, self.p3, self.p4, self.p5))
        self.pdf = ROOT.RooEffProd("corebern5_" +cat + "_model", "Bernstein5 X Core " + cat, self.ZGMCShape, self.bern)
        self.SBpdf = ROOT.RooGenericPdf("SBcorebern5_" +cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
        self.name = "core_bern5_"+ cat
    def checkBond(self):
        tol = 0.001
        par_list = [self.p1, self.p2, self.p3, self.p4, self.p5]
        if any(bondComp(par, tol) for par in par_list):
            print ("The pdf ", self.pdf.GetName(), " needs refit.")
    def reset(self):
         self.p0.setVal(self.init_list[0])
         self.p1.setVal(self.init_list[1])
         self.p2.setVal(self.init_list[1])
         self.p3.setVal(self.init_list[1])
         self.p4.setVal(self.init_list[1])
         self.p5.setVal(self.init_list[1])
                  


class CorePow1Class:
    def __init__(self, x, cat="", p_init = -6., fileName="ZGCoreShape_01jet_NAFCorr", shapeName=""):
        self.init_list = [p_init]
        self.p = ROOT.RooRealVar("cpow1p_" + cat, "core p1 pow1" + cat, p_init, -10., 5.)
        self.pow = ROOT.RooGenericPdf("core_pow1_" + cat, "core_pow1_" + cat, "@0^@1", ROOT.RooArgList(x,self.p))
        file_= ROOT.TFile("/afs/cern.ch/user/f/fanx/CMSSW_12_6_0_patch1/src/HtoZg_fitting/Core_func_shape/"+fileName+".root", "READ")
        w = file_.Get("w")
        self.ZGMCShape = w.pdf(shapeName)
        w.var(x.GetName()).setRange(x.getMin(), x.getMax())
        self.pdf = ROOT.RooEffProd("corepow1_" +cat + "_model", "Pow1 X Core " + cat, self.ZGMCShape, self.pow)
        self.SBpdf = ROOT.RooGenericPdf("SBcorepow1_" +cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
        self.name = "corepow1_"+ cat
    def checkBond(self):
        tol = 0.001
        par_list = [self.p]
        if any(bondComp(par, tol) for par in par_list):
            print ("The pdf ", self.pdf.GetName(), " needs refit.")
    def reset(self):
        self.p.setVal(self.init_list[0])

class CorePow2Class:
    def __init__(self, x, cat="", p1_init = -15., p2_init = -5., f_init = 1, fileName="ZGCoreShape_01jet_NAFCorr", shapeName=""):
        self.init_list = [p1_init, p2_init, f_init]
        self.p1 = ROOT.RooRealVar("cpow2p1_" + cat, "core p1 pow2" + cat, p1_init, -10., 5.)
        self.p2 = ROOT.RooRealVar("cpow2p2_" + cat, "core p2 pow2" + cat, p2_init, -10., 5.)
        self.f1 = ROOT.RooRealVar("cpow2f1_" + cat, "core f1 pow2" + cat, f_init, -10, 10.)
        self.f2 = ROOT.RooRealVar("cpow2f2_" + cat, "core f2 pow2" + cat, f_init, 10., 10.)

        self.pow = ROOT.RooGenericPdf("core_pow1_" + cat, "core_pow1_" + cat, "@3 /((90.^(1+@1) - 25.^(1+@1))/(1+@1)) * (@0-80)^@1 + @4/((90.^(1+@2) - 25.^(1+@2))/(1+@2)) *(@0-80)^@2", ROOT.RooArgList(x,self.p1, self.p2, self.f1, self.f2))
        # self.pow = ROOT.RooGenericPdf("core_pow2_" + cat, "core_pow2_" + cat, "@3 * @0^@1 + @4*@0^@2", ROOT.RooArgList(x,self.p1, self.p2, self.f1, self.f2))
        file_= ROOT.TFile("/afs/cern.ch/user/f/fanx/CMSSW_12_6_0_patch1/src/HtoZg_fitting/Core_func_shape/"+fileName+".root", "READ")
        w = file_.Get("w")
        self.ZGMCShape = w.pdf(shapeName)
        w.var(x.GetName()).setRange(x.getMin(), x.getMax())
        self.pdf = ROOT.RooEffProd("corepow2_" +cat + "_model", "Pow2 X Core " + cat, self.ZGMCShape, self.pow)
        self.SBpdf = ROOT.RooGenericPdf("SBcorepow2_" +cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
        self.name = "corepow2_"+ cat
    def checkBond(self):
        tol = 0.001
        par_list = [self.p1, self.p2, self.f1, self.f2]
        if any(bondComp(par, tol) for par in par_list):
            print ("The pdf ", self.pdf.GetName(), " needs refit.")
    def reset(self):
        self.p1.setVal(self.init_list[0])
        self.p2.setVal(self.init_list[1])
        self.f1.setVal(self.init_list[2])
        self.f2.setVal(self.init_list[2])

        
