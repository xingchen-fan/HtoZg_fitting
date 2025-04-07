import ROOT

def bondComp(par, tol):
    if par.getValV()+tol > par.getMax() or par.getVal()-tol < par.getMin():
        print (par.GetName(), " = ", par.getValV(), " hits the bondary (", par.getMin(), ", ", par.getMax(), ")")
        return True
    else: return False
    
#### Asymmetric Generalized Gaussian
class AGGClass:
    def __init__(self, x, cat="", kappa_init = -0.5, alpha_init = 20.0, zeta_init = 110.0, x_low = 100., x_high = 180.):
        self.init_list = [kappa_init, alpha_init, zeta_init]
        self.kappa = ROOT.RooRealVar("agg_kappa_" + cat, "agg_kappa_" + cat, kappa_init, -100.0, 100.0)
        self.alpha = ROOT.RooRealVar("agg_alpha_" + cat, "agg_alpha_" + cat, alpha_init, 0.0, 1000.0)
        self.zeta = ROOT.RooRealVar("agg_zeta_" + cat, "agg_zeta_" + cat, zeta_init, x_low, x_high)
        self.pdf = ROOT.AsymGenGaussian("agg_"+cat+"_model", "agg_"+cat, x, self.kappa, self.alpha, self.zeta, x_low, x_high)
        self.SBpdf = ROOT.RooGenericPdf("agg_SB_"+cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
        self.name = "agg_"+cat
    def checkBond(self):
        tol = 0.001
        par_list = [self.kappa, self.alpha, self.zeta]
        if any(bondComp(par, tol) for par in par_list):
            print ("The pdf ", self.pdf.GetName(), " needs refit.")
    def reset(self):
        self.kappa.setVal(self.init_list[0])
        self.alpha.setVal(self.init_list[1])
        self.zeta.setVal(self.init_list[2])
    
#### EXModGaus a.k.a analytic EXP1 function
class EXMGClass:
    def __init__(self, x, cat="", mu_init = 110.0, sig_init = 1.0, xsi_init = 0.5, x_low = 100., x_high = 180.):
        self.init_list = [mu_init, sig_init, xsi_init]
        self.mu = ROOT.RooRealVar("exmg_mu_" + cat, "exmg_mu_" + cat, mu_init, 100, 180)
        self.sig = ROOT.RooRealVar("exmg_sig_" + cat, "exmg_sig_" + cat, sig_init, 0, 10.0)
        self.xsi = ROOT.RooRealVar("exmg_xsi_" + cat, "exmg_xsi_" + cat, xsi_init, -5.0, 1000.0) #When written this is lambda, did not want to overwrite lambda functionality\
        self.pdf = ROOT.EXModGaus("exmg_"+cat+"_model", "exmg_"+cat, x, self.mu, self.sig, self.xsi, x_low, x_high)
        self.SBpdf = ROOT.RooGenericPdf("exmg_SB_"+cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
        self.name = "exmg_"+cat
    def checkBond(self):
        tol = 0.001
        par_list = [self.mu, self.sig, self.xsi]
        if any(bondComp(par, tol) for par in par_list):
            print ("The pdf ", self.pdf.GetName(), " needs refit.")
    def reset(self):
        self.mu.setVal(self.init_list[0])
        self.sig.setVal(self.init_list[1])
        self.xsi.setVal(self.init_list[2])

#### Bernstein
class Bern2Class:
    def __init__(self, x, gauss_mu, cat="", p0=10, p1_init=0.3, p2_init=0.3, bond=20, sigma_init=7., step_init=105., fix_sigma = False):
        self.init_list = [p0, p1_init, p2_init, sigma_init, step_init]
        self.p0 = ROOT.RooRealVar("bern2_p0_" + cat, "bern2_p0_" + cat, p0)
        self.p1 = ROOT.RooRealVar("bern2_p1_" + cat, "bern2_p1_" + cat, p1_init ,-bond, bond)
        self.p2 = ROOT.RooRealVar("bern2_p2_" + cat, "bern2_p2_" + cat, p2_init,-bond, bond)
        self.sigma = ROOT.RooRealVar("bern2_sigma_" + cat,"bern2_sigma_" + cat,sigma_init,  0, 15.)
        self.sigma.setConstant(fix_sigma)
        self.stepval = ROOT.RooRealVar("bern2_step_" + cat, "bern2_step_" + cat, step_init, 90., 120.)
        self.pdf = ROOT.RooGaussStepBernstein("bern2_" +cat + "_model", "Bernstein2 (X) gauss " + cat, x,  gauss_mu, self.sigma, self.stepval, ROOT.RooArgList(self.p0,self.p1,self.p2))
        self.SBpdf = ROOT.RooGenericPdf("bern2_SB_" +cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
        self.name = "bern2_"+ cat
    def checkBond(self):
        tol = 0.001
        par_list = [self.p1, self.p2, self.sigma, self.stepval]
        if any(bondComp(par, tol) for par in par_list):
            print ("The pdf ", self.pdf.GetName(), " needs refit.")
    def reset(self):
         self.p0.setVal(self.init_list[0])
         self.p1.setVal(self.init_list[1])
         self.p2.setVal(self.init_list[2])
         self.sigma.setVal(self.init_list[3])
         self.stepval.setVal(self.init_list[4])

class Bern3Class:
    def __init__(self, x, gauss_mu, cat="", p0=10, p1_init=0.3, p2_init=0.3, p3_init=0.3, bond=20, sigma_init=7., step_init=105., fix_sigma = False):
        self.init_list = [p0, p1_init, p2_init, p3_init, sigma_init, step_init]
        self.p0 = ROOT.RooRealVar("bern3_p0_" + cat, "bern3_p1_" + cat, p0)
        self.p1 = ROOT.RooRealVar("bern3_p1_" + cat, "bern3_p1_" + cat, p1_init,-bond, bond)
        self.p2 = ROOT.RooRealVar("bern3_p2_" + cat, "bern3_p2_" + cat, p2_init,-bond, bond)
        self.p3 = ROOT.RooRealVar("bern3_p3_" + cat, "bern3_p3_" + cat, p3_init,-bond, bond)
        self.sigma = ROOT.RooRealVar("bern3_sigma_" + cat,"bern3_sigma_" + cat, sigma_init,  0, 15.)
        self.sigma.setConstant(fix_sigma)
        self.stepval = ROOT.RooRealVar("bern3_step_" + cat, "bern3_step_" + cat, step_init, 90., 120.)
        self.pdf = ROOT.RooGaussStepBernstein("bern3_" +cat + "_model", "Bernstein3 (X) gauss " + cat, x,  gauss_mu, self.sigma, self.stepval, ROOT.RooArgList(self.p0,self.p1,self.p2,self.p3))
        self.SBpdf = ROOT.RooGenericPdf("bern3_SB_" +cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
        self.name = "bern3_"+ cat
    def checkBond(self):
        tol = 0.001
        par_list = [self.p1, self.p2, self.p3, self.sigma, self.stepval]
        if any(bondComp(par, tol) for par in par_list):
            print ("The pdf ", self.pdf.GetName(), " needs refit.")
    def reset(self):
         self.p0.setVal(self.init_list[0])
         self.p1.setVal(self.init_list[1])
         self.p2.setVal(self.init_list[2])
         self.p3.setVal(self.init_list[3])
         self.sigma.setVal(self.init_list[4])
         self.stepval.setVal(self.init_list[5])

class Bern4Class:
    def __init__(self, x, gauss_mu, cat="", p0=10, p1_init=0.3, p2_init=0.3, p3_init=0.3, p4_init=0.3, bond=20, sigma_init=7., step_init=105., fix_sigma = False):
        self.init_list = [p0, p1_init, p2_init, p3_init, p4_init, sigma_init, step_init]
        self.p0 = ROOT.RooRealVar("bern4_p0_" + cat, "bern4_p1_" + cat, p0)
        self.p1 = ROOT.RooRealVar("bern4_p1_" + cat, "bern4_p1_" + cat, p1_init, -bond, bond)
        self.p2 = ROOT.RooRealVar("bern4_p2_" + cat, "bern4_p2_" + cat, p2_init, -bond, bond)
        self.p3 = ROOT.RooRealVar("bern4_p3_" + cat, "bern4_p3_" + cat, p3_init, -bond, bond)
        self.p4 = ROOT.RooRealVar("bern4_p4_" + cat, "bern4_p4_" + cat, p4_init, -bond, bond)
        self.sigma = ROOT.RooRealVar("bern4_sigma_" + cat,"bern4_sigma_" + cat, sigma_init,  0.1, 15.)
        self.sigma.setConstant(fix_sigma)
        self.stepval = ROOT.RooRealVar("bern4_step_" + cat, "bern4_step_" + cat, step_init, 90., 120.)
        self.pdf = ROOT.RooGaussStepBernstein("bern4_" +cat + "_model", "Bernstein4 (X) gauss " + cat, x,  gauss_mu, self.sigma, self.stepval, ROOT.RooArgList(self.p0,self.p1,self.p2, self.p3, self.p4))
        self.SBpdf = ROOT.RooGenericPdf("bern4_SB_" +cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
        self.name = "bern4_"+ cat
    def checkBond(self):
        tol = 0.001
        par_list = [self.p1, self.p2, self.p3, self.p4, self.sigma, self.stepval]
        if any(bondComp(par, tol) for par in par_list):
            print ("The pdf ", self.pdf.GetName(), " needs refit.")
    def reset(self):
         self.p0.setVal(self.init_list[0])
         self.p1.setVal(self.init_list[1])
         self.p2.setVal(self.init_list[2])
         self.p3.setVal(self.init_list[3])
         self.p4.setVal(self.init_list[4])
         self.sigma.setVal(self.init_list[5])
         self.stepval.setVal(self.init_list[6])

class Bern5Class:
    def __init__(self, x, gauss_mu, cat="", p0=1., p1_init=0.3, p2_init=0.3, p3_init=0.3, p4_init=0.3, p5_init=0.3, bond=20, sigma_init=7., step_init=105., fix_sigma = False):
        self.init_list = [p0, p1_init, p2_init, p3_init, p4_init, p5_init, sigma_init, step_init]
        self.p0 = ROOT.RooRealVar("bern5_p0_" + cat, "bern5_p1_" + cat, p0)
        self.p1 = ROOT.RooRealVar("bern5_p1_" + cat, "bern5_p1_" + cat, p1_init,-bond, bond)
        self.p2 = ROOT.RooRealVar("bern5_p2_" + cat, "bern5_p2_" + cat, p2_init,-bond, bond)
        self.p3 = ROOT.RooRealVar("bern5_p3_" + cat, "bern5_p3_" + cat, p3_init,-bond, bond)
        self.p4 = ROOT.RooRealVar("bern5_p4_" + cat, "bern5_p4_" + cat, p4_init,-bond, bond)
        self.p5 = ROOT.RooRealVar("bern5_p5_" + cat, "bern5_p5_" + cat, p5_init,-bond, bond)
        self.sigma = ROOT.RooRealVar("bern5_sigma_" + cat,"sigma_bern5_" + cat,sigma_init,  0.1, 15.)
        self.sigma.setConstant(fix_sigma)
        self.stepval = ROOT.RooRealVar("bern5_step_" + cat, "bern5_step_" + cat, step_init, 90., 120.)
        self.pdf = ROOT.RooGaussStepBernstein("bern5_" +cat + "_model", "Bernstein5 (X) gauss " + cat, x,  gauss_mu, self.sigma, self.stepval, ROOT.RooArgList(self.p0,self.p1,self.p2, self.p3, self.p4, self.p5))
        self.SBpdf = ROOT.RooGenericPdf("bern5_SB_" +cat + "_model", "bern5_SB_" +cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
        self.name = "bern5_"+ cat
    def checkBond(self):
        tol = 0.001
        par_list = [self.p1, self.p2, self.p3, self.p4, self.p5, self.sigma, self.stepval]
        if any(bondComp(par, tol) for par in par_list):
            print ("The pdf ", self.pdf.GetName(), " needs refit.")
    def reset(self):
         self.p0.setVal(self.init_list[0])
         self.p1.setVal(self.init_list[1])
         self.p2.setVal(self.init_list[2])
         self.p3.setVal(self.init_list[3])
         self.p4.setVal(self.init_list[4])
         self.p5.setVal(self.init_list[5])
         self.sigma.setVal(self.init_list[6])
         self.stepval.setVal(self.init_list[7])

#### Pow series
class Pow1Class:
    def __init__(self, x, gauss_mu, cat="", sigma_init = 7., step_init = 108., p_init = -6., p_low = -10., p_high = 0., di_gauss = False, fix_sigma = False):
        self.init_list = [sigma_init, step_init, p_init]
        self.t = ROOT.RooRealVar("pow1_t_" + cat, "t pow1" + cat, step_init, 90., 110.)
        self.p = ROOT.RooRealVar("pow1_p_" + cat, "p1 pow1" + cat, p_init, p_low, p_high)
        self.sigma = ROOT.RooRealVar("pow1_sigma_" + cat,"pow1_sigma_"+cat, sigma_init,  0., 15.)
        self.sigma.setConstant(fix_sigma)
        self.sigma2 = ROOT.RooRealVar("sigma2_pow1_" + cat,"sigma2_pow1_"+cat, sigma_init * 2,  0.1, 15.)
        self.gauss = ROOT.RooGaussian("gaussxpow1_"+cat, "gaussian PDF pow1 " + cat, x, gauss_mu, self.sigma)
        self.gauss2 = ROOT.RooGaussian("gauss2xpow1_"+cat, "gaussian2 PDF pow1 " + cat, x, gauss_mu, self.sigma2)
        self.gc = ROOT.RooRealVar("gc_pow1_" + cat, "gc_pow1_" + cat, 1, 0, 10)
        self.addG = ROOT.RooAddPdf("addG_pow1_" + cat, "addG_pow1_" + cat, self.gauss2, self.gauss, self.gc)
        self.transX = ROOT.RooFormulaVar("transX_"+cat, "@0 - 100", ROOT.RooArgList(x)) 
        self.transt = ROOT.RooFormulaVar("transt_"+cat, "@0 - 100", ROOT.RooArgList(self.t))
        self.step = ROOT.RooGenericPdf("pow1_step" + cat, "pow1_step_" + cat,\
                                        "( ((@0-@1)*153.85 <0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85) ) )*(@0^@2)", ROOT.RooArgList(x,self.t,self.p))
        x.setBins(20000, "cache")
        if di_gauss == True:
            self.pdf = ROOT.RooFFTConvPdf("pow1_" + cat + "_model", "step pow1 (X) gauss" + cat, x, self.step, self.addG)
        else:
            self.pdf = ROOT.RooFFTConvPdf("pow1_" + cat + "_model", "step pow1 (X) gauss" + cat, x, self.step, self.gauss)
        self.pdf.setBufferFraction(0.5)
        self.SBpdf = ROOT.RooGenericPdf("pow1_SB_" +cat + "_model", "pow1_SB_" +cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
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
    def __init__(self, x, gauss_mu, cat="", sigma_init = 7., step_init = 107., p1_init = -5.8, p1_low = -8., p1_high = 0., p2_init = -15, p2_low = -20., p2_high = -5., f1_init = 0.9, f2_init = 0.9, xmax = 165., const_f1 = False, di_gauss = False, fix_sigma = False):
        self.init_list = [sigma_init, step_init, p1_init, p2_init, f1_init, f2_init]
        self.xmax = ROOT.RooRealVar("pow2_xmax_" + cat, "pow2_xmax_" + cat, xmax)
        self.t = ROOT.RooRealVar("pow2_t_" + cat, "t pow2" + cat, step_init, 90., 118.)
        self.p1 = ROOT.RooRealVar("pow2_p1_" + cat, "p1 pow2" + cat, p1_init, p1_low, p1_high)
        self.p2 = ROOT.RooRealVar("pow2_p2_" + cat, "p2 pow2" + cat, p2_init, p2_low, p2_high)
        self.f1 = ROOT.RooRealVar("pow2_f1_" + cat, "f1 pow2" + cat, f1_init, 0., 10.)
        self.f2 = ROOT.RooRealVar("pow2_f2_" + cat, "f2 pow2" + cat, f2_init, 0., 10.)
        self.sigma = ROOT.RooRealVar("pow2_sigma_" + cat,"pow2_sigma_"+cat, sigma_init,  0.1, 15.)
        self.sigma.setConstant(fix_sigma)
        self.sigma2 = ROOT.RooRealVar("pow2_sigma2_" + cat,"pow2_sigma2_"+cat, sigma_init * 2,  0.1, 15.)
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
        self.norm1 = ROOT.RooFormulaVar("pow2_norm1_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.t, self.p1, self.xmax))
        self.norm2 = ROOT.RooFormulaVar("pow2_norm2_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.t, self.p2, self.xmax))
        self.step = ROOT.RooGenericPdf("pow2_step_" + cat, "pow2_step_" + cat,"( ((@0-@1)*153.85 <0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85) ) )*(@5/@6 * @0^@2 + @4 / @7 * @0^@3)", ROOT.RooArgList(x, self.t,self.p1, self.p2, self.f1, self.f2, self.norm1, self.norm2))
        x.setBins(20000, "cache")
        if di_gauss == True:
            self.pdf = ROOT.RooFFTConvPdf("pow2_" + cat + "_model", "step pow2 (X) gauss" + cat, x, self.step, self.addG)
        else:
            self.pdf = ROOT.RooFFTConvPdf("pow2_" + cat + "_model", "step pow2 (X) gauss" + cat, x, self.step, self.gauss)
        self.pdf.setBufferFraction(0.5)
        self.SBpdf = ROOT.RooGenericPdf("pow2_SB_" +cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
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
    def __init__(self, x, gauss_mu, cat="", sigma_init = 7., step_init = 108., p1_init = -5.7, p1_low = -8., p1_high = 0., p2_init = -8., p2_low = -20., p2_high = -5., p3_init = -10., p3_low = -40., p3_high = -5., f1_init = 0.9, f2_init = 0.2, f3_init = 0.3, xmax = 165., const_f1 = False, di_gauss = False, fix_sigma = False):
        self.init_list = [sigma_init, step_init, p1_init, p2_init, p3_init, f1_init, f2_init,  f3_init]
        self.xmax = ROOT.RooRealVar("pow3_xmax_" + cat, "pow3_xmax_" + cat, xmax)
        self.t = ROOT.RooRealVar("pow3_t_" + cat, "t pow3" + cat, step_init, 90., 118.)
        self.p1 = ROOT.RooRealVar("pow3_p1_" + cat, "p1 pow3" + cat, p1_init, p1_low, p1_high)
        self.p2 = ROOT.RooRealVar("pow3_p2_" + cat, "p2 pow3" + cat, p2_init, p2_low, p2_high)
        self.p3 = ROOT.RooRealVar("pow3_p3_" + cat, "p3 pow3" + cat, p3_init, p3_low, p3_high)
        self.f1 = ROOT.RooRealVar("pow3_f1_" + cat, "f1 pow3" + cat, f1_init, 0., 10.)
        self.f2 = ROOT.RooRealVar("pow3_f2_" + cat, "f2 pow3" + cat, f2_init, 0., 10.)
        self.f3 = ROOT.RooRealVar("pow3_f3_" + cat, "f3 pow3" + cat, f3_init, 0., 10.)
        self.offset = ROOT.RooRealVar("os_" + cat, "os_" + cat, 1e-10)
        if const_f1 == True:
            self.f1.setConstant(True)
        self.sigma = ROOT.RooRealVar("pow3_sigma_" + cat,"pow3_sigma_"+cat, sigma_init,  1., 15.)
        self.sigma.setConstant(fix_sigma)
        self.sigma2 = ROOT.RooRealVar("pow3_sigma2_" + cat,"sigma2_pow3_"+cat, sigma_init * 2,  0.1, 15.)
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
        self.norm1 = ROOT.RooFormulaVar("pow3_norm1_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.t, self.p1, self.xmax))
        self.norm2 = ROOT.RooFormulaVar("pow3_norm2_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.t, self.p2, self.xmax))
        self.norm3 = ROOT.RooFormulaVar("pow3_norm3_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.t, self.p3, self.xmax))
        self.step = ROOT.RooGenericPdf("pow3_step_" + cat, "step_pow3" + cat, "( ((@0-@1)*307.69 <0.0) ? 0.0 : (((@0-@1)*307.69 >1.0) ? 1.0 : ((@0-@1)*307.69) ) )*( @5/ @8 *@0^@2 + @6/ @9 * @0^@3 + @7/ @10 * @0^@4)", ROOT.RooArgList(x,self.t,self.p1, self.p2, self.p3, self.f1, self.f2, self.f3, self.norm1, self.norm2, self.norm3))
        x.setBins(40000, "cache")
        if di_gauss == True:
            self.pdf = ROOT.RooFFTConvPdf("pow3_" + cat + "_model", "step pow3 (X) gauss " + cat, x, self.step, self.addG)
        else:
            self.pdf = ROOT.RooFFTConvPdf("pow3_" + cat + "_model", "step pow3 (X) gauss " + cat, x, self.step, self.gauss)
        self.pdf.setBufferFraction(0.5)
        self.SBpdf = ROOT.RooGenericPdf("pow3_SB_" +cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
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
    def __init__(self, x, gauss_mu, cat="", sigma_init = 7., step_init = 108., p_init = -0.02, p_low = -2, p_high = 0., di_gauss = False, fix_sigma = False):
        self.init_list = [sigma_init, step_init, p_init]
        self.sigma = ROOT.RooRealVar("exp1_sigma_" + cat,"sigma_exp1" + cat       ,sigma_init,  0., 15.)
        self.sigma.setConstant(fix_sigma)
        self.sigma2 = ROOT.RooRealVar("exp1_sigma2_" + cat,"sigma2_exp1_"+cat, sigma_init * 2,  0.1, 15.)
        self.t = ROOT.RooRealVar("exp1_t_" + cat, "t exp1 " + cat, step_init, 90., 118.)
        self.p1 = ROOT.RooRealVar("exp1_p1_" + cat, "p1 exp1 " +cat, p_init, p_low, p_high)
        self.step = ROOT.RooGenericPdf("exp1_step_" + cat, "step_exp1 " + cat, "( ((@0-@1)*153.85<0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85) ) )*(TMath::Exp(@0*@2))", ROOT.RooArgList(x,self.t,self.p1))
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
        self.SBpdf = ROOT.RooGenericPdf("exp1_SB_" +cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
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
    def __init__(self, x, gauss_mu, cat="", sigma_init = 7., step_init = 105., p1_init = -0.02, p1_low = -1, p1_high = 0., p2_init = -0.02, p2_low = -1, p2_high = 0., f1_init = 0.5, f2_init = 0.5, xmax = 165., const_f1 = False, di_gauss = False, fix_sigma = False):
        self.init_list = [sigma_init, step_init, p1_init, p2_init, f1_init, f2_init]
        self.xmax = ROOT.RooRealVar("exp2_xmax_" + cat, "xmax_exp2_" + cat, xmax)
        self.sigma = ROOT.RooRealVar("exp2_sigma_" + cat,"sigma_exp2" + cat       ,sigma_init,  1., 15.)
        self.sigma.setConstant(fix_sigma)
        self.sigma2 = ROOT.RooRealVar("exp2_sigma2_" + cat,"sigma2_exp2_"+cat, sigma_init * 2,  0.1, 15.)
        self.offset = ROOT.RooRealVar("os_" + cat, "os_" + cat, 1e-10)
        self.t = ROOT.RooRealVar("exp2_t_" + cat, "t exp2 " + cat, step_init, 90., 115.)
        self.p1 = ROOT.RooRealVar("exp2_p1_" + cat, "p1 exp2 " +cat, p1_init, p1_low, p1_high)
        self.p2 = ROOT.RooRealVar("exp2_p2_" + cat, "p2 exp2 " +cat, p2_init, p2_low, p2_high)
        self.f1 = ROOT.RooRealVar("exp2_f1_" + cat, "f1 exp2 "+ cat, f1_init, 0., 10.)
        self.f2 = ROOT.RooRealVar("exp2_f2_" + cat, "f2 exp2 "+ cat, f2_init, 0., 10.)
        if const_f1 == True:
            self.f1.setConstant(True)
        self.norm1 = ROOT.RooFormulaVar("exp2_norm1_"+cat, "(TMath::Exp(@2*@1) - TMath::Exp(@0*@1))/@1", ROOT.RooArgList(self.t, self.p1, self.xmax))
        self.norm2 = ROOT.RooFormulaVar("exp2_norm2_"+cat, "(TMath::Exp(@2*@1) - TMath::Exp(@0*@1))/@1", ROOT.RooArgList(self.t, self.p2, self.xmax))
        self.step = ROOT.RooGenericPdf("exp2_step_" + cat, "step_exp2_" + cat, "( ((@0-@1)*153.85<0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85) ) )*(@5 / @6 * TMath::Exp(@0*@2) + @4 / @7 * TMath::Exp(@0*@3))", ROOT.RooArgList(x,self.t,self.p1,self.p2,self.f1, self.f2, self.norm1, self.norm2))
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
        self.SBpdf = ROOT.RooGenericPdf("exp2_SB_" +cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
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
    def __init__(self, x, gauss_mu, cat="", sigma_init = 7., step_init = 108., p1_init = -0.02, p1_low = -2, p1_high = 0., p2_init = -0.02, p2_low = -2, p2_high = 0., p3_init = -0.02, p3_low = -2, p3_high = 0., f1_init = 0.33, f2_init = 0.5, f3_init = 0.5, xmax = 165., const_f1 = False, di_gauss = False, fix_sigma = False):
        self.init_list = [sigma_init, step_init, p1_init, p2_init, p3_init, f1_init, f2_init, f3_init]
        self.sigma = ROOT.RooRealVar("exp3_sigma_" + cat,"exp3_sigma_" + cat ,sigma_init,  1., 15.)
        self.sigma.setConstant(fix_sigma)
        self.sigma2 = ROOT.RooRealVar("exp3_sigma2_" + cat,"exp3_sigma2_"+cat, sigma_init * 2,  0.1, 15.)
        self.xmax = ROOT.RooRealVar("exp3_xmax_" + cat, "exp3_xmax_" + cat, xmax)
        self.t = ROOT.RooRealVar("exp3_t_" + cat, "t exp3 " + cat, step_init, 90., 115.)
        self.p1 = ROOT.RooRealVar("exp3_p1_" + cat, "p1 exp3 " +cat, p1_init, p1_low, p1_high)
        self.p2 = ROOT.RooRealVar("exp3_p2_" + cat, "p2 exp3 " +cat, p2_init, p2_low, p2_high)
        self.p3 = ROOT.RooRealVar("exp3_p3_" + cat, "p3 exp3 " +cat, p3_init, p3_low, p3_high)
        self.f1 = ROOT.RooRealVar("exp3_f1_" + cat, "f1 exp3 "+ cat, f1_init, 0, 10.)
        self.f2 = ROOT.RooRealVar("exp3_f2_" + cat, "f2 exp3 "+ cat, f2_init, 0, 10.)
        self.f3 = ROOT.RooRealVar("exp3_f3_" + cat, "f3 exp3 "+ cat, f3_init, 0., 10.)
        if const_f1 == True:
            self.f1.setConstant(True)
        self.norm1 = ROOT.RooFormulaVar("exp3_norm1_"+cat, "(TMath::Exp(@2*@1) - TMath::Exp(@0*@1))/@1", ROOT.RooArgList(self.t, self.p1, self.xmax))
        self.norm2 = ROOT.RooFormulaVar("exp3_norm2_"+cat, "(TMath::Exp(@2*@1) - TMath::Exp(@0*@1))/@1", ROOT.RooArgList(self.t, self.p2, self.xmax))
        self.norm3 = ROOT.RooFormulaVar("exp3_norm3_"+cat, "(TMath::Exp(@2*@1) - TMath::Exp(@0*@1))/@1", ROOT.RooArgList(self.t, self.p3, self.xmax))
        self.step = ROOT.RooGenericPdf("exp3_step_" + cat, "exp3_step_" + cat, "( ((@0-@1)*153.85<0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85) ) )*( @5 / @8 * TMath::Exp(@0*@2) + @6 / @9 * TMath::Exp(@0*@3) + @7 / @10  * TMath::Exp(@0*@4))", ROOT.RooArgList(x,self.t,self.p1,self.p2, self.p3, self.f1, self.f2, self.f3, self.norm1, self.norm2, self.norm3))
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
        self.SBpdf = ROOT.RooGenericPdf("exp3_SB_" +cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
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
    def __init__(self, x, gauss_mu, cat="", sigma_init = 7., step_init = 108., p1 = -7, p2 = -6, f_init = 0.1, xmax = 165., const_f1 = False, di_gauss = False, fix_sigma = False):
        self.init_list = [sigma_init, step_init, p1, p2, f_init]
        self.sigma = ROOT.RooRealVar("lau2_sigma_" + cat,"lau2_sigma_" + cat       ,sigma_init,  1., 40.)
        self.sigma.setConstant(fix_sigma)
        self.sigma2 = ROOT.RooRealVar("lau2_sigma2_" + cat,"lau2_sigma2_"+cat, sigma_init * 2,  0.1, 40.)
        self.offset = ROOT.RooRealVar("lau2_os_" + cat, "lau2_os_" + cat, 1e-10)
        self.xmax = ROOT.RooRealVar("lau2_xmax_" + cat, "lau2_xmax_" + cat, xmax)
        self.t = ROOT.RooRealVar("lau2_t_" + cat, "t lau2 " + cat, step_init, 90., 118.)
        self.f1 = ROOT.RooRealVar("lau2_f1_" + cat, "f1 lau2 " + cat, f_init, 0, 30.)
        self.f2 = ROOT.RooRealVar("lau2_f2_" + cat, "f2 lau2 " + cat, f_init, 0., 30.)
        self.p1 = ROOT.RooRealVar("lau2_p1_" + cat, "p1 lau2 " + cat, p1)
        self.p2 = ROOT.RooRealVar("lau2_p2_" + cat, "p2 lau2 " + cat, p2)
        self.norm1 = ROOT.RooFormulaVar("lau2_norm1_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.t, self.p1, self.xmax))
        self.norm2 = ROOT.RooFormulaVar("lau2_norm2_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.t, self.p2, self.xmax))
        self.step = ROOT.RooGenericPdf("lau2_step_" + cat, "lau2_step_" + cat, "( ((@0-@1)*153.85<0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85) ) )*(@2/@6*@0^(@4) + @3/@7*@0^(@5))", \
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
        self.SBpdf = ROOT.RooGenericPdf("lau2_SB_" +cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
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
    def __init__(self, x, gauss_mu, cat="", sigma_init = 7., step_init = 108., p1 = -8, p2 = -7, p3 = -6, f_init = 0.1, xmax = 165., const_f1 = False, di_gauss = False, fix_sigma = False):
        self.init_list = [sigma_init, step_init, p1, p2, p3, f_init]
        self.sigma = ROOT.RooRealVar("lau3_sigma_" + cat,"lau3_sigma_" + cat       ,sigma_init,  1., 25.)
        self.sigma.setConstant(fix_sigma)
        self.sigma2 = ROOT.RooRealVar("lau3_sigma2_" + cat,"lau3_sigma2_"+cat, sigma_init * 2,  0.1, 25.)
        self.t = ROOT.RooRealVar("lau3_t_" + cat, "t lau3 " + cat, step_init, 90., 118.)
        self.xmin = ROOT.RooRealVar("lau3_xmin_" + cat, "lau3_xmin_" + cat,  100)
        self.xmax = ROOT.RooRealVar("lau3_xmax_" + cat, "lau3_xmax_" + cat, xmax)
        self.offset = ROOT.RooRealVar("os_" + cat, "os_" + cat, 1e-5, 0., 0.01)
        self.f1 = ROOT.RooRealVar("lau3_f1_" + cat, "f1 lau3 " + cat, f_init, 0., 10.)
        self.f2 = ROOT.RooRealVar("lau3_f2_" + cat, "f2 lau3 " + cat, f_init, 0., 10.)
        self.f3 = ROOT.RooRealVar("lau3_f3_" + cat, "f3 lau3 " + cat, f_init, 0., 10.)
        self.p1 = ROOT.RooRealVar("lau3_p1_" + cat, "p1 lau3 " + cat, p1)
        self.p2 = ROOT.RooRealVar("lau3_p2_" + cat, "p2 lau3 " + cat, p2)
        self.p3 = ROOT.RooRealVar("lau3_p3_" + cat, "p3 lau3 " + cat, p3)
        self.norm1 = ROOT.RooFormulaVar("lau3_norm1_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.xmin, self.p1, self.xmax))
        self.norm2 = ROOT.RooFormulaVar("lau3_norm2_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.xmin, self.p2, self.xmax))
        self.norm3 = ROOT.RooFormulaVar("lau3_norm3_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.xmin, self.p3, self.xmax))
        self.step = ROOT.RooGenericPdf("lau3_step_" + cat, "lau3_step_" + cat, "( ((@0-@1)*153.85<0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85) ) )*(@2/@8*@0^(@4) + @3/@9*@0^(@5) + @7/@10 * @0^(@6))", \
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
        self.SBpdf = ROOT.RooGenericPdf("lau3_SB_" +cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
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
    def __init__(self, x, gauss_mu, cat="", sigma_init = 7., step_init = 108., p1 = -8, p2 = -7, p3 = -6, p4 = -5, f_init = 0.1, xmax = 165., const_f1 = False, di_gauss = False, fix_sigma = False):
        self.init_list = [sigma_init, step_init, p1, p2, p3, p4, f_init]
        self.sigma = ROOT.RooRealVar("lau4_sigma_" + cat,"lau4_sigma_" + cat , sigma_init,  1., 25.)
        self.sigma.setConstant(fix_sigma)
        self.sigma2 = ROOT.RooRealVar("lau4_sigma2_" + cat,"lau4_sigma2_"+cat, sigma_init * 2,  0.1, 25.)
        self.t = ROOT.RooRealVar("lau4_t_" + cat, "t lau4 " + cat, step_init, 90., 118.)
        self.xmax = ROOT.RooRealVar("lau4_xmax_" + cat, "lau4_xmax_" + cat, xmax)
        self.offset = ROOT.RooRealVar("os_" + cat, "os_" + cat, 1e-5, 0., 0.01)
        self.f1 = ROOT.RooRealVar("lau4_f1_" + cat, "f1 lau4 " + cat, f_init, 0., 10.)
        self.f2 = ROOT.RooRealVar("lau4_f2_" + cat, "f2 lau4 " + cat, f_init, 0., 10.)
        self.f3 = ROOT.RooRealVar("lau4_f3_" + cat, "f3 lau4 " + cat, f_init, 0., 10.)
        self.f4 = ROOT.RooRealVar("lau4_f4_" + cat, "f4 lau4 " + cat, f_init, 0., 10.)
        self.p1 = ROOT.RooRealVar("lau4_p1_" + cat, "p1 lau4 " + cat, p1)
        self.p2 = ROOT.RooRealVar("lau4_p2_" + cat, "p2 lau4 " + cat, p2)
        self.p3 = ROOT.RooRealVar("lau4_p3_" + cat, "p3 lau4 " + cat, p3)
        self.p4 = ROOT.RooRealVar("lau4_p4_" + cat, "p4 lau4 " + cat, p4)
        self.norm1 = ROOT.RooFormulaVar("lau4_norm1_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.t, self.p1, self.xmax))
        self.norm2 = ROOT.RooFormulaVar("lau4_norm2_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.t, self.p2, self.xmax))
        self.norm3 = ROOT.RooFormulaVar("lau4_norm3_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.t, self.p3, self.xmax))
        self.norm4 = ROOT.RooFormulaVar("lau4_norm4_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.t, self.p4, self.xmax))

        self.step = ROOT.RooGenericPdf("lau4_step_" + cat, "lau4_step_" + cat, "( ((@0-@1)*153.85<0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85) ) )*(@2/@10*@0^(@5) + @3/@11*@0^(@6) + @4/@12*@0^(@7) + @9/@13*@0^(@8))", \
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
        self.SBpdf = ROOT.RooGenericPdf("lau4_SB_" +cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
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
        self.sigma = ROOT.RooRealVar("lau5_sigma_" + cat,"lau5_sigma_" + cat       ,sigma_init,  1., 25.)
        self.xmax = ROOT.RooRealVar("lau5_xmax_" + cat, "lau5_xmax_" + cat, xmax)
        self.t = ROOT.RooRealVar("lau5_t_" + cat, "t lau5 " + cat, step_init, 90., 118.)
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
        self.norm1 = ROOT.RooFormulaVar("lau5_norm1_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.t, self.p1, self.xmax))
        self.norm2 = ROOT.RooFormulaVar("lau5_norm2_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.t, self.p2, self.xmax))
        self.norm3 = ROOT.RooFormulaVar("lau5_norm3_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.t, self.p3, self.xmax))
        self.norm4 = ROOT.RooFormulaVar("lau5_norm4_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.t, self.p4, self.xmax))
        self.norm5 = ROOT.RooFormulaVar("lau5_norm5_"+cat, "(@2^(1+@1) - @0^(1+@1))/(1+@1)", ROOT.RooArgList(self.t, self.p5, self.xmax))
        
        self.step = ROOT.RooGenericPdf("lau5_step_" + cat, "step_" + cat, "( ((@0-@1)*153.85<0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85) ) )*(@2/@12*@0^@7 + @3/@13*@0^@8 + @4/@14*@0^@9 + @5/@15*@0^@10 + @6/@16*@0^@11)", \
                                      ROOT.RooArgList(x,self.t, self.f1, self.f2, self.f3, self.f4, self.f5, self.p1, self.p2, self.p3, self.p4, self.p5, self.norm1, self.norm2, self.norm3, self.norm4, self.norm5))
        self.gauss = ROOT.RooGaussian("gaussxlau5_" + cat, "gaussian PDF lau5 " + cat, x, gauss_mu, self.sigma)
        if const_f1 == True:
            self.f1.setConstant(True)
        x.setBins(20000, "cache")
        self.pdf = ROOT.RooFFTConvPdf("lau5_" +cat+ "_model", "step lau5 (X) gauss " + cat, x, self.step, self.gauss)
        self.pdf.setBufferFraction(0.5)
        self.SBpdf = ROOT.RooGenericPdf("lau5_SB_" +cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
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
        self.m0 = ROOT.RooRealVar("modg_m0_"+cat,"mass peak value [GeV]" , m0,100,130)
        self.vl = ROOT.RooRealVar("modg_nuL_" +cat,"low-end power"       , vl,  0,  15)
        self.vr = ROOT.RooRealVar("modg_nuRange_"+cat,"power range"       , vr, -5,  5)
        self.s0 = ROOT.RooRealVar("modg_NOUSE_sigma0_"+cat,"peak width"       , 0)#,  1., 10.)
        self.sl = ROOT.RooRealVar("modg_sigmaL_"+cat,"low-end width"    , sl, 1., 40)
        self.sh = ROOT.RooRealVar("modg_sigmaH_"+cat,"high-end width"   , sh, 1.,60)
        self.pdf = ROOT.ModGaus("modg_"+cat+"_model","modg_"+cat+"_model", x, self.m0, self.vl, self.vr, self.s0, self.sl, self.sh, lowx, highx)
        self.SBpdf = ROOT.RooGenericPdf("modg_SB_" +cat + "_model", "((@0 < 120)? 1:((@0 > 130)? 1:0)) * @1", ROOT.RooArgList(x, self.pdf))
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

        
