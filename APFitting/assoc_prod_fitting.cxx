//#include <RooFit.h>
//#include <RooFFTConvPdf.h>
////#include "../Utilities/RooGaussStepBernstein.cxx"
//#include "RooFit/RooCrystalBall.h"
//#include "../Utilities/ModGaus.cxx"

using namespace RooFit;
using namespace std;
#include "../Utilities/ModGaus.h"
#include "../Utilities/ModGaus.cxx"

#include "../Utilities/EXModGaus.h"
#include "../Utilities/EXModGaus.cxx"

#include "../Utilities/AsymGenGaussian.h"
#include "../Utilities/AsymGenGaussian.cxx"

#include "../Utilities/RooGaussStepBernstein.h"

double br_left  = 120;
double br_right = 130;

std::vector<string> filename_vector(string prefix, std::vector<string> &filenames){
  std::vector<string> ret_vector = {};
  for(unsigned int idx = 0; idx < filenames.size(); idx++){
    ret_vector.push_back(prefix + filenames[idx]);
  }
  return ret_vector;
}


TH1D* get_hist(string infile, string hist_name){
  TFile *file = new TFile(infile.c_str());
  TCanvas* canvas_hist = static_cast<TCanvas*>(file->Get("canvas"));
  TH1D *hist_return = nullptr;
  //Finds histogram in file
  for(int idx_h=0; idx_h<50000; idx_h++){
    if(hist_return == nullptr){
      hist_return = static_cast<TH1D*>( canvas_hist->FindObject( (hist_name+std::to_string(idx_h)).c_str() ) );
    }

    if(hist_return != nullptr){break;}
  }

  return hist_return;
}

TH1D* get_bkg_hist(string infile){ return get_hist(infile,"bkg_Background_");}
TH1D* get_sig_hist(string infile){ return get_hist(infile,"sig_Signal_");}


  /*
  //Bernstein functions
  //Bern 5
  RooRealVar p0("p0", "p0", 15.);
  RooRealVar b5p1("b5p1", "b5p1", 0.3, -200., 200.);
  RooRealVar b5p2("b5p2", "b5p2", 0.3,-200.,200.);
  RooRealVar b5p3("b5p3", "b5p3", 0.3,-200.,200.);
  RooRealVar b5p4("b5p4", "b5p4", 0.3,-200.,200.);
  RooRealVar b5p5("b5p5", "b5p5", 0.3,-200.,200.);

  RooRealVar sigma_bern5("sigma_bern5","sigma_bern5"       ,4.7647,  1., 20.);
  RooRealVar stepval_bern5("stepval_bern5", "stepval_bern5", 103.46, 100., 120.);
  RooGaussStepBernstein bern5_model("bern5_model", "Bernstein 5(X) gauss", m_lly, meanx, sigma_bern5, stepval_bern5, RooArgList(p0,b5p1,b5p2,b5p3, b5p4, b5p5));
  RooFitResult* ber5_fit = ber5_model.fitTo(hist_total_rdh,SumW2Error(kTRUE) ,Save(kTRUE));

  //Bern 2
  RooGaussStepBernstein bern2_model("bern2_model", "Bernstein 2(X) gauss", m_lly, meanx, sigma_bern5, stepval_bern5, RooArgList(p0,b5p1,b5p2));
  RooFitResult* ber2_fit = ber2_model.fitTo(hist_total_rdh,SumW2Error(kTRUE) ,Save(kTRUE));

  //Bern 3
  RooGaussStepBernstein bern3_model("bern3_model", "Bernstein 3(X) gauss", m_lly, meanx, sigma_bern5, stepval_bern5, RooArgList(p0,b5p1,b5p2, b5p3));
  RooFitResult* ber3_fit = ber3_model.fitTo(hist_total_rdh,SumW2Error(kTRUE) ,Save(kTRUE));

  // Bern 4 
  RooGaussStepBernstein bern4_model("bern4_model", "Bernstein 4(X) gauss", m_lly, meanx, sigma_bern5, stepval_bern5, RooArgList(p0,b5p1,b5p2,b5p3, b5p4));
  RooFitResult* ber4_fit = ber4_model.fitTo(hist_total_rdh,SumW2Error(kTRUE) ,Save(kTRUE));
  //End of Bernstein functions
  */

RooCrystalBall* fit_sig_model(RooRealVar &x, RooDataHist &hist_to_fit, const char *name){
  RooRealVar *a1_dcb       = new RooRealVar(     "a1",   "a1",  0.5,  0.2,   1);
  RooRealVar *a2_dcb       = new RooRealVar(     "a2",   "a2",  0.5,  0.2,   1);
  RooRealVar *n1_dcb       = new RooRealVar(     "n1",   "n1",  25.,  4.0, 50.);
  RooRealVar *n2_dcb       = new RooRealVar(     "n2",   "n2",  25.,  4.0, 50.);
  RooRealVar *sigma_dcb    = new RooRealVar(  "sigma","sigma",    2,  0.1, 3.0);
  RooRealVar *dm_dcb       = new RooRealVar(     "dm",   "dm", -0.1, -1.0, 1.0);
  RooFormulaVar *mean_dcb  = new RooFormulaVar("mean", "mean", "(125+@0)", RooArgList(*dm_dcb));

  RooCrystalBall *pdf  = new RooCrystalBall(name,name,x,*mean_dcb,*sigma_dcb,*a1_dcb,*n1_dcb,*a2_dcb,*n2_dcb);
  pdf -> fitTo(hist_to_fit, SumW2Error(kTRUE) ,Save(kTRUE));
  return pdf;
}

RooGaussian* fit_sig_gauss(RooRealVar &x, RooDataHist &hist_to_fit, const char *name){
  RooRealVar    *dmean = new RooRealVar(   "dm",   "dm",   -0.1, -1.0, 1.0);
  RooRealVar    *sig   = new RooRealVar(   "sig",  "sigma", 2.0,  1.0, 3.0); sig -> setError(0.1);
  RooFormulaVar *mean  = new RooFormulaVar("mean", "mean", "(125+@0)", RooArgList(*dmean));
  RooGaussian   *pdf = new RooGaussian(  name,   "gaussian", x, *mean, *sig);
  
  pdf -> fitTo(hist_to_fit, SumW2Error(kTRUE) ,Save(kTRUE));
  return pdf;
}


ModGaus* fit_modgauss_model(RooRealVar &x, RooDataHist &hist_to_fit, const char *name){
  //Modified Gaussian fitting function
  RooRealVar* m0 = new RooRealVar("m_{0}","mass peak value [GeV]",110, 100, 180); m0 -> setError(1);
  RooRealVar* vl = new RooRealVar("#nu_{L}","low-end power"      ,1.0,  -2,  10); vl -> setError(0.1);
  RooRealVar* vr = new RooRealVar("#Delta#nu","power range"      ,0.1,  -5,   5); vr -> setError(0.1);
  RooRealVar* s0 = new RooRealVar("#sigma_{0}","peak width"      , 15,  5.,  30); s0 -> setError(1);
  RooRealVar* sl = new RooRealVar("#sigma_{L}","low-end width"   ,  5,  1.,  40); sl -> setError(1);
  RooRealVar* sh = new RooRealVar("#sigma_{H}","high-end width"  , 20,  5.,  60); sh -> setError(1);
  ModGaus* modgauss_model = new ModGaus(name,"G_{M}", x, *m0, *vl, *vr, *s0, *sl, *sh, 100, 180);
  modgauss_model -> fitTo(hist_to_fit,SumW2Error(kTRUE) ,Save(kTRUE));
  //End of Modified Gaussian fitting function
  return modgauss_model;
}

//This is the analytical form of exp1 function with fewer degrees of freedom
EXModGaus* fit_exmg_model(RooRealVar &x, RooDataHist &hist_to_fit, const char *name, double low_d=100, double high_d=180){
  //Modified Gaussian fitting function
  //new range is -1 to 1
  double range_d = 2;
  //RooRealVar *range    = new RooRealVar("range","range", range_d);
  //RooRealVar *low      = new RooRealVar(  "low",  "low", low_d);
  //RooRealVar *high     = new RooRealVar("high", "high", high_d);
  //RooFormulaVar* mod_x = new RooFormulaVar("mod_x","mod_x", "2.0/(@2-@1)*(@0 - (@2 + @1)/2.0)",RooArgList(x,low,high) );

  //Naming the vars used in the fit
  RooRealVar* mu     = new RooRealVar("mu",    "peak g",  110, 100, 180);//   0, -1.0, 1.0); mu -> setError(1);
  RooRealVar* sig    = new RooRealVar("sigma", "width g",  10,   0,  80); // 1,   0.,   5); sig -> setError(1);
  RooRealVar* lambda = new RooRealVar("lambda","lambda",  0.5,  -5,  10); //0.5,   0.,   5); lambda -> setError(1);
  EXModGaus* exmg_model = new EXModGaus(name,"Analytic Exp. (N=1)", x, *mu, *sig, *lambda, low_d, high_d);
  exmg_model -> fitTo(hist_to_fit,SumW2Error(kTRUE) ,Save(kTRUE));
  return exmg_model;
}

AsymGenGaussian* fit_agg_model(RooRealVar &x, RooDataHist &hist_to_fit, const char *name, double low_d=100, double high_d=180){
  //Modified Gaussian fitting function
  //new range is -1 to 1
  double range_d = 2;  
  x.setRange("Range1",100,br_left);
  x.setRange("Range2",br_right,180);
 

  //Naming the vars used in the fit
  RooRealVar* kappa = new RooRealVar("kappa", "kappa", -0.5,  -2,    2);
  RooRealVar* alpha = new RooRealVar("alpha", "alpha",   10,   0, 1000);
  RooRealVar* xsi   = new RooRealVar("xsi",   "xsi",    110,  50,  180);
  //Currently this value is not used but is part of the constructor
 // RooRealVar* sigma = new RooRealVar("sigma", "sigma",    10,    0,  60); //0.5,   0.,   5); lambda -> setError(1);

  AsymGenGaussian*  agg_model = new AsymGenGaussian(name,"Asymm. Gen. Gaussian", x, *kappa, *alpha, *xsi, low_d, high_d);
  agg_model -> fitTo(hist_to_fit,SumW2Error(kTRUE) ,Save(kTRUE), Range("Range1,Range2"));
  //agg_model -> fitTo(hist_to_fit,SumW2Error(kTRUE) ,Save(kTRUE), Range("Range2"));

  return agg_model;
}

RooGamma* fit_gam_model(RooRealVar &x, RooDataHist &hist_to_fit, const char *name, double low_d=100, double high_d=180){
  //Naming the vars used in the fit
  RooRealVar* mu     = new RooRealVar("mu",     "peak g",  110,  0, 180);  //  0, -1.0, 1.0); mu -> setError(1);
  RooRealVar* gamma  = new RooRealVar("gamma",   "gamma",  5.0,  0, 1000); //0.5,   0.,   5); gamma -> setError(1);
  RooRealVar* beta   = new RooRealVar("beta",     "beta",  10.0, 0, 1000); //0.5,   0.,   5); beta -> setError(1);
  RooGamma* gamma_model = new RooGamma(name,"gamma_func", x, *gamma, *beta, *mu);//*mu, *sig, *lambda, low_d, high_d);
  gamma_model -> fitTo(hist_to_fit,SumW2Error(kTRUE) ,Save(kTRUE));
  return gamma_model;
}

RooFFTConvPdf* fit_pow3_model(RooRealVar &x, RooDataHist &hist_to_fit, const char *name){
  RooRealVar  *mux    = new RooRealVar( "mux",    "mean of gaussian x", 0);// 110, 0,180);
  RooRealVar  *sigxp3 = new RooRealVar( "sigxp3", "width of gaussian x pow3", 4.47, 1., 10.); sigxp3 -> setError(0.1);
  RooRealVar  *p3t    = new RooRealVar( "t3",     "t pow3", 106, 100., 180.);                 p3t -> setError(1.0);
  RooRealVar  *p3f1   = new RooRealVar( "f13",    "f1 pow3", 0.5, 0., 10.);                   p3f1 -> setError(0.005);
  RooRealVar  *p3p1   = new RooRealVar( "p13",    "p1 pow3", -7., -20., -1.1);                p3p1 -> setError(0.1);
  RooRealVar  *p3p2   = new RooRealVar( "p23",    "p2 pow3", -6.8, -20., -1.1);               p3p2 -> setError(0.1);
  RooGaussian *gxp3   = new RooGaussian("gxp3",   "gaussian PDF pow3", x, *mux, *sigxp3);


  //Two forms of POW3
  RooGenericPdf *st3 = new RooGenericPdf("st3", "step3", "( ((@0-@1)*153.846<0.0) ? 0.0 : (((@0-@1)*153.846 >1.0) ? 1.0 : ((@0-@1)*153.846) ) )*((@3)/((170^(1+@2) - (@1)^(1+@2))/(1+@2))*(@0)^(@2) + (1. - @3)/((170^(1+@4) - (@1)^(1+@4))/(1+@4))*(@0)^(@4))", RooArgList(x,*p3t,*p3p1,*p3f1,*p3p2));
  //RooGenericPdf *st3 = new RooGenericPdf("st3", "step3", "((@3)/((170^(1+@2) - (@1)^(1+@2))/(1+@2))*(@0)^(@2) + (1. - @3)/((170^(1+@4) - (@1)^(1+@4))/(1+@4))*(@0)^(@4))", RooArgList(x,*p3t,*p3p1,*p3f1,*p3p2));
  RooFFTConvPdf *p3_model = new RooFFTConvPdf(name, "Power Law (N=3)", x, *st3, *gxp3);
  p3_model -> setBufferFraction(0.5);
  p3_model -> fitTo(hist_to_fit,SumW2Error(kTRUE) ,Save(kTRUE));

  return p3_model;
}

//RooPoisson* fit_poisson_model(){
//
//}


RooFFTConvPdf* fit_pow1_model(RooRealVar &x, RooDataHist &hist_to_fit, const char *name){
  RooRealVar  *mux        = new RooRealVar("mux",        "mean of gaussian x", 0);
  RooRealVar  *sigxpow1   = new RooRealVar("sigmaxpow1", "width of gaussian x pow1", 7.1179, 3., 10); sigxpow1 -> setError(0.1);
  RooRealVar  *pow1t      = new RooRealVar("t1",  "t pow1",  105.16, 100., 120.);  pow1t -> setError(1.0);
  RooRealVar  *pow1f1     = new RooRealVar("f11", "f1 pow1", 0.26223,   0.,   6.);
  RooRealVar  *pow1p1     = new RooRealVar("p11", "p1 pow1",    -6.5, -12.,  -4.); 
  RooRealVar  *pow1offset = new RooRealVar("offset", "offset pow1", 32., 10., 50.);
  RooGaussian *gaussxpow1 = new RooGaussian("gaussxpow1", "gaussian PDF pow1", x, *mux, *sigxpow1);

  //Power one fitting function
  RooGenericPdf *step = new RooGenericPdf("step", "step", "(((@0-@1)*153.846<0.0) ? 0.0 : (((@0-@1)*153.846 >1.0) ? 1.0 : ((@0-@1)*153.846) )) * (@0^@2)", RooArgList(x,*pow1t,*pow1p1));
  RooFFTConvPdf *pow1_model = new RooFFTConvPdf(name, "Power Law (N=1)", x, *step, *gaussxpow1);
  pow1_model -> setBufferFraction(0.5);
  //RooFitResult* pow1_fit = pow1_model.fitTo(hist_total_rdh,SumW2Error(kTRUE) ,Save(kTRUE));
  pow1_model -> fitTo(hist_to_fit,SumW2Error(kTRUE) ,Save(kTRUE));

  return pow1_model;
}

//Exponential function of order 1
RooFFTConvPdf* fit_exp1_model(RooRealVar &x, RooDataHist &hist_to_fit, const char *name){
  RooRealVar *mux        = new RooRealVar("mux",        "mean of gaussian x", 0);
  RooRealVar *sigma_exp1 = new RooRealVar("sigma_exp1", "sigma_exp1",    7.,  1., 15.);
  RooRealVar *exp1_t     = new RooRealVar("exp_t1",     "t exp1",    106.55, 95., 115.);
  RooRealVar *exp1_p1    = new RooRealVar("exp_p1",     "p1 exp1",    -0.02, -2., 0.0);

  RooGenericPdf *step_exp1 = new RooGenericPdf("step_exp1", "step_exp1", "( ((@0-@1)*153.85<0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85) ) )*(TMath::Exp(@0*@2))", 
                                               RooArgList(x,*exp1_t,*exp1_p1));
  RooGaussian *gauss_exp1  = new RooGaussian("gaussxexp1", "gaussian PDF exp1", x, *mux,  *sigma_exp1);

  RooFFTConvPdf *exp1_model = new RooFFTConvPdf(name, "Num. Exp (N=1)", x, *step_exp1, *gauss_exp1);
  exp1_model -> setBufferFraction(0.5);
  RooFitResult* exp1_fit = exp1_model -> fitTo(hist_to_fit,SumW2Error(kTRUE) ,Save(kTRUE));

  return exp1_model;
}

//Landau 1 function
RooFFTConvPdf* fit_lau1_model(RooRealVar &x, RooDataHist &hist_to_fit, const char *name){
  RooRealVar *mux        = new RooRealVar("mux",        "mu", 0);
  RooRealVar *sigma_lau1 = new RooRealVar("sigma_lau1", "sigma",   7.,  1., 15.);
  RooRealVar *lau1_t     = new RooRealVar("lau1_t1",    "t lau1", 100., 95., 115.);
  RooRealVar *lau1_p1    = new RooRealVar("lau1_p1",    "p1 lau1",  0.5,  0., 1.);
  RooGenericPdf *step_lau1  = new RooGenericPdf("step_lau1", "step_lau1", "( ((@0-@1)*153.85<0.0) ? 0.0 : (((@0-@1)*153.85 >1.0) ? 1.0 : ((@0-@1)*153.85) ) )*(@2*@0^(-8) + (1-@2)*@0^(-7))",
                                                RooArgList(x,*lau1_t,*lau1_p1));
  RooGaussian *gauss_lau1   = new RooGaussian("gaussxlau1",     "gaussian PDF lau1", x, *mux,      *sigma_lau1);
  RooFFTConvPdf *lau1_model = new RooFFTConvPdf("lau1_model", "Laurent Function (N=1)", x, *step_lau1, *gauss_lau1);

  lau1_model -> setBufferFraction(0.5);
  lau1_model -> fitTo(hist_to_fit,SumW2Error(kTRUE) ,Save(kTRUE));
  return lau1_model;
}


//This function will plot just one function on the histogram with a residual plot
void plot_one_func(string infile, string outfile, string func){
  //Sets a couple of constants used. Maybe shift to input variables with default values
  int nbins = 160;
  double m_lly_low = 100;
  double m_lly_high = 180;

  //Debug 1
  cout << "Top of function." << endl;

  //Creates fit x-axis variable and converts histogram to RooDataHist
  RooRealVar m_lly("m_lly" ,"m_lly", m_lly_low, m_lly_high, "GeV");

  //Debug 1.25
  cout << "..." << endl;


  RooDataHist hist_total_rdh("hist_total_rdh", "m_lly_dh", m_lly,Import(*(get_bkg_hist(infile))));
  //hist_total_rdh_reduce = (hist_total_rdh.reduce((to_string(br_right) + ">m_lly||m_lly<" + to_string(br_left)).c_str()));
  TH1D *sig_hist = get_sig_hist(infile); sig_hist -> Scale(0.1);
  RooDataHist hist_signal_rdh("hist_signal_rdh", "m_lly_sdh", m_lly,Import(*(sig_hist)));
  //delete sig_hist;

  //Debug 1.5
  cout << "Really? Here?" << endl;

  //RooCrystalBall* signal_fit = fit_sig_model(m_lly,hist_signal_rdh,"dscb");
  RooGaussian *signal_fit = fit_sig_gauss(m_lly, hist_signal_rdh,"dscb");


  //Calls whichever function specified by the input variable "func"
  RooFFTConvPdf*   fit_model      = nullptr;
  ModGaus*         fit_model_mg   = nullptr;
  EXModGaus*       fit_model_exmg = nullptr;
  RooGamma*        fit_model_gam  = nullptr;
  AsymGenGaussian* fit_model_agg  = nullptr;
  //RooGenericPdf* fit_model_agg  = nullptr;

  if(func=="pow3"){
    fit_model = fit_pow3_model(m_lly,hist_total_rdh,func.c_str());
  } else if (func=="pow1"){
    fit_model = fit_pow1_model(m_lly,hist_total_rdh,func.c_str());   
  } else if (func=="exp1"){
    fit_model = fit_exp1_model(m_lly,hist_total_rdh,func.c_str());   
  } else if (func=="modg"){
    fit_model_mg = fit_modgauss_model(m_lly,hist_total_rdh,func.c_str());
  } else if (func=="exmg"){
    fit_model_exmg = fit_exmg_model(m_lly,hist_total_rdh,func.c_str(),m_lly_low,m_lly_high);
  } else if (func=="gam"){
    fit_model_gam = fit_gam_model(m_lly,hist_total_rdh,func.c_str());
  } else if (func=="agg"){
    fit_model_agg = fit_agg_model(m_lly,hist_total_rdh,func.c_str(),m_lly_low,m_lly_high);
  } else {
    std::cout << "Could not find fit function specified in input. Please try with a correct function name (e.g. pow3, lau1, or bern4)" << std::endl;
    return;
  }


  RooRealVar *rap_coeff = new RooRealVar("rap_coeff","coeff",1);//0.95,0,1.0);
  RooAddPdf *splusb_fit = nullptr;
  if(func=="modg"){
    splusb_fit = new RooAddPdf("splusb_fit", "splusb_fit", *fit_model_mg, *signal_fit, *rap_coeff);//RooArgList(fit_model,signal_fit)); RooArgList(func.c_str(),"dscb") 
  } else if(func=="exmg"){
    splusb_fit = new RooAddPdf("splusb_fit", "splusb_fit", *fit_model_exmg, *signal_fit, *rap_coeff);//RooArgList(fit_model,signal_fit)); RooArgList(func.c_str(),"dscb") 
  } else if(func=="gam"){
    splusb_fit = new RooAddPdf("splusb_fit", "splusb_fit", *fit_model_gam, *signal_fit, *rap_coeff);//RooArgList(fit_model,signal_fit)); RooArgList(func.c_str(),"dscb") 
  } else if(func=="agg"){
    splusb_fit = new RooAddPdf("splusb_fit", "splusb_fit", *fit_model_agg, *signal_fit, *rap_coeff);
  } else{
    splusb_fit = new RooAddPdf("splusb_fit", "splusb_fit", *fit_model, *signal_fit, *rap_coeff);//RooArgList(fit_model,signal_fit)); RooArgList(func.c_str(),"dscb") 
  }

  //hist_total_rdh.add(hist_signal_rdh);
  //splusb_fit -> fitTo(hist_total_rdh,SumW2Error(kTRUE) ,Save(kTRUE));

  //Debug 2
  cout << "model returned." << endl;

  //Define different plots for the data/fit plot and the residual plot
  //Both will be placed on the same canvas  
  RooPlot* data_and_fit  = m_lly.frame(m_lly_low, m_lly_high);
  hist_total_rdh.plotOn( data_and_fit,Binning(nbins), RooFit::Name("hist_total_rdh_fr"));
  splusb_fit -> plotOn( data_and_fit, RooFit::Components("splusb_fit"), LineStyle(kSolid),  LineColor(TColor::GetColor("#0000ff")), RooFit::Name("splusb_fit_fr"));
  splusb_fit -> plotOn( data_and_fit, RooFit::Components(func.c_str()), LineStyle(kDashed), LineColor(TColor::GetColor("#0000ff")), RooFit::Name((func+"_fr").c_str()));
  splusb_fit -> plotOn( data_and_fit, RooFit::Components("dscb"),       LineStyle(kDashed), LineColor(TColor::GetColor("#ff0000")), RooFit::Name("dscb_fr"));
  //signal_fit -> plotOn(data_and_fit, RooFit::Components("dscb"),     LineStyle(kSolid), LineColor(TColor::GetColor("#ff0000")), RooFit::Name("dscb_fr2"));
  //fit_model -> plotOn( data_and_fit, RooFit::Components(func.c_str()), LineStyle(kDashed), LineColor(TColor::GetColor("#df9f1f")), RooFit::Name(func.c_str()));
  //hist_signal_rdh.plotOn(data_and_fit,Binning(nbins), RooFit::Name("hist_signal_rdh"));
  //fit_model_exmg -> plotOn( data_and_fit, RooFit::Components(func.c_str()), LineStyle(kDashed), LineColor(TColor::GetColor("#df9f1f")), RooFit::Name(func.c_str()));


  //Debug 3
  cout << "data_and_fit plot made." << endl;


  //Residual plot created from the data_and_fit frame
  RooHist *residual_rdh = data_and_fit -> residHist("hist_total_rdh_fr","splusb_fit_fr");//func.c_str()); 

  //Debug 4
  cout << "Residual plot made." << endl;


  TString pdf= outfile + ".pdf";
  TCanvas *canvas_plot     = new TCanvas("canvas_plot", "outplot", 1920, 1440);
  TLegend *legend          = new TLegend(0.6,0.8,0.9,0.9);
  TLine   *resid_plot_line = new TLine(m_lly_low,0,m_lly_high,0);

  //Add only entry to legend, make transparent?
  legend->AddEntry("hist_total_rdh_fr"  , "Bkg. MC");//,"L");
  legend->AddEntry("splusb_fit_fr"      , "S + B Fit", "L");//,"L");
  legend->AddEntry((func+"_fr").c_str() , "Bkg. Fit", "L");//,"L");
  legend->AddEntry("dscb_fr"            , "Sig. Fit", "L");//,"L");
  legend->SetFillStyle(0);
  legend->SetLineColorAlpha(kWhite,0);

  //Divide canvas into two pieces
  canvas_plot -> Divide(1,2,0.0,0.0);

  //Set minimum number of events to 0 (quality of life change)
  data_and_fit -> SetMinimum(0.0);

  //A couple of options that should make plotting look nice
  canvas_plot -> Update();
  gStyle      -> SetOptFit();
  
  cout << "plotting on canvas here" << endl;
  double margin_lr  = 0.1;
  double margin_top = 0.05;
  double margin_bot = 0.1;

  double split_point = 0.35;
  double min_x = 0.05;
  double max_x = 1.00;
  double min_y = 0.15;
  
  
  //plot upper plot in upper 7/10ths of plot
  canvas_plot -> cd(1);
  gPad        -> SetPad(min_x, split_point, max_x, 1.0);
  gPad        -> SetMargin(margin_lr, margin_lr, 0.0, margin_top);

  data_and_fit -> GetXaxis() -> SetLabelSize(0);
  data_and_fit -> SetTitle("");
  data_and_fit -> Draw();
  legend       -> Draw();


  //plot lower plot in lower 3/10ths of plot
  canvas_plot -> cd(2);
  gPad        -> SetPad(min_x, min_y,  1, split_point);
  gPad        -> SetMargin(margin_lr, margin_lr, margin_bot, 0);

  //residual_rdh.Draw();
  //residual_rdh -> SetTitle("");
  //residual_rdh -> GetXaxis() -> SetTitle("m_{ll#gamma}");
  residual_rdh -> SetTitle(";m_{ll#gamma};data - fit");
  residual_rdh -> GetXaxis() -> SetTitleSize(0.15);
  residual_rdh -> GetXaxis() -> SetLabelSize(0.1);
  residual_rdh -> GetXaxis() -> SetLimits(m_lly_low, m_lly_high);
  residual_rdh -> GetXaxis() -> SetTitleOffset(1);

  //residual_rdh -> GetYaxis() -> SetTitle("data - fit");
  residual_rdh -> GetYaxis() -> SetTitleSize(0.15);
  residual_rdh -> GetYaxis() -> SetLabelSize(0.10);
  residual_rdh -> GetYaxis() -> SetTickLength(0.01);
  residual_rdh -> GetYaxis() -> SetTitleOffset(1);
  residual_rdh -> Draw("AP");

  resid_plot_line -> SetLineStyle(kDashed);
  resid_plot_line -> SetLineColorAlpha(kBlue,0.75);
  resid_plot_line -> Draw("SAME");


  //Debug 4
  cout << "File plotted, will save and return." << endl;

  canvas_plot -> Print(pdf);

  return;
}



void fit_and_plot(string infile,string outfile){
  //Opens file with signal and bkg histogram
  TFile *file = new TFile(infile.c_str());

  //Defines variables to pull hists from files
  TCanvas* canvas_hist = static_cast<TCanvas*>(file->Get("canvas"));
  TH1D *hist_bkg = nullptr;
  TH1D *hist_sig = nullptr;
  //Finds histogram in file
  for(int idx_h=0; idx_h<50000; idx_h++){

    //Checks if the histogram is a nullptr and looks for each histogram in the file with the particular number value after it
    if(hist_sig == nullptr){
      hist_sig = static_cast<TH1D*>( canvas_hist->FindObject( ("sig_Signal_"+std::to_string(idx_h)).c_str() ) );
    }
    if(hist_bkg == nullptr){
      hist_bkg = static_cast<TH1D*>( canvas_hist->FindObject( ("bkg_Background_"+std::to_string(idx_h)).c_str() ) );
    }

    if(hist_sig != nullptr && hist_bkg != nullptr){break;}
  }

  //Scales signal to appropriate height
  hist_sig -> Scale(0.1);

  //Creates total S + B histogram. Right now just want B
  TH1D *hist_total = static_cast<TH1D*>(hist_bkg -> Clone());
  //hist_total -> Add(hist_sig); 

  //Scales hist and rebins for plotting later
  hist_sig -> Scale(10.0);
  hist_sig = static_cast<TH1D*>(hist_sig -> Rebin(4));


  //Creates fit x-axis variable and converts histogram to RooDataHist
  RooRealVar m_lly("m_lly" ,"m_lly", 100.0, 180.0, "GeV");
  RooDataHist hist_total_rdh("hist_total_rdh", "m_lly_dh",m_lly ,Import(*hist_total));
//  RooDataHist hist_sig_rdh("hist_sig_rdh", "m_lly_dh",m_lly ,Import(*hist_sig));


  RooFFTConvPdf*   fit_model_exp1 = fit_exp1_model(    m_lly, hist_total_rdh, "exp1_model");
  RooFFTConvPdf*   fit_model_pow1 = fit_pow1_model(    m_lly, hist_total_rdh, "pow1_model");   
  RooFFTConvPdf*   fit_model_pow3 = fit_pow3_model(    m_lly, hist_total_rdh, "pow3_model");
  //RooFFTConvPdf*   fit_model_lau1 = fit_lau1_model(    m_lly, hist_total_rdh, "lau1_model");
  ModGaus*         fit_model_mg   = fit_modgauss_model(m_lly, hist_total_rdh, "mg_model");
  EXModGaus*       fit_model_exmg = fit_exmg_model(    m_lly, hist_total_rdh, "exmg_model", 100, 180);
  RooGamma*        fit_model_gam  = fit_gam_model(     m_lly, hist_total_rdh, "gam_model");
  AsymGenGaussian* fit_model_agg  = fit_agg_model(     m_lly, hist_total_rdh, "agg_model");


  RooPlot* xframe4 = m_lly.frame(100,180);
  hist_total_rdh.plotOn(xframe4,Binning(80), RooFit::Name("hist_total_rdh"));
  //hist_sig_rdh.plotOn(xframe4,Binning(80), RooFit::Name("hist_sig_rdh"));

  fit_model_pow3 -> plotOn( xframe4, RooFit::Components("pow3_model"), LineStyle(kSolid), LineColor(TColor::GetColor("#df9f1f")), RooFit::Name("pow3_model_fr"));
  fit_model_pow1 -> plotOn( xframe4, RooFit::Components("pow1_model"), LineStyle(kSolid), LineColor(TColor::GetColor("#eb0161")), RooFit::Name("pow1_model_fr"));
  fit_model_exp1 -> plotOn( xframe4, RooFit::Components("exp1_model"), LineStyle(kSolid), LineColor(TColor::GetColor("#3be363")), RooFit::Name("exp1_model_fr"));
  //fit_model_lau1 -> plotOn( xframe4, RooFit::Components("lau1_model"), LineStyle(kSolid), LineColor(TColor::GetColor("#f1b83e")), RooFit::Name("lau1_model"));
  fit_model_exmg -> plotOn( xframe4, RooFit::Components("exmg_model"), LineStyle(kSolid), LineColor(TColor::GetColor("#fbff00")), RooFit::Name("exmg_model_fr"));
  fit_model_agg  -> plotOn( xframe4, RooFit::Components("agg_model"),  LineStyle(kSolid), LineColor(TColor::GetColor("#6edbff")), RooFit::Name("agg_model_fr"));

  //fit_model_gam  -> plotOn( xframe4, RooFit::Components("gam_model"),  LineStyle(kDashed), LineColor(TColor::GetColor("#f00078")), RooFit::Name("gam_model"));
  //fit_model_mg   -> plotOn( xframe4, RooFit::Components("mg_model"),   LineStyle(kDashed), LineColor(TColor::GetColor("#a800f0")), RooFit::Name("mg_model"));

  //fit_model_modgauss -> plotOn(xframe4, RooFit::Components("modgauss_model"), LineStyle(kDashed), LineColor(TColor::GetColor("#313dfd")), Name("modgauss_model"));
  //fit_ber2_model.plotOn( xframe4, RooFit::Components("bern2_model"),     LineStyle(kDashed), LineColor(TColor::GetColor("#cd1b10")), Name("bern2_model"));
  //fit_ber3_model.plotOn( xframe4, RooFit::Components("bern3_model"),     LineStyle(kDashed), LineColor(TColor::GetColor("#b400d4")), Name("bern3_model"));
  //fit_ber4_model.plotOn( xframe4, RooFit::Components("bern4_model"),     LineStyle(kDashed), LineColor(TColor::GetColor("#22bbc5")), Name("bern4_model"));
  //fit_ber5_model.plotOn( xframe4, RooFit::Components("bern5_model"),     LineStyle(kDashed), LineColor(TColor::GetColor("#009c26")), Name("bern5_model"));

  TString pdf= outfile + ".pdf";
  TCanvas *canvas_plot = new TCanvas("c1", "outplot", 1920, 1080);
  //std::unique_ptr<TLegend> legend = xframe4 -> BuildLegend();
  
  
  TLegend *legend = new TLegend(0.6,0.6,0.9,0.9);
  legend -> AddEntry(xframe4 -> findObject("hist_total_rdh"), "Background MC");
  legend -> AddEntry(xframe4 -> findObject("pow3_model_fr"),  "power 3 polynomial",  "L");
  legend -> AddEntry(xframe4 -> findObject("pow1_model_fr"),  "power 1 polynomial",  "L");
  legend -> AddEntry(xframe4 -> findObject("exp1_model_fr"),  "exp. 1 function",     "L");
  //legend -> AddEntry(xframe4 -> findObject("lau1_model_fr"),  "Landau 1 function",   "L");
  legend -> AddEntry(xframe4 -> findObject("exmg_model_fr"),  "ex. mod. Gaussian",   "L");
  legend -> AddEntry(xframe4 -> findObject("agg_model_fr"),   "Asym. Gen. Gaussian", "L");

/*  legend->AddEntry("pow3_model", "power 3 polynomial","L");
  legend->AddEntry("pow1_model", "power 1 polynomial","L");
  legend->AddEntry("exp1_model", "exp. 1 function",   "L");
  legend->AddEntry("lau1_model", "Landau 1 function", "L");
  legend->AddEntry("exmg_model", "ex. mod. Gaussian", "L");
*/  //legend->AddEntry("gam1_model", "gamma function",    "L");
  //legend->AddEntry("mg1_model" , "mod. gaussian",     "L");

  //legend->AddEntry("modgauss_model","modified gaussian","L");
  //legend->AddEntry("bern2_model","Bernstein 2 polynomial","L");
  //legend->AddEntry("bern3_model","Bernstein 3 polynomial","L");
  //legend->AddEntry("bern4_model","Bernstein 4 polynomial","L");
  //legend->AddEntry("bern5_model","Bernstein 5 polynomial","L");



  xframe4 -> SetMinimum(0.0);

  canvas_plot -> Update();
  gStyle  -> SetOptFit();
  //xframe4 -> getAttText()->SetTextSize(0.03);
  xframe4 -> GetXaxis() -> SetTitle("m_{ll#gamma}"); 
  xframe4 -> Draw();
  //hist_sig -> Draw("HIST SAME");
  legend -> Draw();
  canvas_plot -> Print(pdf);

  delete legend;
  delete canvas_plot;
  delete xframe4;
//  delete ber5_fit; delete ber4_fit; delete ber3_fit; delete ber2_fit;
  delete fit_model_mg;   delete fit_model_exmg; delete fit_model_gam;
  delete fit_model_exp1; delete fit_model_pow1; delete fit_model_pow3; //delete fit_model_lau1;
}


std::vector<std::string> plot_list(std::string name_segment_one, std::string name_segment_two, bool all=true){
  //Loops to cover
  std::vector<std::string> flavor   = {"_ee","_mumu"};
  std::vector<std::string> topology = {"ggF", "VBF", "WH_3l", "ZH_MET", "ttH_had", "ttH_lep"};

  if(!all){
    flavor   = {"_ee","_mumu"};
    topology = {"WH_3l", "ZH_MET", "ttH_had", "ttH_lep"};
  }

  std::vector<std::string> list_o_hists = {};
  for(std::string flav : flavor){
    for(std::string top : topology){
      list_o_hists.push_back(name_segment_one + top + flav + name_segment_two); 
    }
  }
  return list_o_hists;
}



//Main function 
void assoc_prod_fitting(){
  //String with paths and plots of interest
  std::string file_path = "/net/cms27/cms27r0/abarzdukas/drawPico/gitChangesDirectory/AfterMichaelAugustPush/draw_pico/plots/run3_cat_mlly/";
  std::string name_seg_one = "cat_";
  std::string name_seg_two = "_ph_all_lly_m__mlly__wgt_1invfb__lumi_nonorm_lin";
  
  //test fitting code with one plot
  std::string test_file = "cat_VBF_ee_ph_all_lly_m__mlly__wgt_1invfb__lumi_nonorm_lin.root";
  //std::string test_file = "cat_ttH_lep_ee_ph_all_lly_m__mlly__wgt_1invfb__lumi_nonorm_lin.root";

  fit_and_plot(file_path + test_file, "output/test_plot_VBF_ee.pdf");
  //plot_one_func(file_path + test_file, "output/test_plot_VBF_ee.pdf", "pow3");
  //plot_one_func(file_path + test_file, "output/test_plot", "exmg");
  plot_one_func(file_path + test_file, "output/test_plot", "agg");



  //Create list of all the plots
  std::vector<std::string> list_of_plots = plot_list(name_seg_one, name_seg_two);//,false);
  //for(string plot_name : list_of_plots){ fit_and_plot(file_path + plot_name + ".root", "output/" + plot_name); }
  //for(string plot_name : list_of_plots){ plot_one_func(file_path + plot_name + ".root", "output/exp1_fit_" + plot_name, "exp1"); }
  //for(string plot_name : list_of_plots){ plot_one_func(file_path + plot_name + ".root", "output/pow3_fit_" + plot_name, "pow3"); }

  //for(string plot_name : list_of_plots){ plot_one_func(file_path + plot_name + ".root", "output/modg_fit_" + plot_name, "modg"); }
  //for(string plot_name : list_of_plots){ plot_one_func(file_path + plot_name + ".root", "output/exmg_fit_" + plot_name, "exmg"); }
  //for(string plot_name : list_of_plots){ plot_one_func(file_path + plot_name + ".root", "output/gam_fit_" + plot_name, "gam"); }
  //for(string plot_name : list_of_plots){ plot_one_func(file_path + plot_name + ".root", "output/agg_fit_" + plot_name, "agg"); }

}
