#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
//#include "../src/ModGaus.cxx"
// #include "../src/Run2FuncPow1.cxx"
//#include "../src/Run2FuncPow3.cxx"
//#include "../src/SigFunc.cxx"
#include "default.h"
//#include "../src/RooGaussStepBernstein.h"
#include "HZGRooPdfs.h"
#include "RooCrystalBall.h"
//#include "../src/ModGaus.cxx"

using namespace RooFit;
/*
Int_t SumW2Cov(RooNLLVar& nll, RooMinimizer& somemini){
    
    // nll.applyWeightSquared(false);
    // nll.enableOffsetting(true);
    // somemini.setOffsetting(false);
    // somemini.setOffsetting(true);
    // somemini.hesse();
    auto r2 = somemini.save();
    auto cov2 = r2->covarianceMatrix();
    // somemini.setOffsetting(false);
    // somemini.setOffsetting(true);
    nll.applyWeightSquared(true);
    somemini.hesse();
    auto r1 = somemini.save();
    auto cov1 = r1->covarianceMatrix();

    nll.applyWeightSquared(false);
    somemini.hesse();
    ROOT::Math::CholeskyDecompGenDim<double> decomp(cov1.GetNrows(), cov1);
    decomp.Invert(cov1);
    for (int i = 0; i < cov1.GetNrows(); ++i) {
      for (int j = 0; j < i; ++j) {
         cov1(j, i) = cov1(i, j);
      }
   }

    // std::cout << "V = " <<std::endl;
    // cov2.Print();
    // std::cout << "C^-1 = " <<std::endl;
    // cov1.Print();
    // std::cout << std::endl;
    auto invconv =cov1.Similarity(cov2);
    somemini.applyCovarianceMatrix(invconv);
    // std::cout << "no weight NLL = " << r2->covQual() << std::endl << "weighted NLL = " << r1->covQual() << std::endl;
    // res->setCovQual(std::min(r1->covQual(), r2->covQual()));
    // std::cout << "VC^-1V = " <<std::endl;
    // invconv.Print();
    std::cout << std::endl;
    return (std::min(r1->covQual(), r2->covQual()));

}
*/


void range_signi_opt(double lowrange){
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ;

    double xlow = lowrange;
    double xhigh = xlow + 65.;
    int N = 4 * (xhigh - xlow);
    RooRealVar x("x", "m_llg", xlow, xhigh);
    RooRealVar y("y", "photon pT", 15., 80.);
    RooRealVar w("w", "w", -40, 40);
    RooRealVar bdt("bdt", "bdt", -1, 1);
    RooRealVar year("year", "year", 0, 10000);
    RooRealVar lep("lep", "lep", 0, 1);
    RooRealVar ph_eta("ph_eta", "ph_eta", -3, 3);
    RooRealVar nlep("nlep", "nlep", 0, 10);
    RooRealVar njet("njet", "njet", 0, 10);

    RooRealVar nbkg("nbkg","background events",0,1600000);

//Read dataset files
    // RooDataSet *sig_data = RooDataSet::read("data/FullSig_deathvalley_v3_untagged.dat",RooArgList(x, y, bdt, w, year, lep, ph_eta, nlep, njet));
    // RooDataSet *bkg_data = RooDataSet::read("data/SMZg_deathvalley_v3_untagged.dat,data/DY_deathvalley_v3_untagged.dat",RooArgList(x, y, bdt, w, year, lep, ph_eta, nlep, njet));
    // RooDataSet *tot_data = RooDataSet::read("data/SMZg_deathvalley_v3_untagged.dat,data/DY_deathvalley_v3_untagged.dat,data/FullSig_deathvalley_v3_untagged.dat",RooArgList(x, y, bdt, w, year, lep, ph_eta, nlep, njet));

    string sigsamples = "../CMSSW_11_3_4/src/HToZg_combine/rui_datasamples/Signal_untagged4_ratio_extended.dat,../CMSSW_11_3_4/src/HToZg_combine/rui_datasamples/Signal_untagged3_ratio_extended.dat,../CMSSW_11_3_4/src/HToZg_combine/rui_datasamples/Signal_untagged2_ratio_extended.dat,../CMSSW_11_3_4/src/HToZg_combine/rui_datasamples/Signal_untagged1_ratio_extended.dat";
    string bkgsamples = "../CMSSW_11_3_4/src/HToZg_combine/rui_datasamples/DY_untagged1_ratio_extended.dat,../CMSSW_11_3_4/src/HToZg_combine/rui_datasamples/DY_untagged2_ratio_extended.dat,../CMSSW_11_3_4/src/HToZg_combine/rui_datasamples/DY_untagged3_ratio_extended.dat,../CMSSW_11_3_4/src/HToZg_combine/rui_datasamples/DY_untagged4_ratio_extended.dat,../CMSSW_11_3_4/src/HToZg_combine/rui_datasamples/SMZg_untagged1_ratio_extended.dat,../CMSSW_11_3_4/src/HToZg_combine/rui_datasamples/SMZg_untagged2_ratio_extended.dat,../CMSSW_11_3_4/src/HToZg_combine/rui_datasamples/SMZg_untagged3_ratio_extended.dat,../CMSSW_11_3_4/src/HToZg_combine/rui_datasamples/SMZg_untagged4_ratio_extended.dat";


    //string sigsamples = "../CMSSW_11_3_4/src/HToZg_combine/rui_datasamples/Signal_untagged4_run2.dat,../CMSSW_11_3_4/src/HToZg_combine/rui_datasamples/Signal_untagged3_run2.dat,../CMSSW_11_3_4/src/HToZg_combine/rui_datasamples/Signal_untagged2_run2.dat,../CMSSW_11_3_4/src/HToZg_combine/rui_datasamples/Signal_untagged1_run2.dat";
    //string bkgsamples = "../CMSSW_11_3_4/src/HToZg_combine/rui_datasamples/DY_untagged1_run2.dat,../CMSSW_11_3_4/src/HToZg_combine/rui_datasamples/DY_untagged2_run2.dat,../CMSSW_11_3_4/src/HToZg_combine/rui_datasamples/DY_untagged3_run2.dat,../CMSSW_11_3_4/src/HToZg_combine/rui_datasamples/DY_untagged4_run2.dat,../CMSSW_11_3_4/src/HToZg_combine/rui_datasamples/SMZg_untagged1_run2.dat,../CMSSW_11_3_4/src/HToZg_combine/rui_datasamples/SMZg_untagged2_run2.dat,../CMSSW_11_3_4/src/HToZg_combine/rui_datasamples/SMZg_untagged3_run2.dat,../CMSSW_11_3_4/src/HToZg_combine/rui_datasamples/SMZg_untagged4_run2.dat";

    RooDataSet *sig_data = RooDataSet::read(sigsamples.c_str(),RooArgList(x, y, bdt, w, year, lep, ph_eta, nlep, njet));
    RooDataSet *bkg_data = RooDataSet::read(bkgsamples.c_str(),RooArgList(x, y, bdt, w, year, lep, ph_eta, nlep, njet));
    RooDataSet *tot_data = RooDataSet::read((sigsamples + "," + bkgsamples).c_str(), RooArgList(x, y, bdt, w, year, lep, ph_eta, nlep, njet));

    string samplename1 = "untag1_rui";
    string samplename2 = "untag2_rui";
    string samplename3 = "untag3_rui";
    string samplename4 = "untag4_rui";

    RooDataSet u1_sig("u1_sig", "u1_sig", sig_data, RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), MVAbins(samplename1), "w");
    RooDataSet u2_sig("u2_sig", "u2_sig", sig_data, RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), MVAbins(samplename2), "w");
    RooDataSet u3_sig("u3_sig", "u3_sig", sig_data, RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), MVAbins(samplename3), "w");
    RooDataSet u4_sig("u4_sig", "u4_sig", sig_data, RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), MVAbins(samplename4), "w");

    RooDataSet u1_bkg("u1_bkg", "u1_bkg", bkg_data, RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), MVAbins(samplename1), "w");
    RooDataSet u2_bkg("u2_bkg", "u2_bkg", bkg_data, RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), MVAbins(samplename2), "w");
    RooDataSet u3_bkg("u3_bkg", "u3_bkg", bkg_data, RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), MVAbins(samplename3), "w");
    RooDataSet u4_bkg("u4_bkg", "u4_bkg", bkg_data, RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), MVAbins(samplename4), "w");

    RooDataSet u1_tot("u1_tot", "u1_tot", tot_data, RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), MVAbins(samplename1), "w");
    RooDataSet u2_tot("u2_tot", "u2_tot", tot_data, RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), MVAbins(samplename2), "w");
    RooDataSet u3_tot("u3_tot", "u3_tot", tot_data, RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), MVAbins(samplename3), "w");
    RooDataSet u4_tot("u4_tot", "u4_tot", tot_data, RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), MVAbins(samplename4), "w");

    


//Create RooDataHist for binned fit
    x.setRange(xlow, xhigh);
    x.setBins(N);

    RooDataHist hist_u1("hist_u1", "hist_u1", x, u1_tot);
    RooDataHist hist_u2("hist_u2", "hist_u2", x, u2_tot);
    RooDataHist hist_u3("hist_u3", "hist_u3", x, u3_tot);
    RooDataHist hist_u4("hist_u4", "hist_u4", x, u4_tot);

    RooDataHist hist_u1_bkg("hist_u1_bkg", "hist_u1_bkg", x, u1_bkg);
    RooDataHist hist_u2_bkg("hist_u2_bkg", "hist_u2_bkg", x, u2_bkg);
    RooDataHist hist_u3_bkg("hist_u3_bkg", "hist_u3_bkg", x, u3_bkg);
    RooDataHist hist_u4_bkg("hist_u4_bkg", "hist_u4_bkg", x, u4_bkg);

    RooDataHist hist_u1_sig("hist_u1_sig", "hist_u1_sig", x, u1_sig);
    RooDataHist hist_u2_sig("hist_u2_sig", "hist_u2_sig", x, u2_sig);
    RooDataHist hist_u3_sig("hist_u3_sig", "hist_u3_sig", x, u3_sig);
    RooDataHist hist_u4_sig("hist_u4_sig", "hist_u4_sig", x, u4_sig);

    double expbkg_u1 = u1_bkg.sumEntries();
    double expbkg_u2 = u2_bkg.sumEntries();
    double expbkg_u3 = u3_bkg.sumEntries();
    double expbkg_u4 = u4_bkg.sumEntries();

    double expsig_u1 = u1_sig.sumEntries();
    double expsig_u2 = u2_sig.sumEntries();
    double expsig_u3 = u3_sig.sumEntries();
    double expsig_u4 = u4_sig.sumEntries();

//POW funcs definition
    RooRealVar meanx("meanx", "mean of gaussian x", 0);//, -30, 30); meanx.setError(1.0);
    //RooRealVar meany("meany", "mean of gaussian y", 26);
    RooRealVar sigmaxpow1("sigmaxpow1", "width of gaussian x pow1", 7.1179, 3., 10); sigmaxpow1.setError(0.1);
    RooRealVar sigmaxpow3("sigmaxpow3", "width of gaussian x pow3", 4.47, 1., 10.); sigmaxpow3.setError(0.1);
    //RooRealVar sigmay("sigmay", "width of gaussian y", 7);
    RooGaussian gaussxpow1("gaussxpow1", "gaussian PDF pow1", x, meanx, sigmaxpow1);
    RooGaussian gaussxpow3("gaussxpow3", "gaussian PDF pow3", x, meanx, sigmaxpow3);
    //RooGaussian gaussy("gaussy", "gaussian PDF", y, meany, sigmay);
    RooRealVar pow1t("t1", "t pow1", 105.16, 100., 120.); pow1t.setError(1.0);
    RooRealVar pow1f1("f11", "f1 pow1", 0.26223, 0., 6.); 
    RooRealVar pow1p1("p11", "p1 pow1", -6.5, -12., -4.); //pow1p1.setError(0.1);
    RooRealVar pow1offset("offset", "offset pow1", 32., 10., 50.); 

    RooRealVar pow3t("t3", "t pow3", 106, 100., 115.); pow3t.setError(1.0);
    RooRealVar pow3f1("f13", "f1 pow3", 0.5, 0., 10.); pow3f1.setError(0.005);
    RooRealVar pow3p1("p13", "p1 pow3", -7., -20., -1.1); pow3p1.setError(0.1);
    RooRealVar pow3p2("p23", "p2 pow3", -6.8, -20., -1.1); pow3p2.setError(0.1);
    //Run2Func run2model = Run2Func("r2step", "r2step", x, f, p1, p2, t, meanx, sigmax);
    RooGenericPdf step("step", "step", "( ((@0-@1)*153.846<0.0) ? 0.0 : (((@0-@1)*153.846 >1.0) ? 1.0 : ((@0-@1)*153.846) ) ) * (@0^@2)", RooArgList(x,pow1t,pow1p1));
    RooGenericPdf pow1_erf("pow1_erf","pow1_erf", "erf((@0-@1)/10.)*(@0^@2)", RooArgList(x,pow1t,pow1p1));

//DSCB signal model
    RooRealVar mu("mu", "mu", 124.57, 120., 130.);
    RooRealVar alphaL("alphaL", "alphaL", 7.6328e-01, 0.001,2.3);
    RooRealVar alphaR("alphaR", "alphaR", 1.5, 0.001, 2.5);
    RooRealVar nL("nL", "nL", 10, 1., 120);
    RooRealVar nR("nR", "nR",50, 1., 120);
    RooRealVar sigmaL("sigmaL", "sigmaL", 1.67, 0.5, 3.);

    RooCrystalBall sig_model("sig_model","sig_model", x, mu, sigmaL, alphaL, nL, alphaR, nR);
    double srange = (xlow <= 115.? 115.:xlow) ;
    x.setRange("signal",srange, 132);

    double bound = (xlow <= 103? 10.:(xlow < 112? 50.:30.));
    double p0_ini = (xlow <= 103? 0.3:(xlow < 112? 10:12.));
//Bern 5
    RooRealVar p0("p0", "p0", p0_ini);//, -bound, bound );
    RooRealVar b5p1("b5p1", "b5p1", 10, -bound, bound );
    RooRealVar b5p2("b5p2", "b5p2", 10, -bound, bound);
    RooRealVar b5p3("b5p3", "b5p3",10,-bound,bound);
    RooRealVar b5p4("b5p4", "b5p4", 3,-bound,bound);
    RooRealVar b5p5("b5p5", "b5p5", 3,-bound,bound);
   
//Bern5 102 specific setup
    /*
    RooRealVar b5p1("b5p1_u3", "b5p1_u3", 9.9783e+01, 1.,105.);
    RooRealVar b5p2("b5p2_u3", "b5p2_u3", -5.4134,-50.,0.);
    RooRealVar b5p3("b5p3_u3", "b5p3_u3", 1.1891e+01,5.,25.);
    RooRealVar b5p4("b5p4_u3", "b5p4_u3", 1.8975e+01,-20.,22.);
    RooRealVar b5p5("b5p5_u3", "b5p5_u3", -4.7498,-20.,0.);
    */
    RooRealVar sigma_bern5("sigma_bern5","sigma_bern5"       ,3,  0.01, 15.);
    RooRealVar stepval_bern5("stepval_bern5", "stepval_bern5", 103, 90., 115.);
    RooGaussStepBernstein bern5_model("bern5_model", "Bernstein 5(X) gauss", x, meanx, sigma_bern5, stepval_bern5, RooArgList(p0,b5p1,b5p2,b5p3, b5p4, b5p5));

//Bern 2
    RooGaussStepBernstein bern2_model("bern2_model", "Bernstein 2(X) gauss", x, meanx, sigma_bern5, stepval_bern5, RooArgList(p0,b5p1,b5p2));

//Bern 3
    RooGaussStepBernstein bern3_model("bern3_model", "Bernstein 3(X) gauss", x, meanx, sigma_bern5, stepval_bern5, RooArgList(p0,b5p1,b5p2, b5p3));

// Bern 4 
    RooGaussStepBernstein bern4_model("bern4_model", "Bernstein 4(X) gauss", x, meanx, sigma_bern5, stepval_bern5, RooArgList(p0,b5p1,b5p2,b5p3, b5p4));

//Two forms of POW3
    //RooGenericPdf step3("step3", "step3", "(@0 > @1)*((@3)/((200^(1-@2) - 75^(1-@2))/(1-@2))*(@0)^(@2) + (1. - @3)/((200^(1-@4) - 70^(1-@4))/(1-@4))*(@0)^(@4))", RooArgList(x,pow3t,pow3p1,pow3f1,pow3p2));
    // RooGenericPdf step3("step3", "step3", "(@0 > @1)*((@3)*(@0)^(@2) + (1. - @3)*(@0)^(@4))", RooArgList(x,pow3t,pow3p1,pow3f1,pow3p2));
    RooGenericPdf step3("step3", "step3", "( ((@0-@1)*153.846<0.0) ? 0.0 : (((@0-@1)*153.846 >1.0) ? 1.0 : ((@0-@1)*153.846) ) )*((@3)/((170^(1+@2) - (@1)^(1+@2))/(1+@2))*(@0)^(@2) + (1. - @3)/((170^(1+@4) - (@1)^(1+@4))/(1+@4))*(@0)^(@4))", RooArgList(x,pow3t,pow3p1,pow3f1,pow3p2));
    RooGenericPdf pow3_erf("pow3_erf", "pow3_erf", "erf((@0-@1)/10.)*(@3 * @0^@2 + (1. - @3)*@0^@4)", RooArgList(x,pow3t,pow3p1,pow3f1,pow3p2));


//POW1 model
    RooFFTConvPdf pow1_model("pow1_model", "step (X) gauss", x, step, gaussxpow1);
    pow1_model.setBufferFraction(0.5);

//Mod gauss
    RooRealVar m0("m_{0}","mass peak value [GeV]" ,120,110,180);  m0.setError(1);//110, 180
    RooRealVar vl("#nu_{L}","low-end power"       , 0,  -5,  15);  vl.setError(0.1);// 0, 15
    RooRealVar vr("#Delta#nu","power range"       ,2, -10,  15);  vr.setError(0.1);//-40, 15
    RooRealVar s0("#sigma_{0}","peak width"       , 3,  1., 30);  s0.setError(1);
    RooRealVar sl("#sigma_{L}","low-end width"    , 5,  0.1, 40);  sl.setError(1);//-10, 35
    RooRealVar sh("#sigma_{H}","high-end width"   , 20, 5,60);  sh.setError(1);//0, 100
    RooRealVar cc("const","const"   , 0);
    //RooRealVar r5("r5","pt quadratic in peak"   ,defaults[10], -1,5);
    

    //RooRealVar a0("a0", "a0 par",97.0, 92.0, 99.0);
    //RooRealVar p("p", "p par",5.0, 3.0, 9.0);
    //x.setRange(xlow, 180);
    //ModGaus mod_model = ModGaus("fx","G_{M}", x,m0,vl,vr,s0,sl,sh, cc, xlow, xhigh);

//Choose POW3 between default and self-norm
    //auto run2modelpow3 = Run2FuncPow3("pow3_model", "pow3 model", x, pow3p1, pow3p2, pow3t, pow3f1, meanx, sigmaxpow3, xlow, xhigh);

    x.setBins(20000, "cache");
    //auto run2modelpow3 = RooFFTConvPdf("run2modelpow3", "step (X) gauss", x, step3, gaussxpow3);
    // auto run2modelpow3 = RooFFTConvPdf("run2modelpow3", "erf (X) gauss", x, pow3_erf, gaussxpow3);
    //run2modelpow3.setBufferFraction(0.5);


    RooFitResult *res;
    std::cout << "n u1 = " <<  expbkg_u1 << std::endl;
    std::cout << "n u1 sig = " <<  expsig_u1 << std::endl;

//Choose fitter setup
    double someNLL, otherNLL;

    sig_model.fitTo( hist_u1_sig, AsymptoticError(true),Save(true), PrintLevel(-1), Strategy(0), Range("signal"));

    TCanvas *cs = new TCanvas("cs", "cs", 400, 400);
    cs->cd();
    RooPlot *xframesig = x.frame(Title(plot_title(samplename1)));
    hist_u1_sig.plotOn(xframesig);
    sig_model.plotOn(xframesig,LineColor(kRed), Name("sigmodel"));
    xframesig->Draw();
    cs->Draw();

    mu.setVal(124.7);
    mu.setConstant(true);
    sigmaL.setConstant(true);
    nL.setConstant(true);
    nR.setConstant(true);
    alphaL.setConstant(true);
    alphaR.setConstant(true);

    auto c11 = RooRealVar("c11", "c11", expbkg_u1, 0., 3.*expbkg_u1);
    auto c21 = RooRealVar("c21", "c21", expsig_u1, -1000., 1000);
    
    RooNLLVar *nll_pre = new RooNLLVar("nll_pre","-log(L)", bern3_model, hist_u1_bkg);
    // nll->enableOffsetting(false);
    RooMinimizer mini_pre(*nll_pre);
    mini_pre.setPrintLevel(-1);
    // mini.setOffsetting(true);
    mini_pre.setEps(100);
    mini_pre.setStrategy(0);
    mini_pre.minimize("Minuit2","migrad");
    mini_pre.setEps(1);
    mini_pre.setStrategy(0);
    mini_pre.minimize("Minuit2","migrad");
    mini_pre.setEps(0.01);
    mini_pre.setPrintLevel(1);
    mini_pre.setStrategy(0);
    mini_pre.setOffsetting(true);
    mini_pre.minimize("Minuit2","migrad");
    mini_pre.hesse();
    res = mini_pre.save();
    // b5p1.setConstant(true);
    // b5p4.setConstant(true);

    /*
     RooNLLVar *nll = new RooNLLVar("nll","-log(L)", bern3_model, hist_u1_bkg);
     RooMinimizer *mini = new RooMinimizer(*nll);  //mini->setPrintLevel(-1);
     // nll->enableOffsetting(true);
     mini->setPrintLevel(1);
     mini->setEps(0.01);
     mini->setStrategy(0);
     // mini->setOffsetting(true);
     mini->minimize("Minuit2","migrad");
     mini->hesse();
     //int qual = SumW2Cov(*nll, *mini);
     res = mini->save();
     //res->setCovQual(qual); */
    res->Print("V");

    RooPlot *xframe0 = x.frame(Title(plot_title(samplename1)));
    hist_u1_bkg.plotOn(xframe0);
    bern3_model.plotOn(xframe0,LineColor(kRed), Name("bern3"));
    TCanvas *c = new TCanvas("c", "c", 400, 400);
    c->cd();
    xframe0->Draw();
    c->Draw();
    c->SaveAs(((string)("../public/significance_plot/MC_fit_bern3_" + std::to_string(lowrange) + "_u1.pdf")).c_str());
    RooChi2Var chi2("chi2", "chi2", bern3_model, hist_u1_bkg, DataError(RooAbsData::SumW2));
    double chi2_mc = chi2.getVal();
    // RooPlot *plott = b5p1.frame(Range(b5p1.getValV() - 0.001, b5p1.getValV() + 0.001), Bins(100), Title("p1")); 
    // TCanvas *c1 = new TCanvas("c1", "c1", 400, 400);
    // c1->cd();
    // nll->plotOn(plott,ShiftToZero());
    // plott->Draw();
    // c1->Draw();
    
    // RooPlot *plot2 = b5p2.frame(Range(b5p2.getValV() - 0.001, b5p2.getValV() + 0.001), Bins(100), Title("p2")); 
    // TCanvas *c2 = new TCanvas("c2", "c2", 400, 400);
    // c2->cd();
    // nll->plotOn(plot2,ShiftToZero());
    // plot2->Draw();
    // c2->Draw();

    // RooPlot *plot3 = b5p3.frame(Range(b5p3.getValV() - 0.001, b5p3.getValV() + 0.001), Bins(100), Title("p3")); 
    // TCanvas *c3 = new TCanvas("c3", "c3", 400, 400);
    // c3->cd();
    // nll->plotOn(plot3,ShiftToZero());
    // plot3->Draw();
    // c3->Draw();

    // RooPlot *plot4 = b5p4.frame(Range(b5p4.getValV() - 0.001, b5p4.getValV() + 0.001), Bins(100), Title("p4")); 
    // TCanvas *c4 = new TCanvas("c4", "c4", 400, 400);
    // c4->cd();
    // nll->plotOn(plot4, ShiftToZero());
    // plot4->Draw();
    // c4->Draw();

    // RooPlot *plot5 = stepval_bern3.frame(Range(stepval_bern3.getValV() - 0.001, stepval_bern3.getValV() + 0.001), Bins(100), Title("step")); 
    // TCanvas *c5 = new TCanvas("c5", "c5", 400, 400);
    // c5->cd();
    // nll->plotOn(plot5, ShiftToZero());
    // plot5->Draw();
    // c5->Draw();

    // RooPlot *plot6 = sigma_bern3.frame(Range(sigma_bern3.getValV() - 0.001, sigma_bern3.getValV() + 0.001), Bins(100), Title("sigma")); 
    // TCanvas *c6 = new TCanvas("c6", "c6", 400, 400);
    // c6->cd();
    // nll->plotOn(plot6, ShiftToZero());
    // plot6->Draw();
    // c6->Draw();

    // c11.setVal(100 * expbkg_u1);
    auto tot_model = RooAddPdf("tot_model", "tot_model", RooArgList(sig_model, bern3_model), RooArgList(c21, c11));
    //     x.setRange(xlow,xhigh);
    //auto toy_sample = tot_model.generate(RooArgSet(x), NumEvents(10000 * (int)(expbkg_u1+expsig_u1)));
    //RooRealVar toyw("toyw", "toyw", 0.0001);
    //toy_sample->addColumn(toyw);
    //RooDataSet w_toy_sample("w_toy_sample", "w_toy_sample", toy_sample, RooArgSet(x, toyw),"", "toyw");
    x.setBins(N);
    RooDataHist* hist_toy = tot_model.generateBinned(RooArgSet(x), (int)(expbkg_u1+expsig_u1), Asimov(true)); 
    //RooDataHist *hist_toy = new RooDataHist("hist_toy", "hist_toy", x, w_toy_sample);
    std::cout << "nume events = " << hist_toy->sumEntries() << std::endl;

    c21.setVal(2.*expsig_u1);
    c21.setConstant(false);
    RooNLLVar *altnll = new RooNLLVar("altnll","-log(L) alt", tot_model, *hist_toy, Extended(true));
    RooMinimizer *altmini = new RooMinimizer(*altnll);  //mini->setPrintLevel(-1);
    // nll->enableOffsetting(true);
    altmini->setEps(100);
    altmini->setStrategy(0);
    altmini->setPrintLevel(-1);
    altmini->setOffsetting(true);
    altmini->minimize("Minuit2","migrad");
    altmini->setEps(1);
    altmini->setStrategy(0);
    altmini->minimize("Minuit2","migrad");
    altmini->setEps(0.01);
    altmini->setStrategy(0);
    altmini->minimize("Minuit2","migrad");
    otherNLL = (double)altnll->getVal();
    altmini->hesse();
    res = altmini->save();

    res->Print("V");


    RooPlot *xframe = x.frame(Title("Float r"));
    hist_toy->plotOn(xframe);
    tot_model.plotOn(xframe, LineColor(kRed));
    TCanvas *c0 = new TCanvas("c0", "c0", 400, 400);
    c0->cd();
    xframe->Draw();
    c0->Draw();
    c0->SaveAs(((string)("../public/significance_plot/Asimov_fit_bern3_" + std::to_string(lowrange) + "_u1.pdf")).c_str());
    RooChi2Var chi2_asimov("chi2_asimov", "chi2_asimov", tot_model, *hist_toy, DataError(RooAbsData::Poisson));
    double chi2_as = chi2_asimov.getVal();


    c21.setVal(0);
    c21.setConstant(true);
    RooNLLVar *nullnll = new RooNLLVar("nullnll","-log(L)", tot_model, *hist_toy, Extended(true));
    RooMinimizer *nullmini = new RooMinimizer(*nullnll);  //mini->setPrintLevel(-1);
    // nll->enableOffsetting(true);
    nullmini->setEps(100);
    nullmini->setPrintLevel(-1);
    nullmini->setStrategy(0);
    nullmini->setOffsetting(true);
    nullmini->minimize("Minuit2","migrad");
    nullmini->setEps(1);
    nullmini->setStrategy(0);
    nullmini->minimize("Minuit2","migrad");
    nullmini->setStrategy(0);
    nullmini->setEps(0.01);
    nullmini->minimize("Minuit2","migrad");
    someNLL = (double)nullnll->getVal();
    nullmini->hesse();
    res = nullmini->save();
    res->Print("V");
    
    RooPlot *xframe1 = x.frame(Title("r = 0 "));
    //hist_toy->plotOn(xframe1, MarkerColor(0));
    //    tot_model.plotOn(xframe1,VisualizeError(*res, 3), FillColor(kBlue));
    tot_model.plotOn(xframe1,VisualizeError(*res, 3), DrawOption("L"), LineWidth(2), LineColor(kBlue));
    tot_model.plotOn(xframe1,LineColor(kRed),  LineWidth(2));
    TCanvas *c1 = new TCanvas("c1", "c1", 400, 400);
    c1->cd();
    xframe1->GetYaxis()->SetRangeUser(0, 0.008);
    xframe1->Draw();
    c1->Draw();
    c1->SaveAs(((string)("../public/significance_plot/Asimov_fit_bern3_" + std::to_string(lowrange) + "_u1_error_band.pdf")).c_str());

    std::cout << "Some NLL = " << someNLL << std::endl << "Other NLL = " << otherNLL << std::endl \
    << "sqrt(2DeltaNLL) = " << ROOT::Math::sqrt(2*(someNLL - otherNLL)) << std::endl;
    std::cout << "MC Chi2 = " << chi2_mc << ", Asimov Chi2 = " << chi2_as << std::endl;
}
