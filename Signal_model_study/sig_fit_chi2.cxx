#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooProdPdf.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooCrystalBall.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include <iostream>
#include <fstream>

using namespace RooFit;
using namespace std; 

int LOWX = 115;
int HIGHX = 131;
int NBIN = 4*(HIGHX - LOWX);

int calculateConvMatrix(RooMinimizer &minimizer, RooNLLVar &nll){
  std::vector<RooNLLVar *> nllComponents;
  std::unique_ptr<RooArgSet> comps{nll.getComponents()};
  for (auto const &arg : *comps) {
    if (RooNLLVar *nllComp = dynamic_cast<RooNLLVar *>(arg)) {
      nllComponents.push_back(nllComp);
    }
  }

  // Calculated corrected errors for weighted likelihood fits
  std::unique_ptr<RooFitResult> rw{minimizer.save()};
  for (auto &comp : nllComponents) {
    comp->applyWeightSquared(true);
  }
  cout << "Xingchen: Calculating sum-of-weights-squared correction matrix for covariance matrix"
                 << std::endl;
  minimizer.hesse();
  std::unique_ptr<RooFitResult> rw2{minimizer.save()};
  for (auto &comp : nllComponents) {
    comp->applyWeightSquared(false);
  }

  // Apply correction matrix
  const TMatrixDSym &matV = rw->covarianceMatrix();
  TMatrixDSym matC = rw2->covarianceMatrix();
  ROOT::Math::CholeskyDecompGenDim<double> decomp(matC.GetNrows(), matC);
  if (!decomp) {
    cout << "Xingchen ERROR: Cannot apply sum-of-weights correction to covariance matrix: correction "
                              "matrix calculated with weight-squared is singular"
                   << std::endl;
    return -1;
  }

  // replace C by its inverse
  decomp.Invert(matC);
  // the class lies about the matrix being symmetric, so fill in the
  // part above the diagonal
  for (int i = 0; i < matC.GetNrows(); ++i) {
    for (int j = 0; j < i; ++j) {
      matC(j, i) = matC(i, j);
    }
  }
  matC.Similarity(matV);
  // C now contiains V C^-1 V
  // Propagate corrected errors to parameters objects
  minimizer.applyCovarianceMatrix(matC);
  
  return std::min(rw->covQual(), rw2->covQual());
}

vector<double> combine_sigma(double s1, double e1, double s2, double e2, double frac){
  vector<double> result;
  result.push_back(frac*s1 + (1-frac)*s2);
  result.push_back(TMath::Sqrt(frac * frac * e1 * e1 + (1-frac)*(1-frac)*e2*e2));
  return result;
}

void printQual(int qual) {
  cout <<endl << "CovMatrix Quality = ";
  switch(qual) {
  case -1 : cout << "Unknown, matrix was externally provided" ; break ;
  case 0  : cout  << "Not calculated at all" ; break ;
  case 1  : cout  << "Approximation only, not accurate" ; break ;
  case 2  : cout  << "Full matrix, but forced positive-definite" ; break ;
  case 3  : cout  << "Full, accurate covariance matrix" ; break ;
  }
  cout << endl;
}

int fit_hist_DSCB(RooDataHist &h_ul, RooRealVar &x, string def_name, TCanvas *c4, TLegend *leg, bool output){
  bool OUTPUT = output;
  string name = "DSCB " + def_name;
  string ulname = name + "";
  RooPlot *xframe_ul = x.frame(Title(ulname.c_str()));
  h_ul.plotOn(xframe_ul, LineColor(kBlue), DataError(RooAbsData::SumW2));
  RooPlot *xframe = x.frame(Title(def_name.c_str()));
  string outputfile = "fit_results/" + name + ".txt";
  if (OUTPUT) freopen (outputfile.c_str(),"w",stdout);

  RooRealVar mu("mu", "mu", 124.57, 120., 130.);
  RooRealVar alphaL("alphaL", "alphaL", 0.77, 0.01, 5);
  RooRealVar alphaR("alphaR", "alphaR", 1.43, 0.01, 5);
  RooRealVar nL("nL", "nL", 50, 0., 100);
  RooRealVar nR("nR", "nR", 50, 0., 100);
  RooRealVar sigmaL("sigmaL", "sigmaL", 1.66, 0.01, 5);
  RooRealVar sigmaR("sigmaR", "sigmaR", 1.67, 0.01, 5.);

  int np = 6;

  //alphaR.removeMin();
  //alphaR.removeMax();
  //alphaL.removeMin();
  //alphaL.removeMax();
  //nL.removeMin();
  //nL.removeMax();
  //nR.removeMin();
  //nR.removeMax();

//Set alpha or n constant
  // alphaL.setConstant(true);
  // alphaR.setConstant(true);
  // nL.setConstant(true);
  // nR.setConstant(true);
  // mu.setConstant(true);
  gStyle->SetOptStat(0);

  RooCrystalBall sig_model("sig_model","sig_model", x, mu, sigmaL, alphaL, nL, alphaR, nR);

  alphaR.setError(0.01);
  // mu.setError(0.1);
  // sigmaL.setError(0.01);
  // sigmaR.setError(0.01);
  alphaL.setError(0.01);
  nL.setError(1);
  nR.setError(1);

  int qual = -1;
  
  // auto roores = sig_model.chi2FitTo(h_ul, Save(true),PrintLevel(1),Strategy(0));
  // roores->Print("v");


  RooChi2Var *nll1 = new RooChi2Var("nll1","-log(L)",sig_model, h_ul, DataError(RooAbsData::SumW2)) ;
  //  nll->applyWeightSquared(false)
  RooMinimizer *mini = new RooMinimizer(*nll1);  
  mini->setPrintLevel(-1);
  (*mini).setEps(100);
  (*mini).setStrategy(0);
  // (*mini).setOffsetting(true);
  (*mini).minimize("Minuit2","migrad");
  RooFitResult *res = (*mini).save();
  res->Print("v");
  // (*mini).setOffsetting(false);

  auto nll11 = RooChi2Var("nll11","-log(L)",sig_model, h_ul, DataError(RooAbsData::SumW2)) ;
  RooMinimizer *mini1 = new RooMinimizer(nll11);  
  mini1->setPrintLevel(-1);
  (*mini1).setStrategy(0);
  // (*mini1).setOffsetting(true);
  (*mini1).setEps(1);
  (*mini1).minimize("Minuit2","migrad");
  (*mini1).hesse();
  //   qual = calculateConvMatrix(*mini1, nll11);
  res = (*mini1).save();
  //   res->setCovQual(qual);
  res->Print("v");
  //printQual(qual);

  if (nL.getVal() > 90 || nL.getError() > 100) {nL.setVal(50); nL.setError(0); nL.setConstant(true); np--;}
  if (nR.getVal() > 90 || nR.getError() > 100) {nR.setVal(50); nR.setError(0); nR.setConstant(true); np--;}

  if (nL.getVal() < 0.1) {nL.setVal(0); nL.setError(0); nL.setConstant(true); np--;}
  if (nR.getVal() < 0.1) {nR.setVal(0); nR.setError(0); nR.setConstant(true); np--;}
  
  auto nll12 = RooChi2Var("nll12","-log(L)",sig_model, h_ul, DataError(RooAbsData::SumW2)) ;
  RooMinimizer *mini02 = new  RooMinimizer(nll12);  //mini02->setPrintLevel(-1);
  (*mini02).setStrategy(0);
  // (*mini02).setOffsetting(true);
  (*mini02).setEps(0.1);
  (*mini02).minimize("Minuit2","migrad");
  (*mini02).hesse();
//   qual = calculateConvMatrix(*mini02, nll12);
  res = (*mini02).save();
//   res->setCovQual(qual);
  res->Print("v");
  int status1 = res->status();
  int qual1 = qual;
  //printQual(qual);

  double sigmaLVal = sigmaL.getVal();
  sig_model.plotOn(xframe_ul, LineColor(kRed));
  double s1 = sigmaL.getVal();
  double e1 = sigmaL.getError();
  double nl1 = nL.getVal();
  double nl1e = nL.getError();
  double nr1 = nR.getVal();
  double nr1e = nR.getError();
  double al1 = alphaL.getVal();
  double al1e = alphaL.getError();
  double ar1 = alphaR.getVal();
  double ar1e = alphaR.getError();

  //   RooChi2Var chi21("chi21", "chi21", sig_model, h_ul, DataError(RooAbsData::Expected));
  double valChi21 = nll12.getVal();
  
  cout <<endl <<  endl;
  //RooAbsReal* nll1 = sig_model.createNLL(h1);


  
  TCanvas *c2 = new TCanvas("c2", "c2", 600, 600);
  RooDataHist *pdfDataHis1 = sig_model.generateBinned(RooArgSet(x), h_ul.sumEntries(), true);
  TH1D *hpdf1 = (TH1D*)pdfDataHis1->createHistogram("hpdf1",x,Binning(NBIN, LOWX, HIGHX));
  TH1D *hh_ul = (TH1D*)h_ul.createHistogram("hh_ul",x,Binning(NBIN, LOWX, HIGHX));
  TH1D *ratio1 = new TH1D("ratio1", "ratio1", NBIN, LOWX, HIGHX);
  // ratio1 -> Sumw2();
  ratio1 -> Add(hh_ul,hpdf1, 1., -1.);
  for (int i(0); i < NBIN; i++){
    h_ul.get(i);
    double error = h_ul.weightError(RooAbsData::SumW2);
    ratio1 -> SetBinError(i, error);
  }
  c2->cd();
  TPad *pad0 = new TPad("pad0","",0,0.2,1,1);
  TPad *pad00 = new TPad("pad00","",0,0,1,0.2);
  pad0->SetBottomMargin(0.);
  pad00->SetTopMargin(0);
  pad0->Draw();
  pad00->Draw();
  pad0->cd();
  //xframe_re-> GetXaxis()-> SetTitle("");
  xframe_ul->Draw();
  pad00->cd();
  ratio1-> SetTitle("");
  ratio1 -> GetYaxis()-> SetLabelSize(0.1);
  ratio1 -> GetYaxis()-> SetRangeUser(-0.2, 0.2);
  ratio1->SetMarkerStyle(8);
  ratio1->SetMarkerSize(1);
  ratio1 -> GetXaxis()-> SetLabelSize(0.1);
  ratio1 -> GetYaxis()-> SetTitleSize(0.15);
  ratio1 -> GetYaxis()-> SetTitleOffset(0.25);
  ratio1 -> GetYaxis()-> SetTitle("Residual ");
  ratio1 -> GetXaxis()-> SetTitleOffset(-0.3);
  ratio1 -> GetXaxis()-> SetTitleSize(0.15);
  ratio1-> GetXaxis()-> SetTitle("m_{llg}");
  ratio1->Draw("");
  TLine *line0 = new TLine( LOWX, 0, HIGHX, 0);
  line0-> SetLineColor(kRed);
  line0->SetLineStyle(7);
  line0->SetLineWidth(2); 
  line0->Draw("same");
  
  if (OUTPUT) c2->SaveAs(((string)("DSCB/" + name + ".pdf")).c_str());


  sig_model.plotOn(xframe, LineColor(kGreen), Name("DSCB"));
  c4->cd();
  xframe->GetXaxis()->SetRangeUser(LOWX, HIGHX);
  xframe->Draw("same");

  string ul="#splitline{#splitline{DSCB MC fit #sigma = ";
  std::stringstream stream;
  stream << std::fixed << std::setprecision(2) << s1 << "#pm" << e1  <<  "}{" << "nL = " << nl1 << "#pm" << nl1e << " nR = " << nr1 << "#pm" << nr1e << "}}{#splitline{" << 
  "#alphaL = " << al1 << "#pm" << al1e << " #alphaR = " << ar1 << "#pm" << ar1e << "}{Status = " << status1 << " Reduced Chi2 = " << valChi21/(NBIN - np)  << "}}";
  std::string ss1 = stream.str();
  string full_ul = ul + ss1;

  leg->AddEntry("DSCB", full_ul.c_str(), "L");
  if (OUTPUT) fclose (stdout);
  return status1;

}

int fit_hist_CBGauss(RooDataHist &h_ul, RooRealVar &x, string def_name, TCanvas *c4, TLegend *leg, bool output){
  bool OUTPUT = output;
  string name = "CB+Gauss " + def_name;
  string ulname = name + "";
  RooPlot *xframe_ul = x.frame(Title(ulname.c_str()));
  h_ul.plotOn(xframe_ul, LineColor(kBlue), DataError(RooAbsData::Auto));
  RooPlot *xframe = x.frame(Title(name.c_str()));
  string outputfile = "fit_results/" + name + ".txt";
  if (OUTPUT) freopen (outputfile.c_str(),"w",stdout);

  RooRealVar mu("mu", "mu", 124.57, 120., 130.);
  RooRealVar alpha("alpha", "alpha", 7.6328e-01, 0.001,2.3);
  RooRealVar n("n", "n", 20, 0.1, 120);
  RooRealVar sigmaCB("sigmaCB", "sigmaCB", 1.67, 0.1, 5.);
  RooRealVar sigma("sigma", "sigma", 2.5, 0.1, 10.);
  RooRealVar frac("frac","Gauss fraction"       ,0.1,  0., 1.);

  int np = 6;

//Set alpha or n constant
  // alphaL.setConstant(true);
  // alphaR.setConstant(true);
  //nL.setConstant(true);
  //nR.setConstant(true);
  // mu.setConstant(true);
  gStyle->SetOptStat(0);

  RooCrystalBall sig_CB("sig_CB","sig_CB", x, mu, sigmaCB, alpha, n);
  RooGaussian sig_gaus("sig_gaus", "sig_gaus", x, mu, sigma);
  RooAddPdf sig_model("sig_model", "Signal", RooArgList(sig_gaus, sig_CB), frac);

  alpha.setError(0.1);
  mu.setError(0.1);
  sigma.setError(0.01);
  sigmaCB.setError(0.01);
  frac.setError(0.001);
  alpha.setError(0.1);
  n.setError(0.1);

  int qual = -1;

  RooChi2Var *nll1 = new RooChi2Var("nll1","-log(L)",sig_model, h_ul, DataError(RooAbsData::SumW2)) ;
  //  nll->applyWeightSquared(false)
  RooMinimizer *mini = new RooMinimizer(*nll1);
  mini->setPrintLevel(-1);
  (*mini).setEps(100);
  (*mini).setStrategy(0);
  (*mini).minimize("Minuit2","migrad");
  RooFitResult *res = (*mini).save();
  res->Print("v");
  (*mini).setOffsetting(false);

  auto nll11 = RooChi2Var("nll11","-log(L)",sig_model, h_ul, DataError(RooAbsData::SumW2)) ;
  RooMinimizer *mini1 = new RooMinimizer(nll11);
  mini1->setPrintLevel(-1);
  (*mini1).setStrategy(0);
  (*mini1).setEps(1);
  (*mini1).minimize("Minuit2","migrad");
  (*mini1).hesse();
  res = (*mini1).save();
  res->Print("v");
  //printQual(qual);

  if (n.getVal() > 100 || n.getError() > 100) {n.setVal(100); n.setError(0); n.setConstant(true); np--;}
  
  auto nll12 = RooChi2Var("nll12","-log(L)",sig_model, h_ul, DataError(RooAbsData::SumW2)) ;
  RooMinimizer *mini02 = new  RooMinimizer(nll12);
  (*mini02).setStrategy(0);
  (*mini02).setEps(0.1);
  (*mini02).minimize("Minuit2","migrad");
  (*mini02).hesse();
  res = (*mini02).save();
  res->Print("v");
  int status1 = res->status();
  int qual1 = qual;
  //printQual(qual);

  double sigmaVal = sigmaCB.getVal();
  sig_model.plotOn(xframe_ul, LineColor(kRed));
  sig_model.plotOn(xframe, LineColor(kBlue), Name("UL"));
  double Gs1 = sigma.getVal();
  double Ge1 = sigma.getError();
  double CBs1 = sigmaCB.getVal();
  double CBe1 = sigmaCB.getError();
  double Gfrac1 = frac.getVal();
  double Gfrac1e = frac.getError();
  double n1 = n.getVal();
  double n1e = n.getError();
  double a1 = alpha.getVal();
  double a1e = alpha.getError();
  vector<double> combine1 = combine_sigma(Gs1, Ge1, CBs1, CBe1, Gfrac1);

    
  TCanvas *c2 = new TCanvas("c2", "c2", 600, 600);
  RooDataHist *pdfDataHis1 = sig_model.generateBinned(RooArgSet(x), h_ul.sumEntries(), true);
  TH1D *hpdf1 = (TH1D*)pdfDataHis1->createHistogram("hpdf1",x,Binning(NBIN, LOWX, HIGHX));
  TH1D *hh_ul = (TH1D*)h_ul.createHistogram("hh_ul",x,Binning(NBIN, LOWX, HIGHX));
  TH1D *ratio1 = new TH1D("ratio1", "ratio1", NBIN, LOWX, HIGHX);
  // ratio1 -> Sumw2();
  ratio1 -> Add(hh_ul,hpdf1, 1., -1.);
  for (int i(0); i < NBIN; i++){
    h_ul.get(i);
    double error = h_ul.weightError(RooAbsData::SumW2);
    ratio1 -> SetBinError(i, error);
  }
  c2->cd();
  TPad *pad0 = new TPad("pad0","",0,0.2,1,1);
  TPad *pad00 = new TPad("pad00","",0,0,1,0.2);
  pad0->SetBottomMargin(0.);
  pad00->SetTopMargin(0);
  pad0->Draw();
  pad00->Draw();
  pad0->cd();
  //xframe_re-> GetXaxis()-> SetTitle("");
  xframe_ul->Draw();
  pad00->cd();
  ratio1-> SetTitle("");
  ratio1 -> GetYaxis()-> SetLabelSize(0.1);
  ratio1 -> GetYaxis()-> SetRangeUser(-0.2, 0.2);
  ratio1->SetMarkerStyle(8);
  ratio1->SetMarkerSize(1);
  ratio1 -> GetXaxis()-> SetLabelSize(0.1);
  ratio1 -> GetYaxis()-> SetTitleSize(0.15);
  ratio1 -> GetYaxis()-> SetTitleOffset(0.25);
  ratio1 -> GetYaxis()-> SetTitle("Residual ");
  ratio1 -> GetXaxis()-> SetTitleOffset(-0.3);
  ratio1 -> GetXaxis()-> SetTitleSize(0.15);
  ratio1-> GetXaxis()-> SetTitle("m_{llg}");
  ratio1->Draw("");
  TLine *line0 = new TLine( LOWX, 0, HIGHX, 0);
  line0-> SetLineColor(kRed);
  line0->SetLineStyle(7);
  line0->SetLineWidth(2); 
  line0->Draw("same");
  if (OUTPUT) c2->SaveAs(((string)("CBGauss/" + name + ".pdf")).c_str());

  sig_model.plotOn(xframe, LineColor(kRed), Name("CB+Gauss"));
  c4->cd();
  xframe->GetXaxis()->SetRangeUser(LOWX, HIGHX);
  xframe->Draw("same");

  string ul="#splitline{#splitline{CB+Gauss MC fit #sigma = ";
  std::stringstream stream;
  stream << std::fixed << std::setprecision(2) << combine1[0] << "#pm" << combine1[1]  <<  "}{" << "G frac = " << Gfrac1 << "#pm" << Gfrac1e << " n = " << n1 << "#pm" << n1e << "}}{#splitline{" << 
  "#alpha = " << a1 << "#pm" << a1e << "}{Status = " << status1  << "Reduced Chi2 = " << (nll12.getVal())/(NBIN - np) << "}}";
  std::string ss1 = stream.str();
  string full_ul = ul + ss1;


  leg->AddEntry("CB+Gauss", full_ul.c_str(), "L");
  if (OUTPUT) fclose (stdout);
  combine1.clear();
  return status1;

}

int fit_hist_sum4Gaus_order1(RooDataHist &h_ul, RooRealVar &x, string def_name, TCanvas *c4, TLegend *leg, bool output){
  bool OUTPUT = output;
  string name = "Sum4Gauss order1 " + def_name;
  string ulname = name + "";
  RooPlot *xframe_ul = x.frame(Title(ulname.c_str()));
  h_ul.plotOn(xframe_ul, LineColor(kBlue), DataError(RooAbsData::Auto));
  RooPlot *xframe = x.frame(Title(def_name.c_str()));
  string outputfile = "fit_results/" + name + ".txt";
  if (OUTPUT) freopen (outputfile.c_str(),"w",stdout);

  RooRealVar *MH = new RooRealVar("MH", "MH", 125., 123., 127.);
  RooArgList *gaussians_order1 = new RooArgList();
  RooArgList *coeffs_order1 = new RooArgList();
  RooArgSet *listOfPolyVars_ = new RooArgSet();
  // int g = 0;
  for (int g(0); g < 4; g++){
    RooFormulaVar *dMH = new RooFormulaVar("dMH", Form("dMH",g), "@0-125",RooArgList(*MH));
    RooRealVar *dm_p0 = new RooRealVar(Form("dm_g%d_p0",g),Form("dm_g%d_p0",g),0.1,-15.0,15.0);
    RooRealVar *dm_p1 = new RooRealVar(Form("dm_g%d_p1",g),Form("dm_g%d_p1",g),0.05,-0.25,0.25);
    RooRealVar *dm_p2 = new RooRealVar(Form("dm_g%d_p2",g),Form("dm_g%d_p2",g),0.01,-0.01,0.01);
    
    RooPolyVar *dm_order1 = new RooPolyVar(Form("dm_g%d_order1",g),Form("dm_g%d_order1",g),*dMH,RooArgList(*dm_p0,*dm_p1)); //y=a+bx
    RooFormulaVar *mean_order1 = new RooFormulaVar(Form("mean_g%d_order1",g),Form("mean_g%d_order1",g),"((@0+@1))",RooArgList(*MH,*dm_order1));

    RooRealVar *sigma_p0 = new RooRealVar(Form("sigma_g%d_p0",g),Form("sigma_g%d_p0",g),(g+1)*1.0,0.2,8);
    RooRealVar *sigma_p1 = new RooRealVar(Form("sigma_g%d_p1",g),Form("sigma_g%d_p1",g),0.01,-0.05,0.05);
    RooRealVar *sigma_p2 = new RooRealVar(Form("sigma_g%d_p2",g),Form("sigma_g%d_p2",g),0.01,-0.01,0.01);

    RooPolyVar *sigma_order1 = new RooPolyVar(Form("sigma_g%d_order1",g),Form("sigma_g%d_order1",g),*dMH,RooArgList(*sigma_p0,*sigma_p1));
    RooAbsPdf *gaus_order1 = new RooGaussian(Form("gaus_g%d_order1",g),Form("gaus_g%d_order1",g),x,*mean_order1,*sigma_order1);

    gaussians_order1->add(*gaus_order1);
    // listOfPolyVars_->add(*dm_order2);
    listOfPolyVars_->add(*mean_order1);
    listOfPolyVars_->add(*sigma_order1);
    if (g < 3){
      RooRealVar *frac_p0 = new RooRealVar(Form("frac_g%d_p0",g),Form("frac_g%d_p0",g),0.5-0.05*g, 0.01,0.99);
      RooRealVar *frac_p1 = new RooRealVar(Form("frac_g%d_p1",g),Form("frac_g%d_p1",g), 0.0001,-0.01,0.01);
      RooRealVar *frac_p2 = new RooRealVar(Form("frac_g%d_p2",g),Form("frac_g%d_p2",g),0.00001,-0.00001,0.00001);
      RooPolyVar *frac_order1 = new RooPolyVar(Form("frac_g%d_order1",g),Form("frac_g%d_order1",g),*dMH,RooArgList(*frac_p0,*frac_p1));
      RooFormulaVar *frac_constrained_order1 = new RooFormulaVar(Form("frac_g%d_constrained_order1",g),Form("frac_g%d_constrained_order1",g),"(@0>0)*(@0<1)*@0+ (@0>1.0)*0.9999",RooArgList(*frac_order1));
      coeffs_order1->add(*frac_constrained_order1);
    }
    
    
  }
  RooAddPdf sig_model("sig_model","sig_model",*gaussians_order1,*coeffs_order1, true);


  int qual = -1;

  // auto *nll1 = new RooChi2Var("nll1","-log(L)",sig_model, h_ul, DataError(RooAbsData::SumW2)) ;
  // RooMinimizer *mini = new RooMinimizer(*nll1);  mini->setPrintLevel(-1);
  // (*mini).setEps(100);
  // (*mini).setStrategy(0);
  // (*mini).minimize("Minuit2","migrad");
  // RooFitResult *res = (*mini).save();
  // res->Print("v");
  // (*mini).setOffsetting(false);

  // auto nll11 = RooChi2Var("nll11","-log(L)",sig_model, h_ul, DataError(RooAbsData::SumW2)) ;
  // RooMinimizer *mini1 = new RooMinimizer(nll11);  mini1->setPrintLevel(-1);
  // (*mini1).setStrategy(0);
  // (*mini1).setEps(0.1);
  // (*mini1).minimize("Minuit2","migrad");
  // (*mini1).hesse();
  // res = (*mini1).save();
  // res->Print("v");


  // auto nll12 = RooChi2Var("nll12","-log(L)",sig_model, h_ul, DataError(RooAbsData::SumW2)) ;
  // RooMinimizer *mini02 = new  RooMinimizer(nll12);  //mini02->setPrintLevel(-1);
  // (*mini02).setStrategy(0);
  // (*mini02).setEps(0.01);
  // (*mini02).minimize("Minuit2","migrad");
  // (*mini02).hesse();
  // res = (*mini02).save();
  // res->Print("v");
  int status1 = 0; //res->status();
//   int qual1 = qual;
  //printQual(qual);


  // RooChi2Var chi21("chi21", "chi21", sig_model, h_ul, DataError(RooAbsData::Expected));
  double valChi21 = 0;//nll12.getVal();
  
  cout <<endl <<  endl;
  //RooAbsReal* nll1 = sig_model.createNLL(h1);

	// sig_model.plotOn(xframe_ul, LineColor(kRed));
  
  // TCanvas *c2 = new TCanvas("c2", "c2", 600, 600);
  // RooDataHist *pdfDataHis1 = sig_model.generateBinned(RooArgSet(x), h_ul.sumEntries(), true);
  // TH1D *hpdf1 = (TH1D*)pdfDataHis1->createHistogram("hpdf1",x,Binning(NBIN, LOWX, HIGHX));
  // TH1D *hh_ul = (TH1D*)h_ul.createHistogram("hh_ul",x,Binning(NBIN, LOWX, HIGHX));
  // TH1D *ratio1 = new TH1D("ratio1", "ratio1", NBIN, LOWX, HIGHX);
  // // ratio1 -> Sumw2();
  // ratio1 -> Add(hh_ul,hpdf1, 1., -1.);
  // for (int i(0); i < NBIN; i++){
  //   h_ul.get(i);
  //   double error = h_ul.weightError(RooAbsData::SumW2);
  //   ratio1 -> SetBinError(i, error);
  // }
  // c2->cd();
  // TPad *pad0 = new TPad("pad0","",0,0.2,1,1);
  // TPad *pad00 = new TPad("pad00","",0,0,1,0.2);
  // pad0->SetBottomMargin(0.);
  // pad00->SetTopMargin(0);
  // pad0->Draw();
  // pad00->Draw();
  // pad0->cd();
  // //xframe_re-> GetXaxis()-> SetTitle("");
  // xframe_ul->Draw();
  // pad00->cd();
  // ratio1-> SetTitle("");
  // ratio1 -> GetYaxis()-> SetLabelSize(0.1);
  // ratio1 -> GetYaxis()-> SetRangeUser(-0.2, 0.2);
  // ratio1->SetMarkerStyle(8);
  // ratio1->SetMarkerSize(1);
  // ratio1 -> GetXaxis()-> SetLabelSize(0.1);
  // ratio1 -> GetYaxis()-> SetTitleSize(0.15);
  // ratio1 -> GetYaxis()-> SetTitleOffset(0.25);
  // ratio1 -> GetYaxis()-> SetTitle("Residual ");
  // ratio1 -> GetXaxis()-> SetTitleOffset(-0.3);
  // ratio1 -> GetXaxis()-> SetTitleSize(0.15);
  // ratio1-> GetXaxis()-> SetTitle("m_{llg}");
  // ratio1->Draw("");
  // TLine *line0 = new TLine( LOWX, 0, HIGHX, 0);
  // line0-> SetLineColor(kRed);
  // line0->SetLineStyle(7);
  // line0->SetLineWidth(2); 
  // line0->Draw("same");
  
  // if (OUTPUT) c2->SaveAs(((string)("SumGauss/" + name + ".pdf")).c_str());

  sig_model.plotOn(xframe, LineColor(kBlue), Name("SumG"));
  c4->cd();
  xframe->GetXaxis()->SetRangeUser(LOWX, HIGHX);
  xframe->Draw("same");


  // TLatex *txt = new TLatx(0, 0.01, )
  // leg1->SetTextSize(0.02);

  string ul="#splitline{SumGaus MC fit";
  std::stringstream stream;
  stream << std::fixed << std::setprecision(2) << "}{Status = " <<  status1 << " Reduced Chi2 = " << valChi21/(NBIN - 23) <<"}";
  std::string ss1 = stream.str();
  string full_ul = ul + ss1;

  leg->AddEntry("SumG", full_ul.c_str(), "L");
  
  // (*listOfPolyVars_)["sigma_g0_order1"].Print();
  // (*listOfPolyVars_)["mean_g0_order1"].Print();
  // (*listOfPolyVars_)["sigma_g1_order1"].Print();
  // (*listOfPolyVars_)["mean_g1_order1"].Print();
  // (*listOfPolyVars_)["sigma_g2_order1"].Print();
  // (*listOfPolyVars_)["mean_g2_order1"].Print();
  // (*listOfPolyVars_)["sigma_g3_order1"].Print();
  // (*listOfPolyVars_)["mean_g3_order1"].Print();

  if (OUTPUT) fclose (stdout);
  return status1;

}


int fit_hist_sum4Gaus_order0(RooDataHist &h_ul, RooRealVar &x, string def_name, TCanvas *c4, TLegend *leg, bool output){
  bool OUTPUT = output;
  string name = "Sum4Gauss order0 " + def_name;
  string ulname = name + "";
  RooPlot *xframe_ul = x.frame(Title(ulname.c_str()));
  h_ul.plotOn(xframe_ul, LineColor(kBlue), DataError(RooAbsData::Auto));
  RooPlot *xframe = x.frame(Title(def_name.c_str()));
  string outputfile = "fit_results/" + name + ".txt";
  if (OUTPUT) freopen (outputfile.c_str(),"w",stdout);
  
  RooRealVar *MH = new RooRealVar("MH", "MH", 125.);
  RooArgList *gaussians_order0 = new RooArgList();
  RooArgList *coeffs_order0 = new RooArgList();
  RooArgSet *listOfPolyVars_ = new RooArgSet();
  // int g = 0;
  for (int g(0); g < 4; g++){
    RooRealVar *dm_p0 = new RooRealVar(Form("dm_g%d_p0",g),Form("dm_g%d_p0",g),0.1,-15.0,15.0);
    RooFormulaVar *mean_order0 = new RooFormulaVar(Form("mean_g%d_order0",g),Form("mean_g%d_order0",g),"((@0+@1))",RooArgList(*MH,*dm_p0));
    RooRealVar *sigma_order0 = new RooRealVar(Form("sigma_g%d_order0",g),Form("sigma_g%d_order0",g),(g+1)*1.0,0.2,8);
    RooAbsPdf *gaus_order0 = new RooGaussian(Form("gaus_g%d_order0",g),Form("gaus_g%d_order0",g),x,*mean_order0,*sigma_order0);
    
    gaussians_order0->add(*gaus_order0);
    // listOfPolyVars_->add(*dm_order2);
    listOfPolyVars_->add(*mean_order0);
    listOfPolyVars_->add(*sigma_order0);
    if (g < 3){
      RooRealVar *frac_p0 = new RooRealVar(Form("frac_g%d_p0",g),Form("frac_g%d_p0",g),0.5-0.05*g, 0.01,0.99);
      coeffs_order0->add(*frac_p0);
    }
    
    
  }
  RooAddPdf sig_model("sig_model","sig_model",*gaussians_order0,*coeffs_order0, true);
  
  
  int qual = 3;

  auto *nll1 = new RooChi2Var("nll1","-log(L)",sig_model, h_ul, DataError(RooAbsData::SumW2)) ;
  RooMinimizer *mini = new RooMinimizer(*nll1);  mini->setPrintLevel(-1);
  (*mini).setEps(100);
  (*mini).setStrategy(0);
  (*mini).minimize("Minuit2","migrad");
  RooFitResult *res = (*mini).save();
  res->Print("v");

  auto nll11 = RooChi2Var("nll11","-log(L)",sig_model, h_ul, DataError(RooAbsData::SumW2)) ;
  RooMinimizer *mini1 = new RooMinimizer(nll11);  mini1->setPrintLevel(-1);
  (*mini1).setStrategy(0);
  (*mini1).setEps(1);
  (*mini1).minimize("Minuit2","migrad");
  (*mini1).hesse();


  auto nll12 = RooChi2Var("nll12","-log(L)",sig_model, h_ul, DataError(RooAbsData::SumW2)) ;
  RooMinimizer *mini02 = new  RooMinimizer(nll12);  //mini02->setPrintLevel(-1);
  (*mini02).setStrategy(0);
  (*mini02).setEps(0.01);
  (*mini02).minimize("Minuit2","migrad");
  (*mini02).hesse();
  res = (*mini02).save();
  res->Print("v");
  int status1 = res->status();
  int qual1 = qual;
  //printQual(qual);


  // RooChi2Var chi21("chi21", "chi21", sig_model, h_ul, DataError(RooAbsData::Expected));
  double valChi21 =nll12.getVal();
  
  cout <<endl <<  endl;

  sig_model.plotOn(xframe_ul, LineColor(kRed));
  
  TCanvas *c2 = new TCanvas("c2", "c2", 600, 600);
  RooDataHist *pdfDataHis1 = sig_model.generateBinned(RooArgSet(x), h_ul.sumEntries(), true);
  TH1D *hpdf1 = (TH1D*)pdfDataHis1->createHistogram("hpdf1",x,Binning(NBIN, LOWX, HIGHX));
  TH1D *hh_ul = (TH1D*)h_ul.createHistogram("hh_ul",x,Binning(NBIN, LOWX, HIGHX));
  TH1D *ratio1 = new TH1D("ratio1", "ratio1", NBIN, LOWX, HIGHX);
  // ratio1 -> Sumw2();
  ratio1 -> Add(hh_ul,hpdf1, 1., -1.);
  for (int i(0); i < NBIN; i++){
    h_ul.get(i);
    double error = h_ul.weightError(RooAbsData::SumW2);
    ratio1 -> SetBinError(i, error);
  }
  c2->cd();
  TPad *pad0 = new TPad("pad0","",0,0.2,1,1);
  TPad *pad00 = new TPad("pad00","",0,0,1,0.2);
  pad0->SetBottomMargin(0.);
  pad00->SetTopMargin(0);
  pad0->Draw();
  pad00->Draw();
  pad0->cd();
  //xframe_re-> GetXaxis()-> SetTitle("");
  xframe_ul->Draw();
  pad00->cd();
  ratio1-> SetTitle("");
  ratio1 -> GetYaxis()-> SetLabelSize(0.1);
  ratio1 -> GetYaxis()-> SetRangeUser(-0.2, 0.2);
  ratio1->SetMarkerStyle(8);
  ratio1->SetMarkerSize(1);
  ratio1 -> GetXaxis()-> SetLabelSize(0.1);
  ratio1 -> GetYaxis()-> SetTitleSize(0.15);
  ratio1 -> GetYaxis()-> SetTitleOffset(0.25);
  ratio1 -> GetYaxis()-> SetTitle("Residual ");
  ratio1 -> GetXaxis()-> SetTitleOffset(-0.3);
  ratio1 -> GetXaxis()-> SetTitleSize(0.15);
  ratio1-> GetXaxis()-> SetTitle("m_{llg}");
  ratio1->Draw("");
  TLine *line0 = new TLine( LOWX, 0, HIGHX, 0);
  line0-> SetLineColor(kRed);
  line0->SetLineStyle(7);
  line0->SetLineWidth(2); 
  line0->Draw("same");
  
  if (OUTPUT) c2->SaveAs(((string)("SumGauss/" + name + ".pdf")).c_str());

	sig_model.plotOn(xframe, LineColor(kBlue), Name("SumG"));
  c4->cd();
  xframe->GetXaxis()->SetRangeUser(LOWX, HIGHX);
  xframe->Draw("same");


  // TLatex *txt = new TLatx(0, 0.01, )
  // leg1->SetTextSize(0.02);

  string ul="#splitline{SumGaus MC fit";
  std::stringstream stream;
  stream << std::fixed << std::setprecision(2) << "}{Status = " <<  status1 << " Reduced Chi2 = " << valChi21/(NBIN - 11) <<"}";
  std::string ss1 = stream.str();
  string full_ul = ul + ss1;

  leg->AddEntry("SumG", full_ul.c_str(), "L");

  // (*listOfPolyVars_)["sigma_g0_order1"].Print();
  // (*listOfPolyVars_)["mean_g0_order1"].Print();
  // (*listOfPolyVars_)["sigma_g1_order1"].Print();
  // (*listOfPolyVars_)["mean_g1_order1"].Print();
  // (*listOfPolyVars_)["sigma_g2_order1"].Print();
  // (*listOfPolyVars_)["mean_g2_order1"].Print();
  // (*listOfPolyVars_)["sigma_g3_order1"].Print();
  // (*listOfPolyVars_)["mean_g3_order1"].Print();

  if (OUTPUT) fclose (stdout);
  return status1;

}

void sig_fit_chi2(string yr, string cat, string lepton){
  if (cat = "Cat 1" && lepton == "Muon"){
    HIGHX = 130;
    NBIN = 60;
  }
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL) ;
  RooRealVar x("x", "m_llg", LOWX, HIGHX);
  RooRealVar y("y", "photon pT", 15., 1000.);
  RooRealVar w("w", "weight", -40., 40.);
  RooRealVar bdt("bdt", "bdt", -1, 1);
  RooRealVar year("year", "year", 2015, 2019);
  RooRealVar lep("lep", "lep", 0, 1); //0 = electron, 1 = muon
  RooRealVar ph_eta("ph_eta", "ph_eta", -3, 3);
  RooRealVar nlep("nlep", "nlep", 0, 10);
  RooRealVar njet("njet", "njet", 0, 10);
  gStyle->SetOptStat(0);

  RooDataSet* ULsample = RooDataSet::read("../../../../CMSSW_11_3_4/src/HToZg_combine/rui_datasamples/Signal_untagged4_ratio_extended.dat,../../../../CMSSW_11_3_4/src/HToZg_combine/rui_datasamples/Signal_untagged3_ratio_extended.dat,../../../../CMSSW_11_3_4/src/HToZg_combine/rui_datasamples/Signal_untagged2_ratio_extended.dat,../../../../CMSSW_11_3_4/src/HToZg_combine/rui_datasamples/Signal_untagged1_ratio_extended.dat", RooArgList(x, y, bdt, w, year, lep, ph_eta, nlep, njet));


  RooDataSet f_data_sig("f_data_sig", "f_data_sig", ULsample, RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), yr.c_str(),"w");



  RooDataSet * d_el_u1 = (RooDataSet*)f_data_sig.reduce(RooArgSet(x,ph_eta),"lep < 1 && nlep <= 2 && njet < 2 && bdt > 0.66");
  RooDataSet *  d_el_u2 = (RooDataSet*)f_data_sig.reduce(RooArgSet(x,ph_eta),"lep < 1 && nlep <= 2 && njet < 2 && bdt < 0.66 && bdt > 0.38");
  RooDataSet * d_el_u3 = (RooDataSet*)f_data_sig.reduce(RooArgSet(x,ph_eta),"lep < 1 && nlep <= 2 && njet < 2 && bdt < 0.38 && bdt > 0.02");
  RooDataSet * d_el_u4 = (RooDataSet*)f_data_sig.reduce(RooArgSet(x,ph_eta),"lep < 1 && nlep <= 2 && njet < 2 && bdt < 0.02 && bdt > -0.58");

  RooDataSet * d_mu_u1 = (RooDataSet*)f_data_sig.reduce(RooArgSet(x,ph_eta),"lep > 0 && nlep <= 2 && njet < 2 && bdt > 0.66");
  RooDataSet *  d_mu_u2 = (RooDataSet*)f_data_sig.reduce(RooArgSet(x,ph_eta),"lep > 0 && nlep <= 2 && njet < 2 && bdt < 0.66 && bdt > 0.38");
  RooDataSet * d_mu_u3 = (RooDataSet*)f_data_sig.reduce(RooArgSet(x,ph_eta),"lep > 0 && nlep <= 2 && njet < 2 && bdt < 0.38 && bdt > 0.02");
  RooDataSet * d_mu_u4 = (RooDataSet*)f_data_sig.reduce(RooArgSet(x,ph_eta),"lep > 0 && nlep <= 2 && njet < 2 && bdt < 0.02 && bdt > -0.58");

  RooDataSet *d_el_B_u1 = (RooDataSet*)d_el_u1->reduce(RooArgSet(x),"ph_eta > -1.4442 && ph_eta < 1.4442");
  RooDataSet *  d_el_E_u1 = (RooDataSet*)d_el_u1->reduce(RooArgSet(x),"(ph_eta > -2.5 && ph_eta < -1.566) || (ph_eta > 1.566 && ph_eta < 2.5)");
  RooDataSet *  d_el_B_u2 = (RooDataSet*)d_el_u2->reduce(RooArgSet(x),"ph_eta > -1.4442 && ph_eta < 1.4442");
  RooDataSet *  d_el_E_u2 = (RooDataSet*)d_el_u2->reduce(RooArgSet(x),"(ph_eta > -2.5 && ph_eta < -1.566) || (ph_eta > 1.566 && ph_eta < 2.5)");
  RooDataSet *  d_el_B_u3 = (RooDataSet*)d_el_u3->reduce(RooArgSet(x),"ph_eta > -1.4442 && ph_eta < 1.4442");
  RooDataSet *  d_el_E_u3 = (RooDataSet*)d_el_u3->reduce(RooArgSet(x),"(ph_eta > -2.5 && ph_eta < -1.566) || (ph_eta > 1.566 && ph_eta < 2.5)");
  RooDataSet *  d_el_B_u4 = (RooDataSet*)d_el_u4->reduce(RooArgSet(x),"ph_eta > -1.4442 && ph_eta < 1.4442");
  RooDataSet *  d_el_E_u4 = (RooDataSet*)d_el_u4->reduce(RooArgSet(x),"(ph_eta > -2.5 && ph_eta < -1.566) || (ph_eta > 1.566 && ph_eta < 2.5)");

  TCanvas *c4 = new TCanvas("c4", "c4", 1200, 1000);
  TLegend *leg = new TLegend(0.1,0.9,0.45,0.5);
  leg->SetTextSize(0.015);

  RooDataHist *h_fit;
  x.setBins(NBIN);
  if (cat == "Cat 1" && lepton == "Electron") {
    RooDataHist h_fit_("h_fit", "h_fit", RooArgSet(x), *d_el_u1);
    h_fit = (RooDataHist *)h_fit_.Clone();
  }
  else if (cat == "Cat 1" && lepton == "Muon") {
    RooDataHist h_fit_("h_fit", "h_fit", RooArgSet(x), *d_mu_u1);
    h_fit = (RooDataHist *)h_fit_.Clone();
  }
  else if (cat == "Cat 2" && lepton == "Electron") {
    RooDataHist h_fit_("h_fit", "h_fit", RooArgSet(x), *d_el_u2);
    h_fit = (RooDataHist *)h_fit_.Clone();
  }
  else if (cat == "Cat 2" && lepton == "Muon") {
    RooDataHist h_fit_("h_fit", "h_fit", RooArgSet(x), *d_mu_u2);
    h_fit = (RooDataHist *)h_fit_.Clone();
  }
  else if (cat == "Cat 3" && lepton == "Electron") {
    RooDataHist h_fit_("h_fit", "h_fit", RooArgSet(x), *d_el_u3);
    h_fit = (RooDataHist *)h_fit_.Clone();
  }
  else if (cat == "Cat 3" && lepton == "Muon") {
    RooDataHist h_fit_("h_fit", "h_fit", RooArgSet(x), *d_mu_u3);
    h_fit = (RooDataHist *)h_fit_.Clone();
  }
  else if (cat == "Cat 4" && lepton == "Electron"){
    RooDataHist h_fit_("h_fit", "h_fit", RooArgSet(x), *d_el_u4);
    h_fit = (RooDataHist *)h_fit_.Clone();
  }
  else if (cat == "Cat 4" && lepton == "Muon"){
    RooDataHist h_fit_("h_fit", "h_fit", RooArgSet(x), *d_mu_u4);
    h_fit = (RooDataHist *)h_fit_.Clone();
  }

  string some_title = yr + " " + lepton +" GGF " + cat;
  
  int status1 = fit_hist_sum4Gaus_order0(*h_fit, x, some_title.c_str(), c4, leg, true);
  int status2 = fit_hist_DSCB(*h_fit, x, some_title.c_str(), c4, leg, true);
  int status3 = fit_hist_CBGauss(*h_fit, x, some_title.c_str(), c4, leg, true); 
  RooPlot *xframe_hist = x.frame();
  h_fit->plotOn(xframe_hist);
  
  c4->cd();
  xframe_hist->Draw("same");
  leg->Draw("same");
  c4->Draw();
  c4->SaveAs(("plots/" + some_title + ".pdf").c_str());
  
  if (status1==0 && status2==0 && status3==0) std::cout <<"Minimization all good!" << std::endl;
  else std::cout << "Minimization fails!" << std::endl;

}
    
