/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include "Riostream.h" 

#include "EXModGaus.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h> 
#include "TMath.h" 
#include "TF1.h"

ClassImp(EXModGaus); 

 EXModGaus::EXModGaus(const char *name, const char *title, 
                        RooAbsReal& _m,
                        RooAbsReal& _mu,
                        RooAbsReal& _sig,
                        RooAbsReal& _lambda,
		                    double xlow_,
                        double xhigh_) :
   RooAbsPdf(name,title), 
   m("m", "m", this, _m),
   mu("mu", "mu", this, _mu),
   sig("sig", "sig", this, _sig),
   lambda("lambda", "lambda", this, _lambda),
   xlow(xlow_),
   xhigh(xhigh_)
 { 
 } 


 EXModGaus::EXModGaus(const EXModGaus& other, const char* name) :  
   RooAbsPdf(other,name), 
   m("m", this,other.m),
   mu("mu", this,other.mu),
   sig("sig", this,other.sig),
   lambda("lambda", this,other.lambda),
   xlow(other.xlow),
   xhigh(other.xhigh)
 { 
 } 


 // Original
 
// Double_t EXModGaus::mfunc(Double_t* x, Double_t* pr){
/*    Double_t width, power;
    // vh - vl >> vh
    power = pr[1] + pr[2]*(x[0] - pr[6])/(pr[7]-pr[6]);
    if(x[0] <= pr[0]) width = pr[4] + (pr[3] - pr[4])*(x[0] - 105)/(pr[0]-pr[6]);
    else width = pr[3] + (pr[5] - pr[3])*(x[0] - pr[0])/(pr[7] - pr[0]);
    return exp(-1*pow(fabs((x[0]-pr[0])/width), power));
*/

Double_t exmg_mfunc(Double_t* x, Double_t* pr){
  Double_t sqrt2 = 1.414213562;
  Double_t erfc_comp  = TMath::Erfc( (pr[0] + pr[2]*pr[1]*pr[1] - x[0])/(sqrt2*pr[1]) );
  Double_t gauss_comp = TMath::Exp( pr[2]/2.0*(2*pr[0] + pr[2]*pr[1]*pr[1] - 2.0*x[0])); 
  return pr[2]/2.0*gauss_comp*erfc_comp;

}

void EXModGaus::updateNorm(std::vector<Double_t> fitpars) const{
    if (fitpars != current) {
        int nbins = 10000;
        Double_t low = xlow;
        Double_t high = xhigh;
        Double_t _norm = 0;
        TF1 intfunc("func", exmg_mfunc, 0., 200., 8);
        for(int i(0); i<nbins; i++) {
            intfunc.SetParameters(mu, sig, lambda, low, high);
            _norm+= ((high - low)/nbins) * intfunc.Eval(low + (i + 0.5)*(high - low)/nbins);

        }
    norm = _norm;
    //cout << "norm = " << norm << endl;
    }
}

/*
  // vh - vl >> vh
  power = vl + vr*(m - xlow)/(xhigh - xlow);
  if(m <= m0) width = sl + (s0 - sl)*(m - xlow)/(m0-xlow);
  else width = s0 + (sh - s0)*(m - m0)/(xhigh - m0);
  // TF1 modg("modg", mfunc, 105, 170, 6);
  // modg.SetParameters(m0, vl, vr, s0, sl, sh);
  // Double_t selfnorm = modg.Integral(105, 170); 
*/

 Double_t EXModGaus::evaluate() const 
 { 
   // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 
  Double_t width, power;
  std::vector<Double_t> vect;
  vect.push_back(mu); vect.push_back(sig); vect.push_back(lambda);
  updateNorm(vect);
  vect.clear();
  current.clear();
  current.insert(current.end(), {mu, sig, lambda});  
 
  //Defines the return here
  Double_t sqrt2 = 1.414213562;
  Double_t erfc_comp  = TMath::Erfc( (mu + lambda*sig*sig - m)/(sqrt2*sig) );
  Double_t gauss_comp = TMath::Exp( lambda/2.0*(2*mu + lambda*sig*sig - 2*m)); 
  return lambda/2.0*gauss_comp*erfc_comp/norm;
 }
//  return exp(-1*pow(fabs((m-m0)/width), power))/norm;

//  Double_t mfunc(Double_t* x, Double_t* pr){
//     Double_t width, power;
//     // vh - vl >> vh
//     power = pr[1] + pr[2]*(x[0] - 100)/(180-100);
//     if(x[0] <= pr[0]) width = pr[3] + pr[4]*(x[0] - pr[0]);
//     else width = pr[3] + pr[5]*(x[0] - pr[0]);
//     return exp(-1*pow(fabs((x[0]-pr[0])/width), power));

// }
//  Double_t EXModGaus::evaluate() const 
//  { 
//   // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 
//   Double_t width, power;
//   // vh - vl >> vh
//   power = vl + vr*(m - 100)/(180-100);
//   if(m <= m0) width = s0 + sl*(m - m0);
//   else width = s0 + sh*(m - m0);
//   TF1 modg("modg", mfunc, 100, 180, 6);
//   modg.SetParameters(m0, vl, vr, s0, sl, sh);
//   Double_t selfnorm = modg.Integral(100, 180);

//   return exp(-1*pow(fabs((m-m0)/width), power))/selfnorm;
//  } 

//  // Type 0. Please make s0 a constant because it is not used.
//  Double_t mfunc(Double_t* x, Double_t* pr){
//     Double_t width, power;
//     // vh - vl >> vh
//     power = pr[1] + pr[2]*(x[0] - 100)/(180-100);
//     width = pr[3] + pr[4]*(x[0] - 100)/80;
//     return exp(-1*pow(fabs((x[0]-pr[0])/width), power));

// }
//  Double_t EXModGaus::evaluate() const 
//  { 
//    // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 
//   Double_t width, power;
//   // vh - vl >> vh
//   power = vl + vr*(m - 100)/(180-100);
//   width = sl + sh*(m - 100)/(180-100);
//   TF1 modg("modg", mfunc, 100, 180, 5);
//   modg.SetParameters(m0, vl, vr, sl, sh);
//   Double_t selfnorm = modg.Integral(110, 180);
//   return exp(-1*pow(fabs((m-m0)/width), power))/selfnorm;
//  } 

 // Type 1-1
//  Double_t mfunc(Double_t* x, Double_t* pr){
//     Double_t width, power;
//     // vh - vl >> vh
//     power = pr[1] + pr[2]*(x[0] - 100)/(180-100);
//     width = pr[4] + pr[3]*(x[0] - 100)/(180 - 100) + pr[5] * pow((x[0] - 100)/(180 - 100),2);
//     return exp(-1*pow(fabs((x[0]-pr[0])/width), power));

// }
 
//  Double_t EXModGaus::evaluate() const 
//  { 
//   // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 
//  Double_t width, power;
//  // vh - vl >> vh
//  power = vl + vr*(m - 100)/(180-100);
//  width = sl + s0 * (m - 100)/(180-100) + sh * pow((m-100)/(180-100),2);
//  TF1 modg("modg", mfunc, 100, 180, 6);
//   modg.SetParameters(m0, vl, vr, s0, sl, sh);
//   Double_t selfnorm = modg.Integral(100, 180);
//  return exp(-1*pow(fabs((m-m0)/width), power))/selfnorm;
//  } 

 //// Type 1-2
 //Double_t EXModGaus::evaluate() const 
 //{ 
 //  // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 
 // Double_t width, power;
 // // vh - vl >> vh
 // width = vl + vr*(m - 100)/(180-100);
 // power = sl + s0 * (m - 100)/(180-100) + sh * pow((m-100)/(180-100),2);
 // return exp(-1*pow(fabs((m-m0)/width), power));
 //} 