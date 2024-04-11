/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: HZGRooPdfs.h,v 1.1 2013/02/18 11:52:16 gpetrucc Exp $
 * Authors:                                                                  *
 *   Kyle Cranmer (L. Gray for step function piece)
 *                                                                           *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_STEPBERNSTEIN
#define ROO_STEPBERNSTEIN

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"

class RooRealVar;
class RooArgList ;

class RooStepBernstein : public RooAbsPdf {
public:

  RooStepBernstein() ;
  RooStepBernstein(const char *name, const char *title,
		   RooAbsReal& _x, RooAbsReal& _stepThresh,
		   const RooArgList& _coefList) ;

  RooStepBernstein(const RooStepBernstein& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new RooStepBernstein(*this, newname); }
  inline virtual ~RooStepBernstein() { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

private:

  RooRealProxy _x;
  RooRealProxy _stepThresh;
  RooListProxy _coefList ;

  Double_t evaluate() const;

  ClassDef(RooStepBernstein,1) // Bernstein polynomial PDF with step function
};

class RooGaussStepBernstein : public RooAbsPdf {
public:

  RooGaussStepBernstein() ;
  RooGaussStepBernstein(const char *name, const char *title,
			RooAbsReal& _x, RooAbsReal& _mean,
			RooAbsReal& _sigma, RooAbsReal& _stepVal,
			const RooArgList& _coefList) ;

  RooGaussStepBernstein(const RooGaussStepBernstein& other,
const char* name = 0);
  virtual TObject* clone(const char* newname) const
  { return new RooGaussStepBernstein(*this, newname); }
  inline virtual ~RooGaussStepBernstein() { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

private:

  RooRealProxy _x;
  RooRealProxy _mean,_sigma,_stepVal;
  RooListProxy _coefList ;

  Double_t evaluate() const;

  ClassDef(RooGaussStepBernstein,1) // Bernstein polynomial PDF with step function convoluted with gaussian
};

//  Modified Gauss
class ModGaus : public RooAbsPdf {
public:
  ModGaus() {} ;
  ModGaus(const char *name, const char *title,
              RooAbsReal& _m,
              RooAbsReal& _m0,
              RooAbsReal& _vl,
              RooAbsReal& _vr,
              RooAbsReal& _s0,
              RooAbsReal& _sl,
              RooAbsReal& _sh,
              double xlow_,
              double xhigh_);
  ModGaus(const ModGaus& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new ModGaus(*this,newname); }
  inline virtual ~ModGaus() { }
  virtual bool selfNormalized() const {
    return true;
  }
  void updateNorm(std::vector<Double_t> fitpars) const;

protected:

  RooRealProxy m ;
  RooRealProxy m0 ;
  RooRealProxy vl ;
  RooRealProxy vr ;
  RooRealProxy s0 ;
  RooRealProxy sl ;
  RooRealProxy sh ;
  double xlow;
  double xhigh;
  Double_t evaluate() const ;

private:
  mutable std::vector<Double_t> current;
  mutable Double_t norm;

  ClassDef(ModGaus,1) // Your description goes here...
};

class dummy{
    public:
        dummy(float k);
};

class NormPow : public RooAbsPdf {
  public:
    NormPow() {};
    NormPow(const char *name, const char *title,
      RooAbsReal& _m,
      RooAbsReal& _t,
      RooAbsReal& _p,
      double xhigh_);
    NormPow(const NormPow& other, const char* name=0) ;
    virtual TObject* clone(const char* newname) const { return new NormPow(*this,newname); }
    inline virtual ~NormPow() { }
    virtual bool selfNormalized() const {
      return true;
    }
  protected:
    RooRealProxy m ;
    RooRealProxy t ;
    RooRealProxy p ;
    double xhigh;
    Double_t evaluate() const ;

    ClassDef(NormPow,1)
};
#endif
