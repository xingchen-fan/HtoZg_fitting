#ifndef ASYMGENGAUSSIAN
#define ASYMGENGAUSSIAN

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

class AsymGenGaussian : public RooAbsPdf {
public:
  AsymGenGaussian() {} ;
  AsymGenGaussian(const char *name, const char *title,
              RooAbsReal& _m,
              RooAbsReal& _kappa,
              RooAbsReal& _alpha,
              RooAbsReal& _xsi,
//              RooAbsReal& _sigma,
              double xlow_,
              double xhigh_);
  AsymGenGaussian(const AsymGenGaussian& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new AsymGenGaussian(*this,newname); }
  inline virtual ~AsymGenGaussian() { }
  virtual bool selfNormalized() const {
    return true;
  }
  void updateNorm(std::vector<Double_t> fitpars) const;

protected:

  RooRealProxy m ;
  RooRealProxy kappa ;
  RooRealProxy alpha ;
//  RooRealProxy sigma ;
  RooRealProxy xsi ;

  double xlow;
  double xhigh;
  Double_t evaluate() const ;

private:
  mutable std::vector<Double_t> current;
  mutable Double_t norm;
  //Double_t mfunc(Double_t* x, Double_t* pr);

  ClassDef(AsymGenGaussian,1) // Your description goes here...
};

#endif


