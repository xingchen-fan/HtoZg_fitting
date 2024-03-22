#ifndef MODGAUS
#define MODGAUS

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

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

#endif


