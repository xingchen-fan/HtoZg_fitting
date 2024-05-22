#ifndef EXMODGAUS
#define EXMODGAUS

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

class EXModGaus : public RooAbsPdf {
public:
  EXModGaus() {} ;
  EXModGaus(const char *name, const char *title,
              RooAbsReal& _m,
              RooAbsReal& _mu,
              RooAbsReal& _sig,
              RooAbsReal& _lambda,
              double xlow_,
              double xhigh_);
  EXModGaus(const EXModGaus& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new EXModGaus(*this,newname); }
  inline virtual ~EXModGaus() { }
  virtual bool selfNormalized() const {
    return true;
  }
  void updateNorm(std::vector<Double_t> fitpars) const;

protected:

  RooRealProxy m ;
  RooRealProxy mu ;
  RooRealProxy sig ;
  RooRealProxy lambda ;

  double xlow;
  double xhigh;
  Double_t evaluate() const ;

private:
  mutable std::vector<Double_t> current;
  mutable Double_t norm;
  //Double_t EXModGaus::mfunc(Double_t* x, Double_t* pr);

  ClassDef(EXModGaus,1) // Your description goes here...
};

#endif


