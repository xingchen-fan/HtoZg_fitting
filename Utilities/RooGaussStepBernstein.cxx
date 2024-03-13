// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME RooGaussStepBernstein
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "RooGaussStepBernstein.h"
#include "RooGaussStepBernstein.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_RooGaussStepBernstein(void *p = 0);
   static void *newArray_RooGaussStepBernstein(Long_t size, void *p);
   static void delete_RooGaussStepBernstein(void *p);
   static void deleteArray_RooGaussStepBernstein(void *p);
   static void destruct_RooGaussStepBernstein(void *p);
   static void streamer_RooGaussStepBernstein(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RooGaussStepBernstein*)
   {
      ::RooGaussStepBernstein *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RooGaussStepBernstein >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RooGaussStepBernstein", ::RooGaussStepBernstein::Class_Version(), "RooGaussStepBernstein.h", 23,
                  typeid(::RooGaussStepBernstein), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RooGaussStepBernstein::Dictionary, isa_proxy, 16,
                  sizeof(::RooGaussStepBernstein) );
      instance.SetNew(&new_RooGaussStepBernstein);
      instance.SetNewArray(&newArray_RooGaussStepBernstein);
      instance.SetDelete(&delete_RooGaussStepBernstein);
      instance.SetDeleteArray(&deleteArray_RooGaussStepBernstein);
      instance.SetDestructor(&destruct_RooGaussStepBernstein);
      instance.SetStreamerFunc(&streamer_RooGaussStepBernstein);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RooGaussStepBernstein*)
   {
      return GenerateInitInstanceLocal((::RooGaussStepBernstein*)0);
   }
   // Static variable to force the class initialization
     static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RooGaussStepBernstein*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr RooGaussStepBernstein::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *RooGaussStepBernstein::Class_Name()
{
   return "RooGaussStepBernstein";
}

//______________________________________________________________________________
const char *RooGaussStepBernstein::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooGaussStepBernstein*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int RooGaussStepBernstein::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooGaussStepBernstein*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *RooGaussStepBernstein::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooGaussStepBernstein*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *RooGaussStepBernstein::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooGaussStepBernstein*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void RooGaussStepBernstein::Streamer(TBuffer &R__b)
{
   // Stream an object of class RooGaussStepBernstein.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      RooAbsPdf::Streamer(R__b);
      _x.Streamer(R__b);
      _mean.Streamer(R__b);
      _sigma.Streamer(R__b);
      _stepVal.Streamer(R__b);
      _coefList.Streamer(R__b);
      R__b >> rangeratio;
      R__b.CheckByteCount(R__s, R__c, RooGaussStepBernstein::IsA());
   } else {
      R__c = R__b.WriteVersion(RooGaussStepBernstein::IsA(), kTRUE);
      RooAbsPdf::Streamer(R__b);
      _x.Streamer(R__b);
      _mean.Streamer(R__b);
      _sigma.Streamer(R__b);
      _stepVal.Streamer(R__b);
      _coefList.Streamer(R__b);
      R__b << rangeratio;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_RooGaussStepBernstein(void *p) {
      return  p ? new(p) ::RooGaussStepBernstein : new ::RooGaussStepBernstein;
   }
   static void *newArray_RooGaussStepBernstein(Long_t nElements, void *p) {
      return p ? new(p) ::RooGaussStepBernstein[nElements] : new ::RooGaussStepBernstein[nElements];
   }
   // Wrapper around operator delete
   static void delete_RooGaussStepBernstein(void *p) {
      delete ((::RooGaussStepBernstein*)p);
   }
   static void deleteArray_RooGaussStepBernstein(void *p) {
      delete [] ((::RooGaussStepBernstein*)p);
   }
   static void destruct_RooGaussStepBernstein(void *p) {
      typedef ::RooGaussStepBernstein current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_RooGaussStepBernstein(TBuffer &buf, void *obj) {
      ((::RooGaussStepBernstein*)obj)->::RooGaussStepBernstein::Streamer(buf);
   }
} // end of namespace ROOT for class ::RooGaussStepBernstein

namespace {
  void TriggerDictionaryInitialization_RooGaussStepBernstein_Impl() {
    static const char* headers[] = {
"RooGaussStepBernstein.h",
"RooGaussStepBernstein.h",
0
    };
    static const char* includePaths[] = {
"/usr/include/root",
"/afs/cern.ch/user/f/fanx/CMSSW_12_6_0_patch1/src/HtoZg_fitting/Utilities/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "RooGaussStepBernstein dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$RooGaussStepBernstein.h")))  RooGaussStepBernstein;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "RooGaussStepBernstein dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "RooGaussStepBernstein.h"
#include "RooGaussStepBernstein.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"RooGaussStepBernstein", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("RooGaussStepBernstein",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_RooGaussStepBernstein_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_RooGaussStepBernstein_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_RooGaussStepBernstein() {
  TriggerDictionaryInitialization_RooGaussStepBernstein_Impl();
}
