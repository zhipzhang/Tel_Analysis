// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dOdOdIsrcdIMonteCarlodict
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
#include "MonteCarloRunHeader.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_MonteCarloRunHeader(void *p = 0);
   static void *newArray_MonteCarloRunHeader(Long_t size, void *p);
   static void delete_MonteCarloRunHeader(void *p);
   static void deleteArray_MonteCarloRunHeader(void *p);
   static void destruct_MonteCarloRunHeader(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MonteCarloRunHeader*)
   {
      ::MonteCarloRunHeader *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::MonteCarloRunHeader >(0);
      static ::ROOT::TGenericClassInfo 
         instance("MonteCarloRunHeader", ::MonteCarloRunHeader::Class_Version(), "MonteCarloRunHeader.h", 14,
                  typeid(::MonteCarloRunHeader), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::MonteCarloRunHeader::Dictionary, isa_proxy, 4,
                  sizeof(::MonteCarloRunHeader) );
      instance.SetNew(&new_MonteCarloRunHeader);
      instance.SetNewArray(&newArray_MonteCarloRunHeader);
      instance.SetDelete(&delete_MonteCarloRunHeader);
      instance.SetDeleteArray(&deleteArray_MonteCarloRunHeader);
      instance.SetDestructor(&destruct_MonteCarloRunHeader);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MonteCarloRunHeader*)
   {
      return GenerateInitInstanceLocal((::MonteCarloRunHeader*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::MonteCarloRunHeader*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr MonteCarloRunHeader::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *MonteCarloRunHeader::Class_Name()
{
   return "MonteCarloRunHeader";
}

//______________________________________________________________________________
const char *MonteCarloRunHeader::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MonteCarloRunHeader*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int MonteCarloRunHeader::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MonteCarloRunHeader*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MonteCarloRunHeader::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MonteCarloRunHeader*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MonteCarloRunHeader::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MonteCarloRunHeader*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void MonteCarloRunHeader::Streamer(TBuffer &R__b)
{
   // Stream an object of class MonteCarloRunHeader.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(MonteCarloRunHeader::Class(),this);
   } else {
      R__b.WriteClassBuffer(MonteCarloRunHeader::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_MonteCarloRunHeader(void *p) {
      return  p ? new(p) ::MonteCarloRunHeader : new ::MonteCarloRunHeader;
   }
   static void *newArray_MonteCarloRunHeader(Long_t nElements, void *p) {
      return p ? new(p) ::MonteCarloRunHeader[nElements] : new ::MonteCarloRunHeader[nElements];
   }
   // Wrapper around operator delete
   static void delete_MonteCarloRunHeader(void *p) {
      delete ((::MonteCarloRunHeader*)p);
   }
   static void deleteArray_MonteCarloRunHeader(void *p) {
      delete [] ((::MonteCarloRunHeader*)p);
   }
   static void destruct_MonteCarloRunHeader(void *p) {
      typedef ::MonteCarloRunHeader current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::MonteCarloRunHeader

namespace {
  void TriggerDictionaryInitialization_MonteCarlodict_Impl() {
    static const char* headers[] = {
"MonteCarloRunHeader.h",
0
    };
    static const char* includePaths[] = {
"/data/home/zhipz/root/include/root",
"/data/home/zhipz/tools/include/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "MonteCarlodict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$MonteCarloRunHeader.h")))  MonteCarloRunHeader;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "MonteCarlodict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "MonteCarloRunHeader.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"MonteCarloRunHeader", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("MonteCarlodict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_MonteCarlodict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_MonteCarlodict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_MonteCarlodict() {
  TriggerDictionaryInitialization_MonteCarlodict_Impl();
}
