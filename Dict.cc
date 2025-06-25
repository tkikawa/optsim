// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME Dict
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
#include "Gui.hh"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static TClass *Gui_Dictionary();
   static void Gui_TClassManip(TClass*);
   static void delete_Gui(void *p);
   static void deleteArray_Gui(void *p);
   static void destruct_Gui(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Gui*)
   {
      ::Gui *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Gui));
      static ::ROOT::TGenericClassInfo 
         instance("Gui", "Gui.hh", 22,
                  typeid(::Gui), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Gui_Dictionary, isa_proxy, 0,
                  sizeof(::Gui) );
      instance.SetDelete(&delete_Gui);
      instance.SetDeleteArray(&deleteArray_Gui);
      instance.SetDestructor(&destruct_Gui);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Gui*)
   {
      return GenerateInitInstanceLocal((::Gui*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::Gui*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Gui_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::Gui*)0x0)->GetClass();
      Gui_TClassManip(theClass);
   return theClass;
   }

   static void Gui_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrapper around operator delete
   static void delete_Gui(void *p) {
      delete ((::Gui*)p);
   }
   static void deleteArray_Gui(void *p) {
      delete [] ((::Gui*)p);
   }
   static void destruct_Gui(void *p) {
      typedef ::Gui current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::Gui

namespace {
  void TriggerDictionaryInitialization_Dict_Impl() {
    static const char* headers[] = {
"Gui.hh",
0
    };
    static const char* includePaths[] = {
"/cern/root/include/",
"/home/kikawa/optsim_test/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$Gui.hh")))  Gui;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "Dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "Gui.hh"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"Gui", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_Dict() {
  TriggerDictionaryInitialization_Dict_Impl();
}
