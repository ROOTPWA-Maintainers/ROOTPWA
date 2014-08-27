#ifndef TFhh_h
#define TFhh_h
/*!
  \class TFhh
  \brief Relativistic helicity-coupling structure

  The class describes relativistic helicity-coupling amplitudes, as
  described in the paper [PRD78(2008)], equation (6.6). It collects a field
  of terms TLSContrib with different LS. It also provides the relation to the
  non-relativistic Zeemach amplitudes.

  \author Jan.Friedrich@ph.tum.de
*/
//
// Uncomment the following line
// if you want to work in CINT (root.cern.ch)
//
//#define __JCINT__
#ifndef __JCINT__
#define Int_t    long long
#define Double_t double
#define Bool_t   bool
#endif

#include "TLSAmpl.h"

class TFhh {
  
 private:
  char* name_str;
  Int_t J;
  Int_t lambda;
  Int_t nu;
  Int_t even_contraction;
  Int_t Nterms;
  TLSContrib* *LSt;
  Int_t NNRterms;
  TLSNonRel* *NRLSt;
  
 public:
  
  TFhh(){
    J=0; lambda=0; nu=0; Nterms=0; LSt=0;
  };
  
  TFhh(Int_t, Int_t , Int_t, 
       Int_t, Int_t, Int_t, TLSAmpl**,
       Int_t);
  
  TFhh(TFhh*, char);
  TFhh(TFhh*, TFhh*);
  
  Int_t GetNterms() {return Nterms;};
  Int_t IsNuNu()      { if (lambda== nu) return 1; return 0; };
  Int_t IsNuMinusNu() { if (lambda==-nu) return 1; return 0; };
  Int_t GetLambda() { return lambda;};
  Int_t GetNu() { return nu;};
  Int_t GetJ() {return J;};
  Int_t GetEvenContr() {return even_contraction;};
  TLSContrib** GetLStPtr() {return LSt;};
  char* GetName() {return name_str;};
  
  Int_t NonRelLimit();
  Int_t PrintNRG();
  Int_t Print();
  //ClassDef(THCTerm,1);
  
};

class TJSS {
  
 private:

  Int_t JMother; Int_t etaJ;
  Int_t SDecay1; Int_t eta1;
  Int_t SDecay2; Int_t eta2;
  Int_t Smin;
  Int_t Smax;
  Int_t Lmin;
  Int_t Lmax;
  
  Int_t NLSAmpl;
  TLSAmpl* *LSAmplitudes;
  
  Int_t NFhhAmpl;
  TFhh* *FhhAmpl;

  Int_t NFhhIdAmpl;
  TFhh* *FhhIdAmpl;
  
 public:
  
  TJSS() {
    JMother = 0;  etaJ = 1;
    SDecay1 = 0;  eta1 = 1;
    SDecay2 = 0;  eta2 = 1;
    CalcAmpl();
  };
  
  TJSS(Int_t J,  Int_t eJ,
       Int_t S1, Int_t e1,
       Int_t S2, Int_t e2,
       Int_t option) {
    JMother = J;  etaJ=eJ;
    SDecay1 = S1; eta1=e1;
    SDecay2 = S2; eta2=e2;
    CalcAmpl();
    if (option>1) PrintHFILE();
  };
  
  Int_t CalcAmpl();
  Int_t PrintHFILE();
  
  //ClassDef(TJSS,1);
};

#endif
