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

#include "TJwfTensor.h"

class TSpinWaveFunction {
  
 private:
  Int_t J;
  Int_t max_pzm;
  Int_t *pot3;
  Int_t *mi;
  Int_t *M;
  char type;      // 's' Spin, 'l' Orbital Angular Momentum
                  // 'c' Spin conjugated
  TFracNum *coeff;

 public:
  TSpinWaveFunction() {
    J=0;
    max_pzm=1;
    pot3=0; mi=0; M=0;
  }
  
  TSpinWaveFunction(Int_t J, char type);
  TTensorSum* GetTensorSum(char name, Int_t delta);
  TTensorSum* GetSpinCoupledTensorSum(TSpinWaveFunction*, Int_t, Int_t);
  
  Int_t CheckCGFormula();
  
  //ClassDef(TSpinWaveFunction,1);
  
};
