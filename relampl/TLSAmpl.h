#ifndef TLSAmpl_h
#define TLSAmpl_h
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

#include "TSpinWaveFunction.h"

/*!
  \class TLSAmpl
  \brief Relativistic LS-coupling amplitudes


  \author Jan.Friedrich@ph.tum.de
*/
class TLSAmpl {
  
 private:
  Int_t J;
  Int_t L;
  Int_t S;
  Int_t delta;

  Int_t ContractionNumber;
  Int_t cPsI; Int_t cPsC;
  Int_t cCP;  Int_t cPsP;
  Int_t cPO; Int_t cCO; 
  Int_t cPE; Int_t cCE;

  Int_t Nterms;
  TTensorSum *TSScalar;
  
 public:
  
  TLSAmpl() {
    J=0; L=0; S=0; delta=0;
    Nterms=0;
  };
  
  TLSAmpl(Int_t RankS1, Int_t RankS2,
	  Int_t RankL,  Int_t RankJ,
	  Int_t delta,  Int_t S_L,
	  Int_t cPsiInt, Int_t cPsiChi, 
	  Int_t cChiPhi, Int_t cPsiPhi,
	  Int_t cPhiOme, Int_t cChiOme,
	  Int_t cPhiEps, Int_t cChiEps,
	  Int_t cNum);
  
  Int_t GetNterms() {return Nterms;};
  
  Int_t GetJ()     {return J;};
  Int_t GetL()     {return L;};
  Int_t GetS()     {return S;};
  Int_t Getdelta() {return delta;};
  Int_t GetContraction() {return ContractionNumber;};  
  Bool_t CheckContraction(Int_t L_, Int_t S_, Int_t cPsI_, Int_t cPsC_,
			  Int_t cCP_, Int_t cPsP_, Int_t cPO_, Int_t cCO_, 
			  Int_t cPE_, Int_t cCE_) {
    if ( L!=L_ || S!=S_ || cPsI!=cPsI_ || cPsC!=cPsC_ || cCP!=cCP_ || 
	 cPsP!=cPsP_ || cPO!=cPO_ || cCO!=cCO_ || cPE!=cPE_ || cCE!=cCE_ )
      return false;
    return true;
  }  
  TTensorTerm *GetTerm(Int_t i) {return TSScalar->GetTerm(i);};
  
  //ClassDef(TLSAmpl,1);
  
};

/*!
  \class TLSContrib
  \brief Relativistic LS-coupling contributions


  \author Jan.Friedrich@ph.tum.de
*/
class TLSContrib {

 private:
  Int_t J;
  Int_t L;
  Int_t S;
  Int_t cNum;
  Int_t delta;
  TFracNum SpinCG;
  
  Int_t Nterms;
  TFracNum NormFactor;     // Square  of normalisation factor
  TFracNum *termFracNum;   // Squares of prefactors
  Int_t    *termg1pot;     // exponent of gamma_s
  Int_t    *termg2pot;     // exponent of gamma_sigma
  
  Bool_t PureRelativistic;

 public:
  
  TLSContrib() {
    J=0; L=0; S=0; delta=0;
    Nterms=0;
    termFracNum=0;
    termg1pot=0;
    termg2pot=0;
  };
  
  TLSContrib(TLSContrib *, Bool_t);
  TLSContrib(TLSAmpl*, Int_t, TFracNum);
  
  Bool_t SameParameter(TLSContrib* b) {
    if (J==b->J && L==b->L && S==b->S && cNum==b->cNum)
      return true;
    return false;
  };
  Int_t GetNterms() {return Nterms;};
  Int_t Add(TLSContrib*, Bool_t);
  Int_t Print();
  Int_t PrintNR();
  Int_t PrintNRG(TFracNum);
  Bool_t IsPureRelativistic() {return PureRelativistic;};
  Int_t GetJ() {return J;};
  Int_t GetL() {return L;};
  Int_t GetS() {return S;};
  Int_t GetDelta() {return delta;};
  Int_t GetRunningNumber() {return cNum;};
  TFracNum* GetSpinCG() {return &SpinCG;};
  TFracNum* GetNormFactor() {return &NormFactor;};
  
  TFracNum* GetTermFracNum() {return termFracNum;};
  Int_t*    GetTermg1pot()   {return termg1pot;}; // exponent of gamma_s
  Int_t*    GetTermg2pot()   {return termg2pot;}; // exponent of gamma_s
  
  //ClassDef(TLSContrib,1);
  
};


/*!
  \class TLSNonRel
  \brief Non-relativistic LS-coupling contributions


  \author Jan.Friedrich@ph.tum.de
*/
class TLSNonRel {

 private:
  Int_t J;
  Int_t L;
  Int_t S;
  Int_t Nterms;
  TLSContrib* *RelLS;
  TFracNum GnrPrefac;

 public:
  TLSNonRel(TLSContrib *C);
  
  Int_t CheckJLS(TLSContrib *C) {
    return (C->GetJ()==J && C->GetL()==L && C->GetS()==S) ? 1 : 0;
  };
  Int_t Add(TLSContrib *C);
  Int_t GetL(){return L;};
  Int_t GetS(){return S;};
  Int_t Print();
  Int_t PrintG();

  //ClassDef(TLSNonRel,1);
};
#endif
