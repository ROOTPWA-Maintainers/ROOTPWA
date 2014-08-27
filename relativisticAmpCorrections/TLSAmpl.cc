#include <iostream>
#include <string>
#include "TLSAmpl.h"
using namespace std;

extern ClebschGordanBox box;

Int_t debugLSAmpl=1;

TLSAmpl::TLSAmpl(Int_t RankS1, Int_t RankS2,
		 Int_t RankL,  Int_t RankJ,
		 Int_t delta_,  Int_t S_L,
		 Int_t cPsiInt, Int_t cPsiChi, 
		 Int_t cChiPhi, Int_t cPsiPhi,
		 Int_t cPhiOme, Int_t cChiOme,
		 Int_t cPhiEps, Int_t cChiEps,
		 Int_t cNum) {
  
  J=RankJ;
  L=RankL;
  S=S_L;
  delta=delta_;
  ContractionNumber=cNum;
  
  cPsI=cPsiInt; cPsC=cPsiChi;
  cCP=cChiPhi;  cPsP=cPsiPhi;
  cPO=cPhiOme;  cCO=cChiOme; 
  cPE=cPhiEps;  cCE=cChiEps;
  
  Int_t totalR=RankS1+RankS2+RankL+RankJ;
  Int_t contractions=2*(cPsiInt+cPsiChi+cChiPhi+cPsiPhi);
  Int_t even_contraction=1;
  if ( totalR % 2 )  even_contraction=0;

  if ( ( even_contraction==1 && contractions != totalR   ) ||
       ( even_contraction==0 && contractions != totalR-3 ) ) {
    cerr << "Invalid contraction occurred" << endl; 
    return;
  }

  if (debugLSAmpl) {
    cout << "LSAmpl: "<<RankS1<<" "<<RankS2<<" L="  
	 <<RankL<<" "<<RankJ<<" d="<<delta<<" S="<<S_L<<" c: "  
	 <<cPsiInt<<" "<<cPsiChi<<" "<<cChiPhi<<" "<<cPsiPhi << " s: "
	 <<cPhiOme<<" "<<cChiOme<<" "<<cPhiEps<<" "<<cChiEps;
    if (even_contraction) cout << endl;
    else cout<< " (iw)" << endl;  
  }
  
  TSpinWaveFunction WFS1(RankS1, 's');
  TSpinWaveFunction WFS2(RankS2, 's');
  
  TTensorSum *TSS = WFS1.GetSpinCoupledTensorSum(&WFS2, delta, S_L);
  if (debugLSAmpl==2) TSS->Print();

  if (cPsiInt) {
    Int_t ires = TSS->SpinInnerContraction(cPsiInt);
    if (ires==0) {
      if (debugLSAmpl) cout << "Inner contraction is zero." << endl;
      Nterms=0;
      return;
    }
    if (debugLSAmpl==2) TSS->Print();
  }
  
  TSpinWaveFunction WFL (RankL,  'l');
  TTensorSum *TSL = WFL.GetTensorSum('c', 0);
  if (debugLSAmpl==2) TSL->Print();
  
  TTensorSum *TSLS = TSS->LSContraction(TSL, cPsiChi, cChiOme, cChiEps, 'c');
  if (TSLS->GetNterms()==0) {
    if (debugLSAmpl) cout << "LS contraction is zero." << endl;
    Nterms=0;
    return;
  }
  if (debugLSAmpl==2) TSLS->Print();
  
  TSpinWaveFunction WFJ (RankJ,  'c');
  TTensorSum *TSJ = WFJ.GetTensorSum('p', delta);
  if (debugLSAmpl==2) TSJ->Print();
  
  TTensorSum *TSLSJ = TSLS->LSContraction(TSJ, cPsiPhi, cPhiOme, cPhiEps, 'p');

  if (TSLSJ->GetNterms()==0) {
    if (debugLSAmpl) cout << "JS contraction is zero." << endl;
    Nterms=0;
    return;
  }

  if (debugLSAmpl==2) TSLSJ->Print();
  
  TSScalar = TSLSJ->LJContraction(cChiPhi, even_contraction);
  
  if (debugLSAmpl==2) TSLSJ->Print();
  
  Nterms=TSScalar->GetNterms();
  
  if (Nterms==0) {
    if (debugLSAmpl) cout << "LJ contraction is zero." << endl;
    return;
  } 
  
  if (debugLSAmpl) {
    //TSScalar->Print();
    TSScalar->Print('s');
  }
}

Int_t debugLSContrib=0;

TLSContrib::TLSContrib(TLSContrib *b, Bool_t particle_exchange) {
  J=b->J; L=b->L; S=b->S; cNum=b->cNum; delta=b->delta;
  Nterms=b->Nterms;
  NormFactor=b->NormFactor;
  PureRelativistic=b->PureRelativistic;
  termFracNum = new TFracNum[Nterms];
  termg1pot = new Int_t[Nterms];
  termg2pot = new Int_t[Nterms];
  for (Int_t i=0; i<Nterms; i++) {
    // force new power fields by multiplication " *1 "
    termFracNum[i] = TFracNum_One * b->termFracNum[i]; 
    if (particle_exchange) {
      termg1pot[i]   = b->termg2pot[i];
      termg2pot[i]   = b->termg1pot[i];
    }
    else {
      termg1pot[i]   = b->termg1pot[i];
      termg2pot[i]   = b->termg2pot[i];
    }
  }
};

TLSContrib::TLSContrib(TLSAmpl *A, Int_t delta_, TFracNum scfac) {
  J=A->GetJ();
  L=A->GetL();
  S=A->GetS();
  cNum=A->GetContraction();
  delta=delta_;
  SpinCG=scfac;
  
  Bool_t sign_reversal=false;
  if (delta<0 && (L+S-J)%2 ) sign_reversal=true;  
  
  Nterms=A->GetNterms();
  termFracNum = new TFracNum[Nterms];
  termg1pot   = new Int_t[Nterms];
  termg2pot   = new Int_t[Nterms];
  
  NormFactor=TFracNum_Zero;
  if (debugLSContrib) cout <<"("<<J<<")"<<L<<S<<"["<<delta<<"]"<<endl;
  for (Int_t i=0; i<Nterms; i++) {
    TTensorTerm *tt = A->GetTerm(i); //TSScalar->GetTerm(i);
    
    if (debugLSContrib)
      cout << scfac.FracStringSqrt()<< ","
	   << tt->GetPreFac().FracStringSqrt()
	   << " L"<<L<<"S"<<S<<"J"<<J<<"Ampldelta"<<A->Getdelta()
	   <<" delta"<<delta<<", sigrev "<<sign_reversal;
    
    termFracNum[i] = scfac * tt->GetPreFac();
    if (sign_reversal==true) termFracNum[i].FlipSign();
    //TFracNum tSqr = NormFactor * termFracNum[i]; 
    //if ( ! tSqrt.Sqrt() )
    //  cout << " Building LSContrib leads not to squared-fractional numbers,"
    //	   << " results will be wrong!" << endl;

    TFracNum *sum = NormFactor.SumSignedRoots(&(termFracNum[i]));
    if (sum) { NormFactor = *sum; }
    else     { 
      cerr << "TLSContrib: Normalisation not root-fractional,"
	   << " *** results will be wrong *** " << endl; 
    }
    
    //    NormFactor=NormFactor+termFracNum[i]+TFracNum_Two*tSqrt;

    termg1pot[i] = tt->GetGamS();
    termg2pot[i] = tt->GetGamSig();
    if (debugLSContrib) cout << " -> Normfactor: " 
			     << NormFactor.FracStringSqrt() << endl;
    
  }
  if (NormFactor==TFracNum_Zero) {
    PureRelativistic=true;
    NormFactor=TFracNum_One;
  }
  else {
    PureRelativistic=false;
    TFracNum NormInv=NormFactor;
    NormInv.Invert();
    // TFracNum InvAbs=TFracNum_One * NormInv; //Bug: real copy forced by "1 *"
    // InvAbs.Abs();
    for (Int_t i=0; i<Nterms; i++) {
      termFracNum[i]=termFracNum[i]*NormInv; // InvAbs;
    }    
  }
}

Int_t 
TLSContrib::Add(TLSContrib *b, Bool_t particle_exchange) {
  if (J!=b->J || L!=b->L || S!=b->S) {
    cerr <<"TLSContrib::Add : Something is wrong, trying to add different"
	 <<" (J;L,S): ("<<J<<";"<<L<<","<<S<<") != ("
	 <<b->J<<";"<<b->L<<","<<b->S<<")"<<endl;
    return -1;
  }

  //
  // Include normalisation factor and the factor (1/2)**2 in the squared
  // representation of prefactors
  //

  for (Int_t i=0; i<Nterms; i++) {
    termFracNum[i]= TFracNum_Quarter * NormFactor *termFracNum[i];
  }
  
  for (Int_t ib=0; ib<b->Nterms; ib++) {
    Bool_t term_summed=false;
    for (Int_t i=0; i<Nterms; i++) {
      if ( !term_summed && cNum == b->cNum &&
	   ( (particle_exchange==true  && 
	      termg1pot[i]==b->termg2pot[ib] && 
	      termg2pot[i]==b->termg1pot[ib]     ) ||
	     (particle_exchange==false && 
	      termg1pot[i]==b->termg1pot[ib] && 
	      termg2pot[i]==b->termg2pot[ib]     ) )  ) {
	term_summed=true;
	TFracNum bterm=TFracNum_Quarter * b->NormFactor * b->termFracNum[ib];

	if (J%2) bterm.FlipSign();

	TFracNum *sum = bterm.SumSignedRoots(&(termFracNum[i]));
	if (sum) { termFracNum[i] = *sum; }
	else     { 
	  cerr << "TLSContrib: Normalisation not root-fractional,"
	       << " *** results will be wrong *** " << endl; 
	}
      }
    }
    if (!term_summed) {
      Nterms++;
      TFracNum *new_termFracNum = new TFracNum[Nterms];
      Int_t *new_termg1pot = new Int_t[Nterms];
      Int_t *new_termg2pot = new Int_t[Nterms];
      for (Int_t i=0; i<Nterms-1; i++) {
	new_termFracNum[i]=termFracNum[i];
	new_termg1pot[i]=termg1pot[i];
	new_termg2pot[i]=termg2pot[i];
      }
      new_termFracNum[Nterms-1] = 
	TFracNum_Quarter * b->NormFactor * b->termFracNum[ib];
      //if ( ! new_termFracNum[Nterms-1].Sqrt() )
      //  cout << "Square root not possible, this will lead to wrong results:"
      //       << new_termFracNum[Nterms-1].FracString() << endl;
      //new_termFracNum[Nterms-1]=
      //	TFracNum_Half * new_termFracNum[Nterms-1] * b->NormFactor;
      if (J%2) new_termFracNum[Nterms-1].FlipSign();
      if ( particle_exchange ) {
	new_termg1pot[Nterms-1] = b->termg2pot[ib];
	new_termg2pot[Nterms-1] = b->termg1pot[ib];
      }
      else {
	new_termg1pot[Nterms-1] = b->termg1pot[ib];
	new_termg2pot[Nterms-1] = b->termg2pot[ib];
      }
      termFracNum=new_termFracNum;
      termg1pot=new_termg1pot;
      termg2pot=new_termg2pot;
    }
  }

  //
  // Eliminate zero entries
  //
  Int_t non_zeros=0;
  for (Int_t i=0; i<Nterms; i++)
    if (! (termFracNum[i]==TFracNum_Zero) ) non_zeros++;

  if (non_zeros==0) {
    Nterms=0;
    return 0;
  }
  else {
    TFracNum *new_termFracNum = new TFracNum[non_zeros];
    Int_t *new_termg1pot = new Int_t[non_zeros];
    Int_t *new_termg2pot = new Int_t[non_zeros];
    
    Int_t j=0;
    for (Int_t i=0; i<Nterms; i++)
      if (! (termFracNum[i]==TFracNum_Zero) ) {
	new_termFracNum[j]=termFracNum[i];
	new_termg1pot[j]=termg1pot[i];
	new_termg2pot[j]=termg2pot[i];
	j++;
      }
    Nterms=non_zeros;
    termFracNum=new_termFracNum;
    termg1pot=new_termg1pot;
    termg2pot=new_termg2pot;
  }
  
  //
  // Recalculate Normalization Factor
  //
  NormFactor=TFracNum_Zero;
  for (Int_t i=0; i<Nterms; i++) {
    TFracNum *sum = NormFactor.SumSignedRoots(&(termFracNum[i]));
    if (sum) { NormFactor = *sum; }
    else     { 
      cerr << "TLSContrib: Normalisation not root-fractional,"
	   << " *** results will be wrong *** " << endl; 
    }
  }
  
  //
  // Apply normalization
  //
  if (NormFactor==TFracNum_Zero) {
    PureRelativistic=true;
    NormFactor=TFracNum_One;
  }
  else {
    PureRelativistic=false;
    TFracNum NormInv=NormFactor;
    NormInv.Invert();
    for (Int_t i=0; i<Nterms; i++) {
      termFracNum[i]=termFracNum[i]*NormInv;
    }    
  }
  return Nterms;
}

Int_t 
TLSContrib::Print() {
  if (cNum==1) cout <<"g"<<"["<<cNum<<"] (";
  if (cNum==2) cout <<"f"<<"["<<cNum<<"] (";
  if (cNum>=3) cout <<"h"<<"["<<cNum<<"] (";
  cout <<J<<")"<<L<<S<<"( ";
  if (!PureRelativistic) {
    cout << NormFactor.FracStringSqrt() << " ) ( ";
  }
  for (Int_t iT=0; iT<Nterms; iT++){
    cout << termFracNum[iT].FracStringSqrt()<<" ";
    if (termg1pot[iT]) {
      if (termg1pot[iT]==1) cout << " gs";
      else              cout << " gs^" << termg1pot[iT];
    }
    if (termg2pot[iT]) {
      if (termg2pot[iT]==1) cout << " gsig";
      else                cout << " gsig^" << termg2pot[iT];
    }
  }
  cout <<" )"<< endl;
  return 0;
}

Int_t
TLSContrib::PrintNR() {
  cout << NormFactor.FracStringSqrt();
  if (cNum==1) cout <<"*g";
  if (cNum==2) cout <<"*f";
  if (cNum==3) cout <<"*h";
  return 0;
}

Int_t
TLSContrib::PrintNRG(TFracNum m) {
  cout << (NormFactor * m).FracStringSqrt();
  if (cNum==1) cout <<"*g";
  if (cNum==2) cout <<"*f";
  if (cNum==3) cout <<"*h";
  return 0;
}


TLSNonRel::TLSNonRel(TLSContrib *C) {
  J=C->GetJ();
  L=C->GetL();
  S=C->GetS();
  Nterms=1;
  RelLS = new TLSContrib*[1];
  RelLS[0]=C;
  TFracNum *JdL0Sd = box.GetCG(J, L, S);
  //cout << "delta=" << C->GetDelta()
  //     << ",S=" << S
  //     << ", 2L+1/2J+1=" << TFracNum(2*L+1,2*J+1).FracStringSqrt()
  //     << ",CG(JdL0Sd)=" 
  //     << JdL0Sd[CGIndex(L, 0, S, C->GetDelta())].FracStringSqrt()
  //     << ", SpinCG=" << C->GetSpinCG()->FracStringSqrt()
  //     << endl;
  GnrPrefac = TFracNum(2*L+1,2*J+1) * 
    JdL0Sd[CGIndex(L, 0, S, C->GetDelta())] *
    *C->GetSpinCG();
}

Int_t TLSNonRel::Add(TLSContrib *C){
  if ( ! CheckJLS(C) ) {
    cout << "TLSNonRel::Add not appropriate, skipping (check code)!!!" << endl;
    return -1;
  }
  Nterms++;
  TLSContrib **newRelLS = new TLSContrib*[Nterms];  
  for (Int_t i=0; i<Nterms-1; i++) {
    newRelLS[i]=RelLS[i];
  }
  newRelLS[Nterms-1]=C;
  RelLS=newRelLS;
  return 0;
}

Int_t TLSNonRel::Print() {
  cout << " [ " << GnrPrefac.FracStringSqrt() << "  G_"<< L << S <<" ] ";
  for (Int_t i=0; i<Nterms; i++) {
    RelLS[i]->PrintNR();
  }
  cout << endl;
  return 0;
}

Int_t TLSNonRel::PrintG() {
  cout << " [ G_"<< L << S <<" ] ";
  for (Int_t i=0; i<Nterms; i++) {
    TFracNum GnrInv(GnrPrefac);
    GnrInv.Invert();
    RelLS[i]->PrintNRG(GnrInv);
  }
  cout << endl;
  return 0;
}
