#include "TLSAmpl.h"

#include <iostream>
#include <string>


using namespace std;


long debugLSAmpl=1;


TLSAmpl::TLSAmpl(long RankS1, long RankS2,
		long RankL,  long RankJ,
		long delta_,  long S_L,
		long cPsiInt, long cPsiChi,
		long cChiPhi, long cPsiPhi,
		long cPhiOme, long cChiOme,
		long cPhiEps, long cChiEps,
		long cNum) {

	J=RankJ;
	L=RankL;
	S=S_L;
	delta=delta_;
	ContractionNumber=cNum;

	cPsI=cPsiInt; cPsC=cPsiChi;
	cCP=cChiPhi;  cPsP=cPsiPhi;
	cPO=cPhiOme;  cCO=cChiOme;
	cPE=cPhiEps;  cCE=cChiEps;

	long totalR=RankS1+RankS2+RankL+RankJ;
	long contractions=2*(cPsiInt+cPsiChi+cChiPhi+cPsiPhi);
	long even_contraction=1;
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
		long ires = TSS->SpinInnerContraction(cPsiInt);
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
