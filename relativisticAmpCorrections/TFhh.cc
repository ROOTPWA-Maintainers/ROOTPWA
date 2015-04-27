#include <iostream>
#include <fstream>
#include <string>
#include "TFhh.h"
using namespace std;

extern ClebschGordanBox box;

Int_t debugFhh=1;

TFhh::TFhh(Int_t J_, Int_t S1, Int_t S2,
		Int_t lambda_, Int_t nu_,
		const vector<TLSAmpl*>& LSampl,
		bool even_contr_) {

	J=J_;
	lambda= lambda_;
	nu= nu_;
	even_contraction = even_contr_;
	Int_t delta=lambda-nu;
	Nterms=0;

	name_str=new char[15];
	sprintf(name_str, "F_%lld_%lld", lambda, nu);

	for (unsigned int iLS=0; iLS<LSampl.size(); iLS++) {
		if (LSampl[iLS]->Getdelta() ==  delta ||
				LSampl[iLS]->Getdelta() == -delta ) {
			// cout << "iLS=" << iLS << ", delta=" << delta << endl;
			TFracNum *CG3S = box.GetCG(LSampl[iLS]->GetS(), S1, S2);
			TFracNum SpinCouplFac = CG3S[CGIndex(S1, lambda, S2, -nu)];
			if (SpinCouplFac==TFracNum_Zero) {
				if (debugFhh==2)
					cout << "Clebsch-Gordan is zero" << endl;
			}
			else {
				TLSContrib* *newPTR = new TLSContrib*[Nterms+1];
				for (Int_t iC=0; iC<Nterms; iC++) newPTR[iC] = LSt[iC];
				newPTR[Nterms] = new TLSContrib(LSampl[iLS], delta, SpinCouplFac);
				// delete[] LSt;
				Nterms++;
				LSt=newPTR;
			}
		}
	}
	if (debugFhh) Print();
}

TFhh::TFhh(TFhh *sFhh, char flag) {
	if (flag != 'i' && flag !='m'){
		cerr << "TFhh::TFhh unknown flag "<<flag<<endl;
		return;
	}
	if (debugFhh)
		cout << "Initializing from single Amplitude" << endl;
	if ( ( flag =='i' && (sFhh->GetJ()) % 2      ) ||
			( flag =='m' && (sFhh->GetJ()) % 2 == 0 )    ) {
		cout << sFhh->GetName() << "[symm] = 0" << endl;
		Nterms=0;
	}
	else {
		name_str = new char[15];
		sprintf(name_str, "%s[symm]", sFhh->GetName());
		J        = sFhh->GetJ();
		lambda   = sFhh->GetLambda();
		nu       = sFhh->GetNu();
		even_contraction = sFhh->GetEvenContr();
		Nterms   = sFhh->GetNterms();
		LSt      = sFhh->GetLStPtr();
		Print();
	}
}

TFhh::TFhh(TFhh *sFhh, TFhh *xFhh) {
	name_str = new char[15];
	sprintf(name_str, "%s[symm]", sFhh->GetName());
	J        = sFhh->GetJ();
	lambda   = sFhh->GetLambda();
	nu       = sFhh->GetNu();
	even_contraction = sFhh->GetEvenContr();
	Nterms=0;
	LSt   =0;
	if (J != xFhh->GetJ() || even_contraction != xFhh->GetEvenContr() ||
			lambda != xFhh->GetNu() || nu != xFhh->GetLambda() ){
		cerr << "TFhh::TFhh(TFhh *, TFhh*): Something is wrong,"    << endl
			<< " source amplitudes for symmetrization do not match" << endl;
		return;
	}

	// Since some LS-terms cancel out, first a max-length array of
	// LSContrib pointers is filled, and then squeezed to only
	// non-zero contributions

	Int_t Ns = sFhh->GetNterms();
	Int_t Nx = xFhh->GetNterms();
	Nterms = Ns + Nx;
	TLSContrib* *pLSt = new TLSContrib*[Nterms];

	Int_t prevterms=0;

	for (Int_t i=0; i<Nterms; i++) {
		TLSContrib *to_be_added;
		if (i<Ns) to_be_added = sFhh->GetLStPtr()[i];
		else      to_be_added = xFhh->GetLStPtr()[i-Ns];
		Int_t found_in_prevterm=0;
		for (Int_t j=0; j<prevterms; j++) {
			if ( pLSt[j]->SameParameter(to_be_added) ) {
				found_in_prevterm=1;
				if (i<Ns) pLSt[j]->Add(to_be_added, false);
				else      pLSt[j]->Add(to_be_added, true);
			}
		}
		if (!found_in_prevterm) {
			if (i<Ns) pLSt[prevterms] = new TLSContrib(to_be_added, false);
			else      pLSt[prevterms] = new TLSContrib(to_be_added, true);
			prevterms++;
		}
	}

	//
	// Cancel zeros
	//
	Int_t non_zeros=0;
	for (Int_t i=0; i<prevterms; i++) {
		if ( pLSt[i]->GetNterms() != 0 ) non_zeros++;
	}

	Nterms=non_zeros;

	if (Nterms) {
		LSt = new TLSContrib*[Nterms];
		Int_t j=0;
		for (Int_t i=0; i<prevterms; i++) {
			if ( pLSt[i]->GetNterms() != 0 ) {
				LSt[j]=new TLSContrib(pLSt[i], false);
				j++;
			}
		}
	}
	Print();

}


Int_t TFhh::NonRelLimit(){
	NNRterms=0;
	for (Int_t i=0; i<Nterms; i++){
		if (! LSt[i]->IsPureRelativistic()) {
			Int_t jfound=-1;
			for (Int_t j=0; j<NNRterms; j++) {
				if (NRLSt[j]->CheckJLS(LSt[i])){
					jfound=j;
					break;
				}
			}
			if (jfound != -1) {
				NRLSt[jfound]->Add(LSt[i]);
			}
			else {
				NNRterms++;
				TLSNonRel **newNRLS=new TLSNonRel*[NNRterms];
				for (Int_t j=0; j<NNRterms-1; j++) {
					newNRLS[j]=NRLSt[j];
				}
				newNRLS[NNRterms-1] = new TLSNonRel(LSt[i]);
				NRLSt=newNRLS;
			}
		}
	}
	cout << name_str << " (NR) = "  << endl;
	for (Int_t j=0; j<NNRterms; j++) {
		cout << "LS=" << NRLSt[j]->GetL() << NRLSt[j]->GetS() << ": ";
		NRLSt[j]->Print();
	}
	return 0;
}

Int_t TFhh::PrintNRG(){
	cout << name_str << " (NR) = "  << endl;
	for (Int_t j=0; j<NNRterms; j++) {
		cout << "LS=" << NRLSt[j]->GetL() << NRLSt[j]->GetS() << " => ";
		NRLSt[j]->PrintG();
	}
	return 0;
}


Int_t
TFhh::Print() {
	cout << name_str << " =";
	if (even_contraction) cout << endl;
	else cout << " (iw)" << endl;
	for (Int_t iLSt=0; iLSt<Nterms; iLSt++) {
		LSt[iLSt]->Print();
	}
	if (Nterms==0) cout << " 0" << endl;
	return 0;
}
