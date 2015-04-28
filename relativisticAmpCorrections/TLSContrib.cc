#include "TLSContrib.h"

#include <iostream>

using namespace std;

long debugLSContrib=0;

TLSContrib::TLSContrib(TLSContrib *b, bool particle_exchange) {
	J=b->J; L=b->L; S=b->S; cNum=b->cNum; delta=b->delta;
	Nterms=b->Nterms;
	NormFactor=b->NormFactor;
	PureRelativistic=b->PureRelativistic;
	termFracNum = new TFracNum[Nterms];
	termg1pot = new long[Nterms];
	termg2pot = new long[Nterms];
	for (long i=0; i<Nterms; i++) {
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


TLSContrib::TLSContrib(TLSAmpl *A, long delta_, TFracNum scfac) {
	J=A->GetJ();
	L=A->GetL();
	S=A->GetS();
	cNum=A->GetContraction();
	delta=delta_;
	SpinCG=scfac;

	bool sign_reversal=false;
	if (delta<0 && (L+S-J)%2 ) sign_reversal=true;

	Nterms=A->GetNterms();
	termFracNum = new TFracNum[Nterms];
	termg1pot   = new long[Nterms];
	termg2pot   = new long[Nterms];

	NormFactor=TFracNum_Zero;
	if (debugLSContrib) cout <<"("<<J<<")"<<L<<S<<"["<<delta<<"]"<<endl;
	for (long i=0; i<Nterms; i++) {
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
		for (long i=0; i<Nterms; i++) {
			termFracNum[i]=termFracNum[i]*NormInv; // InvAbs;
		}
	}
}


long
TLSContrib::Add(TLSContrib *b, bool particle_exchange) {
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

	for (long i=0; i<Nterms; i++) {
		termFracNum[i]= TFracNum_Quarter * NormFactor *termFracNum[i];
	}

	for (long ib=0; ib<b->Nterms; ib++) {
		bool term_summed=false;
		for (long i=0; i<Nterms; i++) {
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
			long *new_termg1pot = new long[Nterms];
			long *new_termg2pot = new long[Nterms];
			for (long i=0; i<Nterms-1; i++) {
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
	long non_zeros=0;
	for (long i=0; i<Nterms; i++)
		if (! (termFracNum[i]==TFracNum_Zero) ) non_zeros++;

	if (non_zeros==0) {
		Nterms=0;
		return 0;
	}
	else {
		TFracNum *new_termFracNum = new TFracNum[non_zeros];
		long *new_termg1pot = new long[non_zeros];
		long *new_termg2pot = new long[non_zeros];

		long j=0;
		for (long i=0; i<Nterms; i++)
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
	for (long i=0; i<Nterms; i++) {
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
		for (long i=0; i<Nterms; i++) {
			termFracNum[i]=termFracNum[i]*NormInv;
		}
	}
	return Nterms;
}


long
TLSContrib::Print() {
	if (cNum==1) cout <<"g"<<"["<<cNum<<"] (";
	if (cNum==2) cout <<"f"<<"["<<cNum<<"] (";
	if (cNum>=3) cout <<"h"<<"["<<cNum<<"] (";
	cout <<J<<")"<<L<<S<<"( ";
	if (!PureRelativistic) {
		cout << NormFactor.FracStringSqrt() << " ) ( ";
	}
	for (long iT=0; iT<Nterms; iT++){
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


long
TLSContrib::PrintNR() {
	cout << NormFactor.FracStringSqrt();
	if (cNum==1) cout <<"*g";
	if (cNum==2) cout <<"*f";
	if (cNum==3) cout <<"*h";
	return 0;
}


long
TLSContrib::PrintNRG(TFracNum m) {
	cout << (NormFactor * m).FracStringSqrt();
	if (cNum==1) cout <<"*g";
	if (cNum==2) cout <<"*f";
	if (cNum==3) cout <<"*h";
	return 0;
}
