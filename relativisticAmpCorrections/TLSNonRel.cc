#include "TLSNonRel.h"

#include <iostream>

using namespace std;

TLSNonRel::TLSNonRel(TLSContrib *C) {
	J = C->GetJ();
	L = C->GetL();
	S = C->GetS();
	Nterms = 1;
	RelLS = new TLSContrib*[1];
	RelLS[0] = C;
	TFracNum *JdL0Sd = ClebschGordanBox::instance()->GetCG(J, L, S);
	//cout << "delta=" << C->GetDelta()
	//     << ",S=" << S
	//     << ", 2L+1/2J+1=" << TFracNum(2*L+1,2*J+1).FracStringSqrt()
	//     << ",CG(JdL0Sd)="
	//     << JdL0Sd[CGIndex(L, 0, S, C->GetDelta())].FracStringSqrt()
	//     << ", SpinCG=" << C->GetSpinCG()->FracStringSqrt()
	//     << endl;
	GnrPrefac = TFracNum(2 * L + 1, 2 * J + 1)
			* JdL0Sd[ClebschGordanBox::CGIndex(L, 0, S, C->GetDelta())]
			* *C->GetSpinCG();
}

long TLSNonRel::Add(TLSContrib *C) {
	if (!CheckJLS(C)) {
		cout << "TLSNonRel::Add not appropriate, skipping (check code)!!!"
				<< endl;
		return -1;
	}
	Nterms++;
	TLSContrib **newRelLS = new TLSContrib*[Nterms];
	for (long i = 0; i < Nterms - 1; i++) {
		newRelLS[i] = RelLS[i];
	}
	newRelLS[Nterms - 1] = C;
	RelLS = newRelLS;
	return 0;
}

long TLSNonRel::Print() {
	cout << " [ " << GnrPrefac.FracStringSqrt() << "  G_" << L << S << " ] ";
	for (long i = 0; i < Nterms; i++) {
		RelLS[i]->PrintNR();
	}
	cout << endl;
	return 0;
}

long TLSNonRel::PrintG() {
	cout << " [ G_" << L << S << " ] ";
	for (long i = 0; i < Nterms; i++) {
		TFracNum GnrInv(GnrPrefac);
		GnrInv.Invert();
		RelLS[i]->PrintNRG(GnrInv);
	}
	cout << endl;
	return 0;
}
