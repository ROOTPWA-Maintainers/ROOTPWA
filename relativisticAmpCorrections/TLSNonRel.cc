#include "TLSNonRel.h"

#include <iostream>

using namespace std;

TLSNonRel::TLSNonRel(TLSContrib* C)
	: _J(C->GetJ()),
	  _L(C->GetL()),
	  _S(C->GetS()),
	  _RelLS(1, C)
{
	TFracNum *JdL0Sd = ClebschGordanBox::instance()->GetCG(_J, _L, _S);
	//cout << "delta=" << C->GetDelta()
	//     << ",S=" << S
	//     << ", 2L+1/2J+1=" << TFracNum(2*L+1,2*J+1).FracStringSqrt()
	//     << ",CG(JdL0Sd)="
	//     << JdL0Sd[CGIndex(L, 0, S, C->GetDelta())].FracStringSqrt()
	//     << ", SpinCG=" << C->GetSpinCG()->FracStringSqrt()
	//     << endl;
	_GnrPrefac = TFracNum(2 * _L + 1, 2 * _J + 1)
	             * JdL0Sd[ClebschGordanBox::CGIndex(_L, 0, _S, C->GetDelta())]
	             * *C->GetSpinCG();
}

void TLSNonRel::Add(TLSContrib* C) {
	if (!CheckJLS(C)) {
		cout << "TLSNonRel::Add not appropriate. Aborting..."
				<< endl;
		throw;
	}
	_RelLS.push_back(C);
}

long TLSNonRel::Print() {
	cout << " [ " << _GnrPrefac.FracStringSqrt() << "  G_" << _L << _S << " ] ";
	for (size_t i = 0; i < _RelLS.size(); i++) {
		_RelLS[i]->PrintNR();
	}
	cout << endl;
	return 0;
}

long TLSNonRel::PrintG() {
	cout << " [ G_" << _L << _S << " ] ";
	for (size_t i = 0; i < _RelLS.size(); i++) {
		TFracNum GnrInv(_GnrPrefac);
		GnrInv.Invert();
		_RelLS[i]->PrintNRG(GnrInv);
	}
	cout << endl;
	return 0;
}
