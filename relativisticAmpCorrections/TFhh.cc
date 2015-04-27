#include "TFhh.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

extern ClebschGordanBox box;


unsigned int TFhh::_debugLevel = 1;


TFhh::TFhh(const long& J,
           const long& S1,
           const long& S2,
           const long& lambda,
           const long& nu,
           const vector<TLSAmpl*>& LSampl,
           const bool& evenContraction)
	: _J(J),
	  _lambda(lambda),
	  _nu(nu),
	  _evenContraction(evenContraction),
	  _nTerms(0)
{
	long delta = _lambda - _nu;
	{
		stringstream sstr;
		sstr << "F_" << _lambda << "_" << _nu;
		_name_str = sstr.str();

	}

	for (unsigned int iLS = 0; iLS < LSampl.size(); iLS++) {
		if ( (LSampl[iLS]->Getdelta() == delta) or (LSampl[iLS]->Getdelta() == -delta) ) {
			// cout << "iLS=" << iLS << ", delta=" << delta << endl;
			TFracNum* CG3S = box.GetCG(LSampl[iLS]->GetS(), S1, S2);
			TFracNum SpinCouplFac = CG3S[CGIndex(S1, _lambda, S2, -_nu)];
			if (SpinCouplFac == TFracNum_Zero) {
				if (_debugLevel == 2) {
					cout << "Clebsch-Gordan is zero" << endl;
				}
			} else {
				TLSContrib* *newPTR = new TLSContrib*[_nTerms + 1];
				for (long iC = 0; iC < _nTerms; iC++) {
					newPTR[iC] = _LSt[iC];
				}
				newPTR[_nTerms] = new TLSContrib(LSampl[iLS], delta, SpinCouplFac);
				// delete[] LSt;
				_nTerms++;
				_LSt = newPTR;
			}
		}
	}
	if (_debugLevel) {
		Print();
	}
}

TFhh::TFhh(TFhh *sFhh, char flag)
	: _J(sFhh->GetJ()),
	  _lambda(sFhh->GetLambda()),
	  _nu(sFhh->GetNu()),
	  _evenContraction(sFhh->GetEvenContraction()),
	  _nTerms(sFhh->GetNterms()),
	  _LSt(sFhh->GetLStPtr())
{
	if (flag != 'i' && flag != 'm') {
		cerr << "TFhh::TFhh unknown flag " << flag << endl;
		return;
	}
	if (_debugLevel) {
		cout << "Initializing from single Amplitude" << endl;
	}
	if ( ( (flag == 'i') and ((sFhh->GetJ()) % 2)) or ((flag == 'm') and ((sFhh->GetJ()) % 2 == 0)) ) {
		cout << sFhh->GetName() << "[symm] = 0" << endl;
		_nTerms = 0;
	} else {
		{
			stringstream sstr;
			sstr << sFhh->GetName() << "[symm]";
			_name_str = sstr.str();
		}
		if (_debugLevel) {
			Print();
		}
	}
}

TFhh::TFhh(TFhh *sFhh, TFhh *xFhh)
	: _J(sFhh->GetJ()),
	  _lambda(sFhh->GetLambda()),
	  _nu(sFhh->GetNu()),
	  _evenContraction(sFhh->GetEvenContraction()),
	  _nTerms(0),
	  _LSt(0)
{
	{
		stringstream sstr;
		sstr << sFhh->GetName() << "[symm]";
		_name_str = sstr.str();
	}
	if (_J != xFhh->GetJ() or _evenContraction != xFhh->GetEvenContraction()
			or _lambda != xFhh->GetNu() or _nu != xFhh->GetLambda()) {
		cerr << "TFhh::TFhh(TFhh *, TFhh*): Something is wrong," << endl
				<< " source amplitudes for symmetrization do not match" << endl;
		throw;
	}

	// Since some LS-terms cancel out, first a max-length array of
	// LSContrib pointers is filled, and then squeezed to only
	// non-zero contributions

	const long& Ns = sFhh->GetNterms();
	const long& Nx = xFhh->GetNterms();
	_nTerms = Ns + Nx;
	TLSContrib* *pLSt = new TLSContrib*[_nTerms];

	long prevTerms = 0;

	for (long i = 0; i < _nTerms; i++) {
		TLSContrib* toBeAdded;
		if (i < Ns) {
			toBeAdded = sFhh->GetLStPtr()[i];
		} else {
			toBeAdded = xFhh->GetLStPtr()[i - Ns];
		}
		long foundInPrevterm = 0;
		for (long j = 0; j < prevTerms; j++) {
			if (pLSt[j]->SameParameter(toBeAdded)) {
				foundInPrevterm = 1;
				if (i < Ns) {
					pLSt[j]->Add(toBeAdded, false);
				} else {
					pLSt[j]->Add(toBeAdded, true);
				}
			}
		}
		if (not foundInPrevterm) {
			if (i < Ns) {
				pLSt[prevTerms] = new TLSContrib(toBeAdded, false);
			} else {
				pLSt[prevTerms] = new TLSContrib(toBeAdded, true);
			}
			prevTerms++;
		}
	}

	//
	// Cancel zeros
	//
	_nTerms = 0;
	for (long i = 0; i < prevTerms; i++) {
		if (pLSt[i]->GetNterms() != 0) {
			_nTerms++;
		}
	}

	if (_nTerms) {
		_LSt = new TLSContrib*[_nTerms];
		long j = 0;
		for (long i = 0; i < prevTerms; i++) {
			if (pLSt[i]->GetNterms() != 0) {
				_LSt[j] = new TLSContrib(pLSt[i], false);
				j++;
			}
		}
	}
	Print();
}

void TFhh::NonRelLimit() {
	_NNRterms = 0;
	for (long i = 0; i < _nTerms; i++) {
		if (not _LSt[i]->IsPureRelativistic()) {
			long jfound = -1;
			for (long j = 0; j < _NNRterms; j++) {
				if (_NRLSt[j]->CheckJLS(_LSt[i])) {
					jfound = j;
					break;
				}
			}
			if (jfound != -1) {
				_NRLSt[jfound]->Add(_LSt[i]);
			} else {
				_NNRterms++;
				TLSNonRel **newNRLS = new TLSNonRel*[_NNRterms];
				for (long j = 0; j < _NNRterms - 1; j++) {
					newNRLS[j] = _NRLSt[j];
				}
				newNRLS[_NNRterms - 1] = new TLSNonRel(_LSt[i]);
				_NRLSt = newNRLS;
			}
		}
	}
	cout << _name_str << " (NR) = " << endl;
	for (long j = 0; j < _NNRterms; j++) {
		cout << "LS=" << _NRLSt[j]->GetL() << _NRLSt[j]->GetS() << ": ";
		_NRLSt[j]->Print();
	}
}

void TFhh::PrintNRG() const
{
	cout << _name_str << " (NR) = " << endl;
	for (long j = 0; j < _NNRterms; j++) {
		cout << "LS=" << _NRLSt[j]->GetL() << _NRLSt[j]->GetS() << " => ";
		_NRLSt[j]->PrintG();
	}
}

void TFhh::Print() const
{
	cout << _name_str << " =";
	if (_evenContraction) {
		cout << endl;
	} else {
		cout << " (iw)" << endl;
	}
	for (long iLSt = 0; iLSt < _nTerms; iLSt++) {
		_LSt[iLSt]->Print();
	}
	if (_nTerms == 0) {
		cout << " 0" << endl;
	}
}
