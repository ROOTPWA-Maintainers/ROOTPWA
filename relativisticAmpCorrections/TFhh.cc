#include "TFhh.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "ClebschGordanBox.h"

#include <reportingUtils.hpp>

using namespace std;


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
	  _LSt(),
	  _NRLSt()
{
	long delta = _lambda - _nu;
	{
		stringstream sstr;
		sstr << "F_" << _lambda << "_" << _nu;
		_name_str = sstr.str();

	}

	for (size_t iLS = 0; iLS < LSampl.size(); iLS++) {
		if ( (LSampl[iLS]->Getdelta() == delta) or (LSampl[iLS]->Getdelta() == -delta) ) {
			// cout << "iLS=" << iLS << ", delta=" << delta << endl;
			TFracNum* CG3S = ClebschGordanBox::instance()->GetCG(LSampl[iLS]->GetS(), S1, S2);
			TFracNum SpinCouplFac = CG3S[ClebschGordanBox::CGIndex(S1, _lambda, S2, -_nu)];
			if (SpinCouplFac == TFracNum::Zero) {
				if (_debugLevel == 2) {
					cout << "Clebsch-Gordan is zero" << endl;
				}
			} else {
				_LSt.push_back(new TLSContrib(LSampl[iLS], delta, SpinCouplFac));
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
	  _LSt(sFhh->GetLSt()),
	  _NRLSt()
{
	if (flag != 'i' && flag != 'm') {
		printErr << "TFhh::TFhh unknown flag " << flag << endl;
		throw;
	}
	if (_debugLevel) {
		cout << "Initializing from single Amplitude" << endl;
	}
	if ( ( (flag == 'i') and ((sFhh->GetJ()) % 2))      or
	     ( (flag == 'm') and ((sFhh->GetJ()) % 2 == 0)))
	{
		cout << sFhh->GetName() << "[symm] = 0" << endl;
		_LSt = vector<TLSContrib*>();
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
	  _LSt(),
	  _NRLSt()
{
	{
		stringstream sstr;
		sstr << sFhh->GetName() << "[symm]";
		_name_str = sstr.str();
	}
	if (_J != xFhh->GetJ() or _evenContraction != xFhh->GetEvenContraction()
			or _lambda != xFhh->GetNu() or _nu != xFhh->GetLambda()) {
		printErr << "TFhh::TFhh(TFhh *, TFhh*): Something is wrong," << endl
		         << " source amplitudes for symmetrization do not match" << endl;
		throw;
	}

	// Since some LS-terms cancel out, first a max-length array of
	// LSContrib pointers is filled, and then squeezed to only
	// non-zero contributions

	const long& Ns = sFhh->GetNterms();
	const long& Nx = xFhh->GetNterms();
	vector<TLSContrib*> pLSt;
	for (long i = 0; i < (Ns+Nx); i++) {
		TLSContrib* toBeAdded;
		if (i < Ns) {
			toBeAdded = sFhh->GetLSt()[i];
		} else {
			toBeAdded = xFhh->GetLSt()[i - Ns];
		}
		bool foundInPrevterm = false;
		for (size_t j = 0; j < pLSt.size(); j++) {
			if (pLSt[j]->SameParameter(toBeAdded)) {
				foundInPrevterm = true;
				if (i < Ns) {
					pLSt[j]->Add(*toBeAdded, false);
				} else {
					pLSt[j]->Add(*toBeAdded, true);
				}
			}
		}
		if (not foundInPrevterm) {
			if (i < Ns) {
				pLSt.push_back(new TLSContrib(toBeAdded, false));
			} else {
				pLSt.push_back(new TLSContrib(toBeAdded, true));
			}
		}
	}

	for(size_t i = 0; i < pLSt.size(); ++i) {
		if(pLSt[i]->GetNterms() != 0) {
			_LSt.push_back(new TLSContrib(pLSt[i], false));
		}
	}

	Print();
}

void TFhh::NonRelLimit() {
	for (size_t i = 0; i < _LSt.size(); i++) {
		if (not _LSt[i]->IsPureRelativistic()) {
			long jfound = -1;
			for (size_t j = 0; j < _NRLSt.size(); j++) {
				if (_NRLSt[j]->CheckJLS(_LSt[i])) {
					jfound = j;
					break;
				}
			}
			if (jfound != -1) {
				_NRLSt[jfound]->Add(_LSt[i]);
			} else {
				_NRLSt.push_back(new TLSNonRel(_LSt[i]));
			}
		}
	}
	cout << _name_str << " (NR) = " << endl;
	for (size_t j = 0; j < _NRLSt.size(); j++) {
		cout << "LS=" << _NRLSt[j]->GetL() << _NRLSt[j]->GetS() << ": ";
		_NRLSt[j]->Print();
	}
}

void TFhh::PrintNRG() const
{
	cout << _name_str << " (NR) = " << endl;
	for (size_t j = 0; j < _NRLSt.size(); j++) {
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
	for (size_t iLSt = 0; iLSt < _LSt.size(); iLSt++) {
		_LSt[iLSt]->Print();
	}
	if (GetNterms() == 0) {
		cout << " 0" << endl;
	}
}
