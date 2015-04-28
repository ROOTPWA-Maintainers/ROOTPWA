#include "TLSContrib.h"

#include <iostream>

using namespace std;

bool TLSContrib::_debug = false;

TLSContrib::TLSContrib(const TLSContrib* b, const bool& particleExchange)
	: _J(b->_J),
	  _L(b->_L),
	  _S(b->_S),
	  _cNum(b->_cNum),
	  _delta(b->_delta),
	  _NormFactor(b->_NormFactor),
	  _factors(b->getFactors()),
	  _pureRelativistic(b->_pureRelativistic)
{
	for (size_t i = 0; i < _factors.size(); i++) {
		if (particleExchange) {
			_factors[i].swapExponents();
		}
	}
}
;

TLSContrib::TLSContrib(/*const*/ TLSAmpl* A,
                       const long& delta,
                       const TFracNum& scfac)
	: _J(A->GetJ()),
	  _L(A->GetL()),
	  _S(A->GetS()),
	  _cNum(A->GetContraction()),
	  _delta(delta),
	  _SpinCG(scfac),
	  _factors(A->GetNterms())
{

	const bool signReversal = (_delta < 0 and (_L + _S - _J) % 2);

	_NormFactor = TFracNum_Zero;
	if (_debug) {
		cout << "(" << _J << ")" << _L << _S << "[" << _delta << "]" << endl;
	}
	for (size_t i = 0; i < _factors.size(); i++) {
		TTensorTerm *tt = A->GetTerm(i); //TSScalar->GetTerm(i);

		if (_debug) {
			cout << scfac.FracStringSqrt() << ","
			     << tt->GetPreFac().FracStringSqrt() << " L" << _L << "S"
			     << _S << "J" << _J << "Ampldelta" << A->Getdelta()
			     << " delta" << _delta << ", sigrev " << signReversal;
		}


		_factors[i].squareOfPrefactor = scfac * tt->GetPreFac();
		if (signReversal == true) {
			_factors[i].squareOfPrefactor.FlipSign();
		}

		//TFracNum tSqr = NormFactor * termFracNum[i];
		//if ( ! tSqrt.Sqrt() )
		//  cout << " Building LSContrib leads not to squared-fractional numbers,"
		//	   << " results will be wrong!" << endl;

		TFracNum* sum = _NormFactor.SumSignedRoots(_factors[i].squareOfPrefactor);
		if (sum) {
			_NormFactor = *sum;
		} else {
			cerr << "TLSContrib: Normalisation not root-fractional,"
					<< " *** results will be wrong *** " << endl;
		}

		//    NormFactor=NormFactor+termFracNum[i]+TFracNum_Two*tSqrt;

		_factors[i].exponentOfGammaS     = tt->GetGamS();
		_factors[i].exponentOfGammaSigma = tt->GetGamSig();
		if (_debug) {
			cout << " -> Normfactor: " << _NormFactor.FracStringSqrt() << endl;
		}
	}
	if (_NormFactor == TFracNum_Zero) {
		_pureRelativistic = true;
		_NormFactor = TFracNum_One;
	} else {
		_pureRelativistic = false;
		TFracNum NormInv = _NormFactor;
		NormInv.Invert();
		// TFracNum InvAbs=TFracNum_One * NormInv; //Bug: real copy forced by "1 *"
		// InvAbs.Abs();
		for (size_t i = 0; i < GetNterms(); i++) {
			_factors[i].squareOfPrefactor = NormInv * _factors[i].squareOfPrefactor;
		}
	}
}

void TLSContrib::Add(TLSContrib *b, bool particleExchange) {
	if (_J != b->_J or _L != b->_L or _S != b->_S) {
		cerr << "TLSContrib::Add : Something is wrong, trying to add different"
		     << " (J;L,S): (" << _J << ";" << _L << "," << _S << ") != ("
		     << b->_J << ";" << b->_L << "," << b->_S << ")" << endl;
		throw;
	}

	//
	// Include normalisation factor and the factor (1/2)**2 in the squared
	// representation of prefactors
	//

	for (size_t i = 0; i < _factors.size(); i++) {
		_factors[i].squareOfPrefactor = TFracNum_Quarter * _NormFactor * _factors[i].squareOfPrefactor;
	}

	for (size_t ib = 0; ib < b->GetNterms(); ib++) {
		bool termSummed = false;
		for (size_t i = 0; i < _factors.size(); i++) {
			if (not termSummed and _cNum == b->_cNum                                            and
			    (
			      (particleExchange                                                             and
			       _factors[i].exponentOfGammaS     == b->getFactors()[ib].exponentOfGammaSigma and
			       _factors[i].exponentOfGammaSigma == b->getFactors()[ib].exponentOfGammaS
			      )                                                                             or
			      (
			       not particleExchange                                                         and
			       _factors[i].exponentOfGammaS     == b->getFactors()[ib].exponentOfGammaS     and
			       _factors[i].exponentOfGammaSigma == b->getFactors()[ib].exponentOfGammaSigma
			      )
			    ) )
			{
				termSummed = true;
				TFracNum bterm = TFracNum_Quarter * b->_NormFactor * b->getFactors()[ib].squareOfPrefactor;

				if (_J % 2) {
					bterm.FlipSign();
				}

				TFracNum* sum = bterm.SumSignedRoots(_factors[i].squareOfPrefactor);
				if (sum) {
					_factors[i].squareOfPrefactor = *sum;
				} else {
					cerr << "TLSContrib: Normalisation not root-fractional,"
					     << " *** results will be wrong *** " << endl;
				}
			}
		}
		if (not termSummed) {

			factors factor(b->getFactors()[ib]);
			factor.squareOfPrefactor = TFracNum_Quarter * b->_NormFactor * factor.squareOfPrefactor;
			if (_J % 2) {
				if(not factor.squareOfPrefactor.FlipSign()) {
					cerr << "Flipping sign failed. Aborting..." << endl;
					throw;
				}
			}
			if(particleExchange) {
				factor.swapExponents();
			}
			_factors.push_back(factor);
		}
	}

	//
	// Eliminate zero entries
	//
	vector<int> indicesToRemove;
	for (size_t i = _factors.size()-1; i >= 0; --i) {
		if (_factors[i].squareOfPrefactor == TFracNum_Zero) {
			indicesToRemove.push_back(i);
		}
	}

	if (indicesToRemove.empty()) {
		return;
	} else {
		for(size_t i = 0; i < indicesToRemove.size(); ++i) {
			_factors.erase(_factors.begin()+i);
		}
	}

	//
	// Recalculate Normalization Factor
	//
	_NormFactor = TFracNum_Zero;
	for (size_t i = 0; i < _factors.size(); i++) {
		TFracNum* sum = _NormFactor.SumSignedRoots(_factors[i].squareOfPrefactor);
		if (sum) {
			_NormFactor = *sum;
		} else {
			cerr << "TLSContrib: Normalisation not root-fractional,"
					<< " *** results will be wrong *** " << endl;
		}
	}

	//
	// Apply normalization
	//
	if (_NormFactor == TFracNum_Zero) {
		_pureRelativistic = true;
		_NormFactor = TFracNum_One;
	} else {
		_pureRelativistic = false;
		TFracNum NormInv = _NormFactor;
		NormInv.Invert();
		for (size_t i = 0; i < _factors.size(); i++) {
			_factors[i].squareOfPrefactor = NormInv * _factors[i].squareOfPrefactor;
		}
	}
}

void TLSContrib::Print() const
{
	if (_cNum == 1) {
		cout << "g" << "[" << _cNum << "] (";
	}
	if (_cNum == 2) {
		cout << "f" << "[" << _cNum << "] (";
	}
	if (_cNum >= 3) {
		cout << "h" << "[" << _cNum << "] (";
	}
	cout << _J << ")" << _L << _S << "( ";
	if (!_pureRelativistic) {
		cout << _NormFactor.FracStringSqrt() << " ) ( ";
	}
	for (size_t iT = 0; iT < _factors.size(); iT++) {
		const factors& factor = _factors[iT];
		cout << factor.squareOfPrefactor.FracStringSqrt() << " ";
		if (factor.exponentOfGammaS) {
			if (factor.exponentOfGammaS == 1) {
				cout << " gs";
			} else {
				cout << " gs^" << factor.exponentOfGammaS;
			}
		}
		if (factor.exponentOfGammaSigma) {
			if (factor.exponentOfGammaSigma == 1) {
				cout << " gsig";
			} else {
				cout << " gsig^" << factor.exponentOfGammaSigma;
			}
		}
	}
	cout << " )" << endl;
}


void TLSContrib::PrintNR() const
{
	cout << _NormFactor.FracStringSqrt();
	if (_cNum == 1) {
		cout << "*g";
	}
	if (_cNum == 2) {
		cout << "*f";
	}
	if (_cNum == 3) {
		cout << "*h";
	}
}


void TLSContrib::PrintNRG(TFracNum m) const
{
	cout << (_NormFactor * m).FracStringSqrt();
	if (_cNum == 1) {
		cout << "*g";
	}
	if (_cNum == 2) {
		cout << "*f";
	}
	if (_cNum == 3) {
		cout << "*h";
	}
}
