#include "TLSContrib.h"

#include <iostream>

#include <reportingUtils.hpp>

using namespace std;


bool TLSContrib::_debug = false;


TLSContrib::TLSContrib(const TLSContrib* b, const bool& particleExchange)
	: _J(b->_J),
	  _L(b->_L),
	  _S(b->_S),
	  _cNum(b->_cNum),
	  _delta(b->_delta),
	  _NormFactor(b->_NormFactor),
	  _polynomialTerms(b->getPolynomialTerms()),
	  _pureRelativistic(b->_pureRelativistic)
{
	for (size_t i = 0; i < _polynomialTerms.size(); i++) {
		if (particleExchange) {
			_polynomialTerms[i].swapExponents();
		}
	}
}


TLSContrib::TLSContrib(const TLSAmpl* A,
                       const long& delta,
                       const TFracNum& scfac)
	: _J(A->GetJ()),
	  _L(A->GetL()),
	  _S(A->GetS()),
	  _cNum(A->GetContraction()),
	  _delta(delta),
	  _SpinCG(scfac),
	  _polynomialTerms(A->GetNterms())
{

	const bool signReversal = (_delta < 0 and (_L + _S - _J) % 2);

	_NormFactor = TFracNum::Zero;
	if (_debug) {
		cout << "(" << _J << ")" << _L << _S << "[" << _delta << "]" << endl;
	}
	for (size_t i = 0; i < _polynomialTerms.size(); i++) {
		const TTensorTerm& tt = A->GetTerm(i); //TSScalar->GetTerm(i);

		if (_debug) {
			cout << scfac.FracStringSqrt() << ","
			     << tt.GetPreFac().FracStringSqrt() << " L" << _L << "S"
			     << _S << "J" << _J << "Ampldelta" << A->Getdelta()
			     << " delta" << _delta << ", sigrev " << signReversal;
		}


		_polynomialTerms[i].squareOfPrefactor = scfac * tt.GetPreFac();
		if (signReversal == true) {
			_polynomialTerms[i].squareOfPrefactor.FlipSign();
		}

		//TFracNum tSqr = NormFactor * termFracNum[i];
		//if ( ! tSqrt.Sqrt() )
		//  cout << " Building LSContrib leads not to squared-fractional numbers,"
		//	   << " results will be wrong!" << endl;

		_NormFactor = _NormFactor.SumSignedRoots(_polynomialTerms[i].squareOfPrefactor);

		//    NormFactor=NormFactor+termFracNum[i]+TFracNum_Two*tSqrt;

		_polynomialTerms[i].exponentOfGammaS     = tt.GetGamS();
		_polynomialTerms[i].exponentOfGammaSigma = tt.GetGamSig();
		if (_debug) {
			cout << " -> Normfactor: " << _NormFactor.FracStringSqrt() << endl;
		}
	}
	if (_NormFactor == TFracNum::Zero) {
		_pureRelativistic = true;
		_NormFactor = TFracNum::One;
	} else {
		_pureRelativistic = false;
		TFracNum NormInv = _NormFactor;
		NormInv.Invert();
		// TFracNum InvAbs=TFracNum_One * NormInv; //Bug: real copy forced by "1 *"
		// InvAbs.Abs();
		for (size_t i = 0; i < GetNterms(); i++) {
			_polynomialTerms[i].squareOfPrefactor *= NormInv; // InvAbs
		}
	}
}


void TLSContrib::Add(const TLSContrib& rhs, bool particleExchange) {
	if (_J != rhs._J or _L != rhs._L or _S != rhs._S) {
		printErr << "TLSContrib::Add : Something is wrong, trying to add different"
		         << " (J;L,S): (" << _J << ";" << _L << "," << _S << ") != ("
		         << rhs._J << ";" << rhs._L << "," << rhs._S << ")" << endl;
		throw;
	}

	//
	// Include normalisation factor and the factor (1/2)**2 in the squared
	// representation of prefactors
	//

	for (size_t i = 0; i < _polynomialTerms.size(); i++) {
		_polynomialTerms[i].squareOfPrefactor *= TFracNum::Quarter * _NormFactor;
	}

	for (size_t ib = 0; ib < rhs.GetNterms(); ib++) {
		bool termSummed = false;
		for (size_t i = 0; i < _polynomialTerms.size(); i++) {
			if (not termSummed and _cNum == rhs._cNum                                            and
			    (
			      (particleExchange                                                              and
			       _polynomialTerms[i].exponentOfGammaS     == rhs.getPolynomialTerms()[ib].exponentOfGammaSigma and
			       _polynomialTerms[i].exponentOfGammaSigma == rhs.getPolynomialTerms()[ib].exponentOfGammaS
			      )                                                                              or
			      (
			       not particleExchange                                                          and
			       _polynomialTerms[i].exponentOfGammaS     == rhs.getPolynomialTerms()[ib].exponentOfGammaS     and
			       _polynomialTerms[i].exponentOfGammaSigma == rhs.getPolynomialTerms()[ib].exponentOfGammaSigma
			      )
			    ) )
			{
				termSummed = true;
				TFracNum bterm = TFracNum::Quarter * rhs._NormFactor * rhs.getPolynomialTerms()[ib].squareOfPrefactor;

				if (_J % 2) {
					bterm.FlipSign();
				}

				_polynomialTerms[i].squareOfPrefactor = bterm.SumSignedRoots(_polynomialTerms[i].squareOfPrefactor);

			}
		}
		if (not termSummed) {

			polynomialTerms factor(rhs.getPolynomialTerms()[ib]);
			factor.squareOfPrefactor *= TFracNum::Quarter * rhs._NormFactor;
			if (_J % 2) {
				if(not factor.squareOfPrefactor.FlipSign()) {
					printErr << "Flipping sign failed. Aborting..." << endl;
					throw;
				}
			}
			if(particleExchange) {
				factor.swapExponents();
			}
			_polynomialTerms.push_back(factor);
		}
	}

	//
	// Eliminate zero entries
	//
	for (size_t i = _polynomialTerms.size(); i > 0; --i) {
		if (_polynomialTerms[i-1].squareOfPrefactor == TFracNum::Zero) {
			_polynomialTerms.erase(_polynomialTerms.begin() + (i-1));
		}
	}

	//
	// Recalculate Normalization Factor
	//
	_NormFactor = TFracNum::Zero;
	for (size_t i = 0; i < _polynomialTerms.size(); i++) {
		_NormFactor = _NormFactor.SumSignedRoots(_polynomialTerms[i].squareOfPrefactor);
	}

	//
	// Apply normalization
	//
	if (_NormFactor == TFracNum::Zero) {
		_pureRelativistic = true;
		_NormFactor = TFracNum::One;
	} else {
		_pureRelativistic = false;
		TFracNum NormInv = _NormFactor;
		NormInv.Invert();
		for (size_t i = 0; i < _polynomialTerms.size(); i++) {
			_polynomialTerms[i].squareOfPrefactor *= NormInv;
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
	if (not _pureRelativistic) {
		cout << _NormFactor.FracStringSqrt() << " ) ( ";
	}
	if(_polynomialTerms.empty()) {
		cout << "0";
	}
	for (size_t iT = 0; iT < _polynomialTerms.size(); iT++) {
		const polynomialTerms& factor = _polynomialTerms[iT];
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
