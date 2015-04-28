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
	  _Nterms(b->_Nterms),
	  _NormFactor(b->_NormFactor),
	  _termFracNum(new TFracNum[_Nterms]),
	  _termg1pot(new long[_Nterms]),
	  _termg2pot(new long[_Nterms]),
	  _pureRelativistic(b->_pureRelativistic)
{
	for (long i = 0; i < _Nterms; i++) {
		// force new power fields by multiplication " *1 "
		_termFracNum[i] = TFracNum_One * b->_termFracNum[i];
		if (particleExchange) {
			_termg1pot[i] = b->_termg2pot[i];
			_termg2pot[i] = b->_termg1pot[i];
		} else {
			_termg1pot[i] = b->_termg1pot[i];
			_termg2pot[i] = b->_termg2pot[i];
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
	  _Nterms(A->GetNterms()),
	  _termFracNum(new TFracNum[_Nterms]),
	  _termg1pot(new long[_Nterms]),
	  _termg2pot(new long[_Nterms])
{
	bool sign_reversal = false;
	if (_delta < 0 && (_L + _S - _J) % 2) {
		sign_reversal = true;
	}

	_Nterms = A->GetNterms();
	_termFracNum = new TFracNum[_Nterms];
	_termg1pot = new long[_Nterms];
	_termg2pot = new long[_Nterms];

	_NormFactor = TFracNum_Zero;
	if (_debug) {
		cout << "(" << _J << ")" << _L << _S << "[" << _delta << "]" << endl;
	}
	for (long i = 0; i < _Nterms; i++) {
		TTensorTerm *tt = A->GetTerm(i); //TSScalar->GetTerm(i);

		if (_debug) {
			cout << scfac.FracStringSqrt() << ","
					<< tt->GetPreFac().FracStringSqrt() << " L" << _L << "S"
					<< _S << "J" << _J << "Ampldelta" << A->Getdelta()
					<< " delta" << _delta << ", sigrev " << sign_reversal;
		}

		_termFracNum[i] = scfac * tt->GetPreFac();
		if (sign_reversal == true) {
			_termFracNum[i].FlipSign();
		}
		//TFracNum tSqr = NormFactor * termFracNum[i];
		//if ( ! tSqrt.Sqrt() )
		//  cout << " Building LSContrib leads not to squared-fractional numbers,"
		//	   << " results will be wrong!" << endl;

		TFracNum *sum = _NormFactor.SumSignedRoots(&(_termFracNum[i]));
		if (sum) {
			_NormFactor = *sum;
		} else {
			cerr << "TLSContrib: Normalisation not root-fractional,"
					<< " *** results will be wrong *** " << endl;
		}

		//    NormFactor=NormFactor+termFracNum[i]+TFracNum_Two*tSqrt;

		_termg1pot[i] = tt->GetGamS();
		_termg2pot[i] = tt->GetGamSig();
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
		for (long i = 0; i < _Nterms; i++) {
			_termFracNum[i] = _termFracNum[i] * NormInv; // InvAbs;
		}
	}
}

long TLSContrib::Add(TLSContrib *b, bool particle_exchange) {
	if (_J != b->_J || _L != b->_L || _S != b->_S) {
		cerr << "TLSContrib::Add : Something is wrong, trying to add different"
		     << " (J;L,S): (" << _J << ";" << _L << "," << _S << ") != ("
		     << b->_J << ";" << b->_L << "," << b->_S << ")" << endl;
		return -1;
	}

	//
	// Include normalisation factor and the factor (1/2)**2 in the squared
	// representation of prefactors
	//

	for (long i = 0; i < _Nterms; i++) {
		_termFracNum[i] = TFracNum_Quarter * _NormFactor * _termFracNum[i];
	}

	for (long ib = 0; ib < b->_Nterms; ib++) {
		bool termSummed = false;
		for (long i = 0; i < _Nterms; i++) {
			if (!termSummed && _cNum == b->_cNum   and
			    ((particle_exchange == true         and
			    _termg1pot[i] == b->_termg2pot[ib]  and
			    _termg2pot[i] == b->_termg1pot[ib]) or
			    (particle_exchange == false         and
			    _termg1pot[i] == b->_termg1pot[ib]  and
			    _termg2pot[i] == b->_termg2pot[ib])))
			{
				termSummed = true;
				TFracNum bterm = TFracNum_Quarter * b->_NormFactor * b->_termFracNum[ib];

				if (_J % 2) {
					bterm.FlipSign();
				}

				TFracNum* sum = bterm.SumSignedRoots(&(_termFracNum[i]));
				if (sum) {
					_termFracNum[i] = *sum;
				} else {
					cerr << "TLSContrib: Normalisation not root-fractional,"
					     << " *** results will be wrong *** " << endl;
				}
			}
		}
		if (!termSummed) {
			_Nterms++;
			TFracNum* new_termFracNum = new TFracNum[_Nterms];
			long* new_termg1pot = new long[_Nterms];
			long* new_termg2pot = new long[_Nterms];
			for (long i = 0; i < _Nterms - 1; i++) {
				new_termFracNum[i] = _termFracNum[i];
				new_termg1pot[i] = _termg1pot[i];
				new_termg2pot[i] = _termg2pot[i];
			}
			new_termFracNum[_Nterms - 1] = TFracNum_Quarter * b->_NormFactor * b->_termFracNum[ib];
			//if ( ! new_termFracNum[Nterms-1].Sqrt() )
			//  cout << "Square root not possible, this will lead to wrong results:"
			//       << new_termFracNum[Nterms-1].FracString() << endl;
			//new_termFracNum[Nterms-1]=
			//	TFracNum_Half * new_termFracNum[Nterms-1] * b->NormFactor;
			if (_J % 2) {
				new_termFracNum[_Nterms - 1].FlipSign();
			}
			if (particle_exchange) {
				new_termg1pot[_Nterms - 1] = b->_termg2pot[ib];
				new_termg2pot[_Nterms - 1] = b->_termg1pot[ib];
			} else {
				new_termg1pot[_Nterms - 1] = b->_termg1pot[ib];
				new_termg2pot[_Nterms - 1] = b->_termg2pot[ib];
			}
			_termFracNum = new_termFracNum;
			_termg1pot = new_termg1pot;
			_termg2pot = new_termg2pot;
		}
	}

	//
	// Eliminate zero entries
	//
	long non_zeros = 0;
	for (long i = 0; i < _Nterms; i++) {
		if (!(_termFracNum[i] == TFracNum_Zero)) {
			non_zeros++;
		}
	}

	if (non_zeros == 0) {
		_Nterms = 0;
		return 0;
	} else {
		TFracNum *new_termFracNum = new TFracNum[non_zeros];
		long *new_termg1pot = new long[non_zeros];
		long *new_termg2pot = new long[non_zeros];

		long j = 0;
		for (long i = 0; i < _Nterms; i++) {
			if (!(_termFracNum[i] == TFracNum_Zero)) {
				new_termFracNum[j] = _termFracNum[i];
				new_termg1pot[j] = _termg1pot[i];
				new_termg2pot[j] = _termg2pot[i];
				j++;
			}
		}
		_Nterms = non_zeros;
		_termFracNum = new_termFracNum;
		_termg1pot = new_termg1pot;
		_termg2pot = new_termg2pot;
	}

	//
	// Recalculate Normalization Factor
	//
	_NormFactor = TFracNum_Zero;
	for (long i = 0; i < _Nterms; i++) {
		TFracNum *sum = _NormFactor.SumSignedRoots(&(_termFracNum[i]));
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
		for (long i = 0; i < _Nterms; i++) {
			_termFracNum[i] = _termFracNum[i] * NormInv;
		}
	}
	return _Nterms;
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
	for (long iT = 0; iT < _Nterms; iT++) {
		cout << _termFracNum[iT].FracStringSqrt() << " ";
		if (_termg1pot[iT]) {
			if (_termg1pot[iT] == 1) {
				cout << " gs";
			} else {
				cout << " gs^" << _termg1pot[iT];
			}
		}
		if (_termg2pot[iT]) {
			if (_termg2pot[iT] == 1) {
				cout << " gsig";
			} else {
				cout << " gsig^" << _termg2pot[iT];
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
