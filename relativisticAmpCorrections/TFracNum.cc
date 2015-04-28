#include <iostream>
#include <string>
#include <stdio.h>
#include <string.h>
#include "TFracNum.h"
using namespace std;

bool TFracNum::debugFracNum = 0;

const char* SQUAREROOT_CHAR = "#";

TFracNum::TFracNum(long mN, long mD, long* N, long* D, long s)
	: _maxPrimNom(mN),
	  _maxPrimDen(mD),
	  _NOM(N),
	  _DEN(D),
	  _signPrefac(s)
{
	if (debugFracNum) {
		cout << "N:";
		for (long i = 0; i < _maxPrimNom; i++) {
			cout << _NOM[i] << ",";
		}
		cout << endl;
		cout << "D:";
		for (long i = 0; i < _maxPrimDen; i++) {
			cout << _DEN[i] << ",";
		}
		cout << endl;
		cout << "NOM pointer: " << _NOM;
		cout << "| ";
		for (long i = 0; i < mN; i++) {
			cout << _NOM[i] << ",";
		}
		cout << endl;
		cout << "DEN pointer: " << _DEN;
		cout << "| ";
		for (long i = 0; i < mD; i++) {
			cout << _DEN[i] << ",";
		}
		cout << endl;
	}
	SetINTs();
}

TFracNum::TFracNum(const long& N, const long& D, const string& s)
	: _maxPrimNom(0),
	  _maxPrimDen(0),
	  _NOM(0),
	  _DEN(0),
	  _signPrefac(1),
	  _NOM_INT(0),
	  _DEN_INT(0),
	  _dValue(0.)
{
	if (s == "factorial") {
		if (debugFracNum) {
			cout << s << endl;
		}
		if (N == D) {
			return;
		}
		long Low = N;
		long High = D;
		if (N > D) {
			Low = D;
			High = N;
		}
		long prim_vec[NPRIMFIELD];
		for (long i = 0; i < NPRIMFIELD; i++) {
			prim_vec[i] = 0;
		}
		long maxPrim = 0;
		for (long fac = Low + 1; fac <= High; fac++) {
			long rest = fac;
			long fmax = 0;
			while (rest != 1 && fmax < NPRIMFIELD) {
				while (rest % PRIMES[fmax] == 0) {
					prim_vec[fmax]++;
					rest /= PRIMES[fmax];
				}
				fmax++;
			}
			if (fmax > maxPrim) {
				maxPrim = fmax;
			}
		}
		if (N < D) {
			_maxPrimDen = maxPrim;
			_DEN = new long[_maxPrimDen];
			for (long jp = 0; jp < _maxPrimDen; jp++) {
				_DEN[jp] = prim_vec[jp];
			}
		} else {
			_maxPrimNom = maxPrim;
			_NOM = new long[_maxPrimNom];
			for (long jp = 0; jp < _maxPrimNom; jp++) {
				_NOM[jp] = prim_vec[jp];
			}
		}
	}
	this->SetINTs();
}

TFracNum::TFracNum(long inom, long iden) {

	if (debugFracNum) {
		cout << "Initializing with " << inom << "," << iden << endl;
	}
	_NOM_INT = inom;
	if (inom < 0) {
		_NOM_INT = -inom;
	}
	_DEN_INT = iden;
	if (iden < 0) {
		_DEN_INT = -iden;
	}

	if (inom == 0) {
		if (iden == 0) {
			_maxPrimNom = 0;
			_maxPrimDen = 0;
			_NOM = 0;
			_DEN = 0;
			_signPrefac = -6666;
			return;
		}
		_maxPrimNom = 0;
		_maxPrimDen = 0;
		_NOM = 0;
		_DEN = 0;
		_signPrefac = 0;
		_NOM_INT = 0;
		_DEN_INT = 1;
		_dValue = 0;
		return;
	}

	if (iden == 0) {
		_signPrefac = -7777;
		return;
	}

	_signPrefac = 1;
	if (inom < 0) {
		_signPrefac *= -1;
		inom *= -1;
	}
	if (iden < 0) {
		_signPrefac *= -1;
		iden *= -1;
	}

	if ((inom > 1000           or iden > 1000)             and // catch the const initialisations
	    (inom > MAXPRIMSQUARED or iden > MAXPRIMSQUARED))      // in the header file
	{
		if (inom > MAXPRIMSQUARED or iden > MAXPRIMSQUARED) {
			cerr << "MAXPRIMSQUARED reached!!! NOM=" << inom << ", DEN=" << iden << endl;
		}
		_maxPrimNom = -inom;
		_maxPrimDen = -iden;
		_NOM = 0;
		_DEN = 0;
		_NOM_INT = inom;
		_DEN_INT = iden;
		_dValue = _signPrefac * double(inom) / double(iden);
	} else {
		_maxPrimNom = 0;
		long rest_nom = inom;
		long prim_vec_nom[NPRIMFIELD];
		for (long i = 0; i < NPRIMFIELD; i++) {
			prim_vec_nom[i] = 0;
		}
		while (rest_nom != 1 && _maxPrimNom < NPRIMFIELD) {
			while (rest_nom % PRIMES[_maxPrimNom] == 0) {
				prim_vec_nom[_maxPrimNom]++;
				rest_nom /= PRIMES[_maxPrimNom];
			}
			_maxPrimNom++;
		}
		if (rest_nom != 1) {
			_maxPrimNom = -inom;
			_maxPrimDen = -iden;
			_NOM = 0;
			_DEN = 0;
		} else {
			_maxPrimDen = 0;
			long rest_den = iden;
			long prim_vec_den[NPRIMFIELD];
			for (long i = 0; i < NPRIMFIELD; i++) {
				prim_vec_den[i] = 0;
			}
			while (rest_den != 1 && _maxPrimDen < NPRIMFIELD) {
				while (rest_den % PRIMES[_maxPrimDen] == 0) {
					prim_vec_den[_maxPrimDen]++;
					rest_den /= PRIMES[_maxPrimDen];
				}
				_maxPrimDen++;
			}
			if (rest_den != 1) {
				_maxPrimNom = -inom;
				_maxPrimDen = -iden;
				_NOM = 0;
				_DEN = 0;
			} else {
				long maxPrim = _maxPrimNom;
				if (_maxPrimDen > maxPrim) {
					maxPrim = _maxPrimDen;
				}
				for (long ip = 0; ip < maxPrim; ip++) {
					if (prim_vec_nom[ip] != 0 && prim_vec_den[ip] != 0) {
						if (prim_vec_den[ip] > prim_vec_nom[ip]) {
							prim_vec_den[ip] -= prim_vec_nom[ip];
							prim_vec_nom[ip] = 0;
						} else {
							prim_vec_nom[ip] -= prim_vec_den[ip];
							prim_vec_den[ip] = 0;
						}
					}
				}

				_maxPrimNom = 0;
				_maxPrimDen = 0;
				for (long ip = 0; ip < NPRIMFIELD; ip++) {
					if (prim_vec_nom[ip] != 0) {
						_maxPrimNom = ip + 1;
					}
					if (prim_vec_den[ip] != 0) {
						_maxPrimDen = ip + 1;
					}
				}

				if (_maxPrimNom) {
					_NOM = new long[_maxPrimNom];
					for (long jp = 0; jp < _maxPrimNom; jp++) {
						_NOM[jp] = prim_vec_nom[jp];
					}
				} else {
					_NOM = 0;
				}
				if (_maxPrimDen) {
					_DEN = new long[_maxPrimDen];
					for (long jp = 0; jp < _maxPrimDen; jp++) {
						_DEN[jp] = prim_vec_den[jp];
					}
				} else {
					_DEN = 0;
				}
			}
		}
		this->SetINTs();
	}
}

bool TFracNum::SetINTs() {
	if (_signPrefac == 0) {
		_NOM_INT = 0;
		_DEN_INT = 1;
		_dValue = 0;
	} else {
		_NOM_INT = 1;
		_DEN_INT = 1;
		if (_maxPrimNom < 0 && _maxPrimDen < 0) {
			_NOM_INT = -_maxPrimNom;
			_DEN_INT = -_maxPrimDen;
		} else {

			long ip = _maxPrimNom - 1;
			while (ip >= 0 && _NOM[ip] == 0) {
				_maxPrimNom = ip;
				ip--;
			}
			//if (maxPrimNom==0) NOM=0;

			ip = _maxPrimDen - 1;
			while (ip >= 0 && _DEN[ip] == 0) {
				_maxPrimDen = ip;
				ip--;
			}

			long ipn = 0;
			while (ipn < _maxPrimNom) {
				for (long jj = 0; jj < _NOM[ipn]; jj++) {
					_NOM_INT *= PRIMES[ipn];
				}
				ipn++;
			}
			long ipd = 0;
			while (ipd < _maxPrimDen) {
				for (long jj = 0; jj < _DEN[ipd]; jj++) {
					_DEN_INT *= PRIMES[ipd];
				}
				ipd++;
			}
		}
		_dValue = double(_signPrefac) * double(_NOM_INT) / double(_DEN_INT);
	}
	return true;
}

long TFracNum::DenomCommonDivisor(const TFracNum &b) const {
	long maxPD = _maxPrimDen;
	if (maxPD > b._maxPrimDen) {
		maxPD = b._maxPrimDen;
	}
	long comdiv = 1;
	for (long i = 0; i < maxPD; i++) {
		long ppot = _DEN[i];
		if (b._DEN[i] < ppot) {
			ppot = b._DEN[i];
		}
		while (ppot-- > 0) {
			comdiv *= PRIMES[i];
		}
	}
	return comdiv;
}

TFracNum*
TFracNum::SumSignedRoots(TFracNum *b) {
	TFracNum mixed = (*this) * (*b);
	TFracNum aa = *this;
	TFracNum bb = *b;
	TFracNum *res = new TFracNum();
	if (mixed.Sqrt()) {
		bool flipsign = (aa.Dval() + bb.Dval() < 0);
		aa.Abs();
		bb.Abs();
		*res = aa + bb + TFracNum_Two * mixed;
		if (flipsign) {
			res->FlipSign();
		}
		return res;
	}
	cerr << "Error in TFracNum::SumSignedRoots()" << endl << "this:" << Dval() << endl;
	cerr << *this;
	cerr << "b:" << b->Dval() << endl;
	cerr << b;
	return 0;
}

bool TFracNum::Sqrt() {
	if (_signPrefac == 0 or _NOM_INT == 0) {
		return true;
	}
	if (debugFracNum) {
		long sqrt_ok = 1;
		for (long i = 0; i < _maxPrimNom; i++) {
			if (_NOM[i] % 2) {
				sqrt_ok = 0;
				break;
			}
		}
		if (sqrt_ok == 1) {
			for (long i = 0; i < _maxPrimDen; i++) {
				if (_DEN[i] % 2) {
					sqrt_ok = 0;
					break;
				}
			}
		}
		if (sqrt_ok == 0) {
			cout << "square root not possible for this fracnum :(" << endl;
			cout << *this;
		}
	}
	for (long i = 0; i < _maxPrimNom; i++) {
		if (_NOM[i] % 2) {
			return false;
		}
	}
	for (long i = 0; i < _maxPrimDen; i++) {
		if (_DEN[i] % 2) {
			return false;
		}
	}
	for (long i = 0; i < _maxPrimNom; i++) {
		_NOM[i] /= 2;
	}
	for (long i = 0; i < _maxPrimDen; i++) {
		_DEN[i] /= 2;
	}
	SetINTs();
	return true;
}
;

bool TFracNum::FlipSign() {
	_signPrefac *= -1;
	_dValue *= -1.0;
	return true;
}

bool TFracNum::Abs() {
	if (_signPrefac == 0) {
		return true;
	}
	_signPrefac = 1;
	if (_NOM_INT < 0) {
		_NOM_INT *= -1;
	}
	if (_dValue < 0) {
		_dValue *= -1.0;
	}
	return true;
}

bool TFracNum::Invert() {
	if (_signPrefac == -7777) {
		_maxPrimNom = 0;
		_maxPrimDen = 0;
		_NOM = 0;
		_DEN = 0;
		_signPrefac = -6666;
		SetINTs();
		return false;
	}
	if (_NOM_INT == 0) {
		_maxPrimNom = 0;
		_maxPrimDen = 0;
		_NOM = 0;
		_DEN = 0;
		_signPrefac = -7777;
		SetINTs();
		return false;
	}
	long MPN = _maxPrimNom;
	_maxPrimNom = _maxPrimDen;
	_maxPrimDen = MPN;
	long* NOMPTR = _NOM;
	_NOM = _DEN;
	_DEN = NOMPTR;
	SetINTs();
	return true;
}


bool TFracNum::operator==(const TFracNum &b) const {
	if (_signPrefac == 0 && b._signPrefac == 0) {
		return true;
	}
	if (_signPrefac != b._signPrefac) {
		return false;
	}
	if (_maxPrimNom != b._maxPrimNom) {
		return false;
	}
	if (_maxPrimDen != b._maxPrimDen) {
		return false;
	}
	for (long i = 0; i < _maxPrimNom; i++) {
		if (_NOM[i] != b._NOM[i]) {
			return false;
		}
	}
	for (long i = 0; i < _maxPrimDen; i++) {
		if (_DEN[i] != b._DEN[i]) {
			return false;
		}
	}
	return true;
}
;

bool TFracNum::PrintDifference(const TFracNum &b) const {
	if (_signPrefac == 0 && b._signPrefac == 0) {
		cout << "Both zero, they are equal." << endl;
		return true;
	}
	if (_signPrefac != b._signPrefac) {
		cout << "Different sign: " << _signPrefac << "!=" << b._signPrefac
				<< endl;
		return false;
	}
	if (_maxPrimNom != b._maxPrimNom) {
		cout << "Different maxPrimNom: " << _maxPrimNom << "!=" << b._maxPrimNom
				<< endl;
		return false;
	}
	if (_maxPrimDen != b._maxPrimDen) {
		cout << "Different maxPrimDen: " << _maxPrimDen << "!=" << b._maxPrimDen
				<< endl;
		return false;
	}
	for (long i = 0; i < _maxPrimNom; i++) {
		if (_NOM[i] != b._NOM[i]) {
			cout << "Different numerator contribution at prime " << i << ": "
					<< _NOM[i] << "!=" << b._NOM[i] << endl;
			return false;
		}
	}
	for (long i = 0; i < _maxPrimDen; i++) {
		if (_DEN[i] != b._DEN[i]) {
			cout << "Different denominator contribution at prime " << i << ": "
					<< _DEN[i] << "!=" << b._DEN[i] << endl;
			return false;
		}
	}
	cout << "Well, they're simply equal!" << endl;
	return true;
}
;

char*
TFracNum::HeaderString() {
	char* hstr = new char[30];
	if (_signPrefac == 0) {
		sprintf(hstr, "{0,1}");
		return hstr;
	}
	if (_signPrefac == -7777) {
		sprintf(hstr, "{1,0}");
		return hstr;
	}
	if (_signPrefac == -6666) {
		sprintf(hstr, "{0,0}");
		return hstr;
	}
	if (_signPrefac == 1) {
		sprintf(hstr, "{%ld,%ld}", _NOM_INT, _DEN_INT);
	} else {
		sprintf(hstr, "{%ld,%ld}", -_NOM_INT, _DEN_INT);
	}
	return hstr;
}
;

bool TFracNum::operator>(const TFracNum &b) const {
	if (_dValue > b._dValue) {
		return true;
	}
	return false;
}
;

TFracNum TFracNum::operator+(const TFracNum &b) const {
	long den_cdiv = DenomCommonDivisor(b);
	long bdc = b._DEN_INT / den_cdiv;
	long adc = _DEN_INT / den_cdiv;
	return TFracNum(_signPrefac * _NOM_INT * bdc + b._signPrefac * b._NOM_INT * adc, _DEN_INT * bdc);
}

TFracNum TFracNum::operator*(const TFracNum &b) const {
	// if one of the two numbers is undetermined,
	// the product is also undetermined
	if (_signPrefac == -6666 || b._signPrefac == -6666) {
		return TFracNum(0, 0, 0, 0, -6666);
	}

	// if one of the two numbers contains division by zero,
	// and the other nominator is zero, the product is undetermined
	if ((_signPrefac == -7777 and b._signPrefac == 0) or
	    (_signPrefac == 0     and b._signPrefac == -7777))
	{
		return TFracNum(0, 0, 0, 0, -6666);
	}

	// other cases with division by zero; product is also infinity
	if ((_signPrefac == -7777 or b._signPrefac == -7777)) {
		return TFracNum(0, 0, 0, 0, -7777);
	}

	if (_signPrefac * b._signPrefac == 0) {
		return TFracNum(0, 0, 0, 0, 0);
	}

	long maxPrimNom_ = _maxPrimNom;
	if (b._maxPrimNom > maxPrimNom_) {
		maxPrimNom_ = b._maxPrimNom;
	}
	long maxPrimDen_ = _maxPrimDen;
	if (b._maxPrimDen > maxPrimDen_) {
		maxPrimDen_ = b._maxPrimDen;
	}
	long maxPrim = maxPrimNom_;
	if (maxPrimDen_ > maxPrim) {
		maxPrim = maxPrimDen_;
	}

	long prim_vec_nom[maxPrim];
	long prim_vec_den[maxPrim];

	for (long ip = 0; ip < maxPrim; ip++) {
		prim_vec_nom[ip] = 0;
		prim_vec_den[ip] = 0;
		if (_maxPrimNom > ip) {
			prim_vec_nom[ip] += _NOM[ip];
		}
		if (b._maxPrimNom > ip) {
			prim_vec_nom[ip] += b._NOM[ip];
		}
		if (_maxPrimDen > ip) {
			prim_vec_den[ip] += _DEN[ip];
		}
		if (b._maxPrimDen > ip) {
			prim_vec_den[ip] += b._DEN[ip];
		}
	}

	for (long ip = 0; ip < maxPrim; ip++) {
		if (prim_vec_nom[ip] != 0 && prim_vec_den[ip] != 0) {
			if (prim_vec_den[ip] > prim_vec_nom[ip]) {
				prim_vec_den[ip] -= prim_vec_nom[ip];
				prim_vec_nom[ip] = 0;
			} else {
				prim_vec_nom[ip] -= prim_vec_den[ip];
				prim_vec_den[ip] = 0;
			}
		}
	}

	for (long ip = 0; ip < maxPrim; ip++) {
		if (prim_vec_nom[ip] != 0) {
			maxPrimNom_ = ip + 1;
		}
		if (prim_vec_den[ip] != 0) {
			maxPrimDen_ = ip + 1;
		}
	}

	long* NOM_ = new long[maxPrimNom_];
	for (long jp = 0; jp < maxPrimNom_; jp++) {
		NOM_[jp] = prim_vec_nom[jp];
	}

	long* DEN_ = new long[maxPrimDen_];
	for (long jp = 0; jp < maxPrimDen_; jp++) {
		DEN_[jp] = prim_vec_den[jp];
	}

	if (debugFracNum) {
		cout << " Initial with maxN=" << maxPrimNom_ << ", maxD=" << maxPrimDen_
		     << ", " << NOM_ << ", " << DEN_ << ", "
		     << _signPrefac * b._signPrefac
		     << endl
		     << "NOM:";
		for (long i = 0; i < maxPrimNom_; i++) {
			cout << NOM_[i] << ",";
		}
		cout << endl
		     << "DEN:";
		for (long i = 0; i < maxPrimDen_; i++) {
			cout << DEN_[i] << ",";
		}
		cout << endl;
	}
	return TFracNum(maxPrimNom_, maxPrimDen_, NOM_, DEN_, _signPrefac * b._signPrefac);
}

std::ostream& TFracNum::Print(std::ostream& out) const {
	if (debugFracNum) {
		out << "nom prime list: " << _maxPrimNom << ",pointer " << _NOM << endl;
		out << "den prime list: " << _maxPrimDen << ",pointer " << _DEN << endl;
		out << "NOM:";
		for (long i = 0; i < _maxPrimNom; i++) {
			out << _NOM[i] << ",";
		}
		out << endl;
		out << "DEN:";
		for (long i = 0; i < _maxPrimDen; i++) {
			out << _DEN[i] << ",";
		}
		out << endl;
	}
	out << "sign_prefac=" << _signPrefac << endl;
	if (_maxPrimNom < 0) {
		if (_signPrefac < 0) {
			out << "-";
		}
		out << -_maxPrimNom << "/" << -_maxPrimDen << endl;
		return out;
	}
	long integrity = 1;
	for (long i = 0; i < _maxPrimNom; i++) {
		if (_NOM[i] < 0 || _NOM[i] > 1000) {
			integrity = 0;
		}
	}
	for (long i = 0; i < _maxPrimDen; i++) {
		if (_DEN[i] < 0 || _DEN[i] > 1000) {
			integrity = 0;
		}
	}
	if (integrity == 0) {
		return out;
	}

	long nom = 1;
	long ipn = 0;
	if (_signPrefac < 0) {
		out << "-NOM = ";
	} else {
		out << " NOM = ";
	}
	long FirstTerm = 1;
	while (ipn < _maxPrimNom) {
		if (_NOM[ipn] != 0) {
			out << PRIMES[ipn] << "^" << _NOM[ipn];
			FirstTerm = 0;
		}
		for (long jj = 0; jj < _NOM[ipn]; jj++) {
			nom *= PRIMES[ipn];
		}
		ipn++;
		if (!FirstTerm && ipn < _maxPrimNom && _NOM[ipn] != 0) {
			out << " * ";
		}
	}
	out << " = " << nom << endl;

	long den = 1;
	long ipd = 0;
	out << " DEN = ";
	FirstTerm = 1;
	while (ipd < _maxPrimDen) {
		if (_DEN[ipd] != 0) {
			out << PRIMES[ipd] << "^" << _DEN[ipd];
			FirstTerm = 0;
		}
		for (long jj = 0; jj < _DEN[ipd]; jj++)
			den *= PRIMES[ipd];
		ipd++;
		if (!FirstTerm && ipd < _maxPrimDen && _DEN[ipd] != 0) {
			out << " * ";
		}
	}
	out << " = " << den << endl;
	out << "NOM_INT=" << _NOM_INT << endl;
	out << "DEN_INT=" << _DEN_INT << endl;
	out << "dvalue=" << _dValue << endl;
	return out;
}


const char*
TFracNum::FracString() {
	char *formstr = new char[50];
	char *fstr = new char[100];
	if (_NOM_INT == 0) {
		sprintf(fstr, "0");
	} else if (_DEN_INT == 1) {
		sprintf(formstr, "%%c%s", IOUTSTRING);
		sprintf(fstr, formstr, _signPrefac < 0 ? '-' : '+', _NOM_INT);
	} else {
		sprintf(formstr, "%%c%s/%s", IOUTSTRING, IOUTSTRING);
		sprintf(fstr, formstr, _signPrefac < 0 ? '-' : '+', _NOM_INT, _DEN_INT);
	}
	return fstr;
}
;

const char NULLSTRING[1] = ""; // workaround CINT warning when sprintf(s,"");

const char*
TFracNum::FracStringSqrt() const {
	char *formstr = new char[50];
	char *fstr = new char[200];
	if (_NOM_INT == 0) {
		sprintf(fstr, "0");
		return fstr;
	}
	long ipn = 0;
	long SQRT_NOM_INT = 1;
	long NOM_INT_REST = 1;
	while (ipn < _maxPrimNom) {
		for (long jj = 0; jj < _NOM[ipn] / 2; jj++) {
			SQRT_NOM_INT *= PRIMES[ipn];
		}
		if (_NOM[ipn] % 2) {
			NOM_INT_REST *= PRIMES[ipn];
		}
		ipn++;
	}
	long ipd = 0;
	long SQRT_DEN_INT = 1;
	long DEN_INT_REST = 1;
	while (ipd < _maxPrimDen) {
		for (long jj = 0; jj < _DEN[ipd] / 2; jj++) {
			SQRT_DEN_INT *= PRIMES[ipd];
		}
		if (_DEN[ipd] % 2) {
			DEN_INT_REST *= PRIMES[ipd];
		}
		ipd++;
	}

	char *sqrtstr = new char[100];
	bool one1 = false;
	bool one2 = false;
	if (SQRT_DEN_INT == 1) {
		if (SQRT_NOM_INT == 1) {
			sprintf(sqrtstr, "%s", NULLSTRING);
			one1 = true;
		} else {
			sprintf(formstr, "%s", IOUTSTRING);
			sprintf(sqrtstr, formstr, SQRT_NOM_INT);
		}
	} else {
		sprintf(formstr, "%s/%s", IOUTSTRING, IOUTSTRING);
		sprintf(sqrtstr, formstr, SQRT_NOM_INT, SQRT_DEN_INT);
	}

	char *reststr = new char[100];
	if (DEN_INT_REST == 1) {
		if (NOM_INT_REST == 1) {
			sprintf(reststr, "%s", NULLSTRING);
			one2 = true;
		} else {
			sprintf(formstr, "%s%s", SQUAREROOT_CHAR, IOUTSTRING);
			sprintf(reststr, formstr, NOM_INT_REST);
		}
	} else {
		sprintf(formstr, "%s%s/%s", SQUAREROOT_CHAR, IOUTSTRING, IOUTSTRING);
		sprintf(reststr, formstr, NOM_INT_REST, DEN_INT_REST);
	}

	if (one1 && one2) {
		sprintf(sqrtstr, "1");
	}
	sprintf(fstr, "%c%s%s", _signPrefac < 0 ? '-' : '+', sqrtstr, reststr);
	return fstr;
}
;

TFracNum a_to_J(long J, long m) {
	long kappa = (J - m) % 2;
	cout << "kappa=" << kappa << endl;
	long nom_ptr[1] = { 1 };
	TFracNum twofac(kappa, 0, nom_ptr, 0, 1);
	TFracNum fac1(J + m, 2 * J, "factorial");
	TFracNum fac2(J - m, 1, "factorial");
	return twofac * fac1 * fac2;
}

TFracNum am0_to_J(long J, long m, long m0) {
	long nom_ptr[1] = { m0 };
	TFracNum twofac(1, 0, nom_ptr, 0, 1);
	TFracNum fac1(J + m, 2 * J, "factorial");
	TFracNum fac2(J - m, 1, "factorial");
	return twofac * fac1 * fac2;
}

TFracNum c_sub_ell(long ell) {
	if (ell == 0) {
		return TFracNum(0, 0, 0, 0, 1);
	}
	long nom_ptr[1] = { ell };
	TFracNum two_to_ell(1, 0, nom_ptr, 0, 1);
	TFracNum fac1(ell, 1, "factorial");
	TFracNum fac2(ell, 2 * ell, "factorial");
	return two_to_ell * fac1 * fac2;
}

TFracNum cm0_sub_ell(long ell, long m0) {
	if (ell == 0) {
		return TFracNum(0, 0, 0, 0, 1);
	}
	long nom_ptr[1] = { (ell + m0) / 2 };
	TFracNum two_to_ell(1, 0, nom_ptr, 0, 1);
	TFracNum fac1(ell, 1, "factorial");
	TFracNum fac2(ell, 2 * ell, "factorial");
	return two_to_ell * fac1 * fac2;
}

TFracNum cm0_sub_ell_2(long ell, long m0) {
	//return  am0_to_J(ell, 0, m0);
	if (ell == 0) {
		return TFracNum(0, 0, 0, 0, 1);
	}
	long nom_ptr[1] = { (ell + m0) };
	TFracNum two_to_ell(1, 0, nom_ptr, 0, 1);
	TFracNum fac1a(ell, 1, "factorial");
	TFracNum fac1b(ell, 1, "factorial");
	TFracNum fac2a(ell, 2 * ell, "factorial");
	TFracNum fac2b(ell, 2 * ell, "factorial");
	return two_to_ell * fac1a * fac1b * fac2a * fac2b;
}

