#include "TFracNum.h"

#include <iostream>
#include <string>
#include <stdio.h>
#include <string.h>

#include <reportingUtils.hpp>

using namespace std;
using rpwa::operator<<;

bool TFracNum::debugFracNum = false;


const TFracNum TFracNum::Zero = TFracNum(0, 1);
const TFracNum TFracNum::One = TFracNum(1, 1);
const TFracNum TFracNum::Two = TFracNum(2, 1);
const TFracNum TFracNum::mTwo = TFracNum(-2, 1);
const TFracNum TFracNum::Quarter = TFracNum(1, 4);


const char* SQUAREROOT_CHAR = "#";


TFracNum::TFracNum(const vector<long>& N, const vector<long>& D, long s)
	: _NOM(N),
	  _DEN(D),
	  _signPrefac(s),
	  _numerator(0),
	  _nomCacheRebuildRequired(true),
	  _denominator(0),
	  _denCacheRebuildRequired(true),
	  _value(0.),
	  _valueCacheRebuildRequired(true)
{
	if (debugFracNum) {
		printDebug << "NOM: " << _NOM << endl;
		printDebug << "DEN: " << _DEN << endl;
	}
}


TFracNum::TFracNum(const long& N, const long& D, const string& s)
	: _NOM(),
	  _DEN(),
	  _signPrefac(1),
	  _numerator(0),
	  _nomCacheRebuildRequired(true),
	  _denominator(0),
	  _denCacheRebuildRequired(true),
	  _value(0.),
	  _valueCacheRebuildRequired(true)
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
		for (long fac = Low + 1; fac <= High; fac++) {
			size_t rest = fac;
			size_t fmax = 0;
			while (rest != 1 && fmax < NPRIMFIELD) {
				while (rest % PRIMES[fmax] == 0) {
					if(N < D) {
						if(fmax >= _DEN.size()) {
							_DEN.resize(fmax+1, 0);
						}
						_DEN[fmax]++;
					} else {
						if(fmax >= _NOM.size()) {
							_NOM.resize(fmax+1, 0);
						}
						_NOM[fmax]++;
					}
					rest /= PRIMES[fmax];
				}
				fmax++;
			}
		}
	}
}


TFracNum::TFracNum(long inom, long iden)
	: _NOM(),
	  _DEN(),
	  _signPrefac(0),
	  _numerator(0),
	  _nomCacheRebuildRequired(true),
	  _denominator(0),
	  _denCacheRebuildRequired(true),
	  _value(0.),
	  _valueCacheRebuildRequired(true)
{

	if (debugFracNum) {
		cout << "Initializing with " << inom << "," << iden << endl;
	}

	if (inom == 0) {
		if (iden == 0) {
			_signPrefac = -6666;
			return;
		}
		_value = 0.;
		_valueCacheRebuildRequired = false;
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
			//TODO: check if throw is correct
			printErr << "MAXPRIMSQUARED reached!!! NOM=" << inom << ", DEN=" << iden << endl;
			throw;
		}
	} else {
		long rest_nom = inom;
		size_t maxPrimNom = 0;
		while (rest_nom != 1 && maxPrimNom < NPRIMFIELD) {
			while (rest_nom % PRIMES[maxPrimNom] == 0) {
				if(maxPrimNom >= _NOM.size()) {
					_NOM.resize(maxPrimNom+1, 0);
				}
				_NOM[maxPrimNom]++;
				rest_nom /= PRIMES[maxPrimNom];
			}
			maxPrimNom++;
		}
		if (rest_nom != 1) {
			//TODO: check what this does
			printErr << "(" << inom << "/" << iden << ") overflow?? Aborting..." << endl;
			cout << endl << (*this) << endl;
			throw;
		} else {
			size_t maxPrimDen = 0;
			long rest_den = iden;
			while (rest_den != 1 && maxPrimDen < NPRIMFIELD) {
				while (rest_den % PRIMES[maxPrimDen] == 0) {
					if(maxPrimDen >= _DEN.size()) {
						_DEN.resize(maxPrimDen+1, 0);
					}
					_DEN[maxPrimDen]++;
					rest_den /= PRIMES[maxPrimDen];
				}
				maxPrimDen++;
			}
			if (rest_den != 1) {
				//TODO: check what this does
				printErr << "(" << inom << "/" << iden << ") overflow?? Aborting..." << endl;
				cout << endl << (*this) << endl;
				throw;
			} else {
				const size_t minPrim = min(_NOM.size(), _DEN.size());
				for (size_t ip = 0; ip < minPrim; ip++) {
					if (_NOM[ip] != 0 and _DEN[ip] != 0) {
						if (_DEN[ip] > _NOM[ip]) {
							_DEN[ip] -= _NOM[ip];
							_NOM[ip] = 0;
						} else {
							_NOM[ip] -= _DEN[ip];
							_DEN[ip] = 0;
						}
					}
				}
				TFracNum::removeZerosFromVector(_NOM);
				TFracNum::removeZerosFromVector(_DEN);
			}
		}
	}
}


const long& TFracNum::GetNumerator() const
{
	if(_nomCacheRebuildRequired) {
		if(_signPrefac == 0) {
			_numerator = 0;
		} else {
			_numerator = getNumberFromFactorization(_NOM);
		}
		_nomCacheRebuildRequired = false;
	}
	return _numerator;
}


const long& TFracNum::GetDenominator() const
{
	if(_denCacheRebuildRequired) {
		if(_signPrefac == 0) {
			_denominator = 1;
		} else {
			_denominator = getNumberFromFactorization(_DEN);
		}
		_denCacheRebuildRequired = false;
	}
	return _denominator;
}


const double& TFracNum::Dval() const
{
	if(_valueCacheRebuildRequired) {
		if(_signPrefac == 0) {
			_value = 0.;
		} else {
			_value = ((double)_signPrefac) * ((double)GetNumerator()) / ((double)GetDenominator());
		}
		_valueCacheRebuildRequired = false;
	}
	return _value;
}


void TFracNum::resetAllCaches() const
{
	_numerator   = 0;
	_denominator = 0;
	_value       = 0.;
	_nomCacheRebuildRequired   = true;
	_denCacheRebuildRequired   = true;
	_valueCacheRebuildRequired = true;
}


long TFracNum::DenomCommonDivisor(const TFracNum &b) const {
	size_t minPD = min(_DEN.size(), b._DEN.size());
	long comdiv = 1;
	for (size_t i = 0; i < minPD; i++) {
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


TFracNum
TFracNum::SumSignedRoots(const TFracNum& b) const
{
	TFracNum mixed = (*this) * b;
	TFracNum aa = *this;
	TFracNum bb = b;
	if (mixed.Sqrt()) {
		bool flipsign = (aa.Dval() + bb.Dval() < 0);
		aa.Abs();
		bb.Abs();
		TFracNum res = aa + bb + TFracNum::Two * mixed;
		if (flipsign) {
			res.FlipSign();
		}
		return res;
	}
	printErr << "Error in TFracNum::SumSignedRoots()" << endl << "this:" << Dval() << endl
	         << *this << endl
	         << "b:" << b.Dval() << endl
	         << b;
	throw;
}


bool TFracNum::Sqrt() {
	if (_signPrefac == 0 or GetNumerator() == 0) {
		return true;
	}
	// TODO: move this to the loops further down (or remove it completely?)
	if (debugFracNum) {
		long sqrt_ok = 1;
		for (size_t i = 0; i < _NOM.size(); i++) {
			if (_NOM[i] % 2) {
				sqrt_ok = 0;
				break;
			}
		}
		if (sqrt_ok == 1) {
			for (size_t i = 0; i < _DEN.size(); i++) {
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
	for (size_t i = 0; i < _NOM.size(); i++) {
		if (_NOM[i] % 2) {
			return false;
		}
	}
	for (size_t i = 0; i < _DEN.size(); i++) {
		if (_DEN[i] % 2) {
			return false;
		}
	}
	for (size_t i = 0; i < _NOM.size(); i++) {
		_NOM[i] /= 2;
	}
	for (size_t i = 0; i < _DEN.size(); i++) {
		_DEN[i] /= 2;
	}
	resetAllCaches();
	return true;
}


bool TFracNum::FlipSign() {
	_signPrefac *= -1;
	if(not _valueCacheRebuildRequired) {
		_value      *= -1.;
	}
	return true;
}


bool TFracNum::Abs() {
	if (_signPrefac == 0) {
		return true;
	}
	_signPrefac = 1;
	if(not _nomCacheRebuildRequired) {
		_numerator =  abs(_numerator);
	}
	if(not _valueCacheRebuildRequired) {
		_value     = fabs(_value);
	}
	return true;
}


bool TFracNum::Invert() {
	if (_signPrefac == -7777) {
		_NOM = vector<long>();
		_DEN = vector<long>();
		_signPrefac = -6666;
		resetAllCaches();
		return false;
	}
	if (GetNumerator() == 0) {
		_NOM = vector<long>();
		_DEN = vector<long>();
		_signPrefac = -7777;
		resetAllCaches();
		return false;
	}
	vector<long> oldNOM = _NOM;
	_NOM = _DEN;
	_DEN = oldNOM;
	resetAllCaches();
	return true;
}


bool TFracNum::operator==(const TFracNum &b) const {
	if (_signPrefac == 0 && b._signPrefac == 0) {
		return true;
	}
	if (_signPrefac != b._signPrefac) {
		return false;
	}
	if (_NOM.size() != b._NOM.size()) {
		return false;
	}
	if (_DEN.size() != b._DEN.size()) {
		return false;
	}
	for (size_t i = 0; i < _NOM.size(); i++) {
		if (_NOM[i] != b._NOM[i]) {
			return false;
		}
	}
	for (size_t i = 0; i < _DEN.size(); i++) {
		if (_DEN[i] != b._DEN[i]) {
			return false;
		}
	}
	return true;
}


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
	if (_NOM.size() != b._NOM.size()) {
		cout << "Different maxPrimNom: " << _NOM.size() << "!=" << b._NOM.size()
		     << endl;
		return false;
	}
	if (_DEN.size() != b._DEN.size()) {
		cout << "Different maxPrimDen: " << _DEN.size() << "!=" << b._DEN.size()
		     << endl;
		return false;
	}
	for (size_t i = 0; i < _NOM.size(); i++) {
		if (_NOM[i] != b._NOM[i]) {
			cout << "Different numerator contribution at prime " << i << ": "
			     << _NOM[i] << "!=" << b._NOM[i] << endl;
			return false;
		}
	}
	for (size_t i = 0; i < _DEN.size(); i++) {
		if (_DEN[i] != b._DEN[i]) {
			cout << "Different denominator contribution at prime " << i << ": "
			     << _DEN[i] << "!=" << b._DEN[i] << endl;
			return false;
		}
	}
	cout << "Well, they're simply equal!" << endl;
	return true;
}


const char*
TFracNum::HeaderString() const
{
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
		sprintf(hstr, "{%ld,%ld}", GetNumerator(), GetDenominator());
	} else {
		sprintf(hstr, "{%ld,%ld}", -GetNumerator(), GetDenominator());
	}
	return hstr;
}


bool TFracNum::operator>(const TFracNum& rhs) const {
	if (Dval() > rhs.Dval()) {
		return true;
	}
	return false;
}



TFracNum& TFracNum::operator+=(const TFracNum &rhs) {
	long den_cdiv = DenomCommonDivisor(rhs);
	long bdc = rhs.GetDenominator() / den_cdiv;
	long adc = GetDenominator() / den_cdiv;
	*this = TFracNum(_signPrefac * GetNumerator() * bdc + rhs._signPrefac * rhs.GetNumerator() * adc, GetDenominator() * bdc);
	return *this;
}


TFracNum& TFracNum::operator*=(const TFracNum& rhs) {

	// if one of the two numbers is undetermined,
	// the product is also undetermined
	if (_signPrefac == -6666 or rhs._signPrefac == -6666) {
		*this = TFracNum(vector<long>(), vector<long>(), -6666);
		return *this;
	}

	// if one of the two numbers contains division by zero,
	// and the other nominator is zero, the product is undetermined
	if ((_signPrefac == -7777 and rhs._signPrefac == 0) or
	    (_signPrefac == 0     and rhs._signPrefac == -7777))
	{
		*this = TFracNum(vector<long>(), vector<long>(), -6666);
		return *this;
	}

	// other cases with division by zero; product is also infinity
	if ((_signPrefac == -7777 or rhs._signPrefac == -7777)) {
		*this = TFracNum(vector<long>(), vector<long>(), -7777);
		return *this;
	}

	if (_signPrefac * rhs._signPrefac == 0) {
		*this = TFracNum::Zero;
		return *this;
	}

	if(_NOM.size() < rhs._NOM.size()) {
		_NOM.resize(rhs._NOM.size(), 0);
	}
	for(size_t i = 0; i < rhs._NOM.size(); ++i) {
		_NOM[i] += rhs._NOM[i];
	}

	if(_DEN.size() < rhs._DEN.size()) {
		_DEN.resize(rhs._DEN.size(), 0);
	}
	for(size_t i = 0; i < rhs._DEN.size(); ++i) {
		_DEN[i] += rhs._DEN[i];
	}

	for(size_t i = 0; i < min(_NOM.size(), _DEN.size()); ++i) {
		if (_NOM[i] != 0 and _DEN[i] != 0) {
			if (_DEN[i] > _NOM[i]) {
				_DEN[i] -= _NOM[i];
				_NOM[i] = 0;
			} else {
				_NOM[i] -= _DEN[i];
				_DEN[i] = 0;
			}
		}
	}

	_signPrefac *= rhs._signPrefac;
	resetAllCaches();
	return *this;
}


std::ostream& TFracNum::Print(std::ostream& out) const {
	if (debugFracNum) {
		out << "nom prime list: " << _NOM.size() << ",pointer " << _NOM << endl;
		out << "den prime list: " << _DEN.size() << ",pointer " << _DEN << endl;
		out << "NOM:";
		for (size_t i = 0; i < _NOM.size(); i++) {
			out << _NOM[i] << ",";
		}
		out << endl;
		out << "DEN:";
		for (size_t i = 0; i < _DEN.size(); i++) {
			out << _DEN[i] << ",";
		}
		out << endl;
	}
	out << "sign_prefac=" << _signPrefac << endl;

	bool integrity = true;
	for (size_t i = 0; i < _NOM.size(); i++) {
		if (_NOM[i] < 0 or _NOM[i] > 1000) {
			integrity = false;
		}
	}
	for (size_t i = 0; i < _DEN.size(); i++) {
		if (_DEN[i] < 0 or _DEN[i] > 1000) {
			integrity = false;
		}
	}
	if (not integrity) {
		return out;
	}

	long nom = 1;
	size_t ipn = 0;
	if (_signPrefac < 0) {
		out << "-NOM = ";
	} else {
		out << " NOM = ";
	}
	long FirstTerm = 1;
	while (ipn < _NOM.size()) {
		if (_NOM[ipn] != 0) {
			out << PRIMES[ipn] << "^" << _NOM[ipn];
			FirstTerm = 0;
		}
		for (long jj = 0; jj < _NOM[ipn]; jj++) {
			nom *= PRIMES[ipn];
		}
		ipn++;
		if (!FirstTerm and ipn < _NOM.size() and _NOM[ipn] != 0) {
			out << " * ";
		}
	}
	out << " = " << nom << endl;

	long den = 1;
	size_t ipd = 0;
	out << " DEN = ";
	FirstTerm = 1;
	while (ipd < _DEN.size()) {
		if (_DEN[ipd] != 0) {
			out << PRIMES[ipd] << "^" << _DEN[ipd];
			FirstTerm = 0;
		}
		for (long jj = 0; jj < _DEN[ipd]; jj++)
			den *= PRIMES[ipd];
		ipd++;
		if (!FirstTerm and ipd < _DEN.size() and _DEN[ipd] != 0) {
			out << " * ";
		}
	}
	out << " = " << den << endl;
	out << "NOM_INT=" << GetNumerator() << endl;
	out << "DEN_INT=" << GetDenominator() << endl;
	out << "dvalue=" << Dval() << endl;
	return out;
}


//TODO: fix this horrible sprintf mess
const char*
TFracNum::FracString() const
{
	char *formstr = new char[50];
	char *fstr = new char[100];
	const long& numerator = GetNumerator();
	const long& denominator = GetDenominator();
	if (numerator == 0) {
		sprintf(fstr, "0");
	} else if (denominator == 1) {
		sprintf(formstr, "%%c%s", IOUTSTRING);
		sprintf(fstr, formstr, _signPrefac < 0 ? '-' : '+', numerator);
	} else {
		sprintf(formstr, "%%c%s/%s", IOUTSTRING, IOUTSTRING);
		sprintf(fstr, formstr, _signPrefac < 0 ? '-' : '+', numerator, denominator);
	}
	return fstr;
}


const char NULLSTRING[1] = ""; // workaround CINT warning when sprintf(s,"");


//TODO: fix this horrible sprintf mess
const char*
TFracNum::FracStringSqrt() const {
	char *formstr = new char[50];
	char *fstr = new char[200];
	if (GetNumerator() == 0) {
		sprintf(fstr, "0");
		return fstr;
	}
	size_t ipn = 0;
	long SQRT_NOM_INT = 1;
	long NOM_INT_REST = 1;
	while (ipn < _NOM.size()) {
		for (long jj = 0; jj < _NOM[ipn] / 2; jj++) {
			SQRT_NOM_INT *= PRIMES[ipn];
		}
		if (_NOM[ipn] % 2) {
			NOM_INT_REST *= PRIMES[ipn];
		}
		ipn++;
	}
	size_t ipd = 0;
	long SQRT_DEN_INT = 1;
	long DEN_INT_REST = 1;
	while (ipd < _DEN.size()) {
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


void TFracNum::removeZerosFromVector(vector<long>& vector)
{
	if(vector.empty()) {
		return;
	}
	size_t neededSize = vector.size() - 1;
	for( ; neededSize >= 0; --neededSize) {
		if(vector[neededSize] != 0) {
			++neededSize;
			break;
		}
	}
	vector.resize(neededSize);
}


long TFracNum::getNumberFromFactorization(const vector<long>& vector)
{
	long retval = 1;
	for(size_t i = 0; i < vector.size(); ++i) {
		for(long j = 0; j < vector[i]; ++j) {
			retval *= PRIMES[i];
		}
	}
	return retval;
}


// TODO: check if this can be deleted
#if(0)
TFracNum TFracNum::a_to_J(long J, long m) {
	long kappa = (J - m) % 2;
	cout << "kappa=" << kappa << endl;
	long nom_ptr[1] = { 1 };
	TFracNum twofac(kappa, 0, nom_ptr, 0, 1);
	TFracNum fac1(J + m, 2 * J, "factorial");
	TFracNum fac2(J - m, 1, "factorial");
	return twofac * fac1 * fac2;
}
#endif

TFracNum TFracNum::am0_to_J(long J, long m, long m0) {
	TFracNum twofac(vector<long>(1, m0),
	                vector<long>(),
	                1);
	TFracNum fac1(J + m, 2 * J, "factorial");
	TFracNum fac2(J - m, 1, "factorial");
	return twofac * fac1 * fac2;
}


TFracNum TFracNum::c_sub_ell(long ell) {
	if (ell == 0) {
		return TFracNum::One;
	}
	TFracNum two_to_ell(vector<long>(1, ell),
	                    vector<long>(),
	                    1);
	TFracNum fac1(ell, 1, "factorial");
	TFracNum fac2(ell, 2 * ell, "factorial");
	return two_to_ell * fac1 * fac2;
}


TFracNum TFracNum::cm0_sub_ell(long ell, long m0) {
	if (ell == 0) {
		return TFracNum::One;
	}
	TFracNum two_to_ell(vector<long>(1, (ell + m0) / 2),
	                    vector<long>(),
	                    1);
	TFracNum fac1(ell, 1, "factorial");
	TFracNum fac2(ell, 2 * ell, "factorial");
	return two_to_ell * fac1 * fac2;
}


TFracNum TFracNum::cm0_sub_ell_2(long ell, long m0) {
	//return  am0_to_J(ell, 0, m0);
	if (ell == 0) {
		return TFracNum::One;
	}
	TFracNum two_to_ell(vector<long>(1, (ell + m0) ),
	                    vector<long>(),
	                    1);
	TFracNum fac1a(ell, 1, "factorial");
	TFracNum fac1b(ell, 1, "factorial");
	TFracNum fac2a(ell, 2 * ell, "factorial");
	TFracNum fac2b(ell, 2 * ell, "factorial");
	return two_to_ell * fac1a * fac1b * fac2a * fac2b;
}
