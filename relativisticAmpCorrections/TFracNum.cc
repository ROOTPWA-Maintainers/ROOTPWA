#include <iostream>
#include <string>
#include <stdio.h>
#include <string.h>
#include "TFracNum.h"
using namespace std;

long debugFracNum = 0;
char factorial[] = { 'f', 'a', 'c', 't', 'o', 'r', 'i', 'a', 'l' };

const char* SQUAREROOT_CHAR = "#";

TFracNum::TFracNum(long mN, long mD, long* N, long* D, long s) {
	if (debugFracNum) {
		cout << "N:";
		for (long i = 0; i < mN; i++)
			cout << N[i] << ",";
		cout << endl;
		cout << "D:";
		for (long i = 0; i < mD; i++)
			cout << D[i] << ",";
		cout << endl;
	}
	maxPrimNom = mN;
	maxPrimDen = mD;
	NOM = N;
	DEN = D;
	sign_prefac = s;
	if (debugFracNum) {
		cout << "NOM pointer: " << NOM;
		cout << "| ";
		for (long i = 0; i < mN; i++)
			cout << NOM[i] << ",";
		cout << endl;
		cout << "DEN pointer: " << DEN;
		cout << "| ";
		for (long i = 0; i < mD; i++)
			cout << DEN[i] << ",";
		cout << endl;
	}
	this->SetINTs();
}

TFracNum::TFracNum(long N, long D, const char* s) {
	maxPrimNom = 0;
	maxPrimDen = 0;
	NOM = 0;
	DEN = 0;
	sign_prefac = 1;
	if (strcmp(s, "factorial") == 0) {
		if (debugFracNum)
			cout << s << endl;
		if (N == D)
			return;
		long Low = N;
		long High = D;
		if (N > D) {
			Low = D;
			High = N;
		}
		long prim_vec[NPRIMFIELD];
		for (long i = 0; i < NPRIMFIELD; i++)
			prim_vec[i] = 0;
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
			if (fmax > maxPrim)
				maxPrim = fmax;
		}
		if (N < D) {
			maxPrimDen = maxPrim;
			DEN = new long[maxPrimDen];
			for (long jp = 0; jp < maxPrimDen; jp++)
				DEN[jp] = prim_vec[jp];
		} else {
			maxPrimNom = maxPrim;
			NOM = new long[maxPrimNom];
			for (long jp = 0; jp < maxPrimNom; jp++)
				NOM[jp] = prim_vec[jp];
		}
	}
	this->SetINTs();
}

TFracNum::TFracNum(long inom, long iden) {

	if (debugFracNum) {
		cout << "Initializing with " << inom << "," << iden << endl;
	}
	NOM_INT = inom;
	if (inom < 0)
		NOM_INT = -inom;
	DEN_INT = iden;
	if (iden < 0)
		DEN_INT = -iden;

	if (inom == 0) {
		if (iden == 0) {
			maxPrimNom = 0;
			maxPrimDen = 0;
			NOM = 0;
			DEN = 0;
			sign_prefac = -6666;
			return;
		}
		maxPrimNom = 0;
		maxPrimDen = 0;
		NOM = 0;
		DEN = 0;
		sign_prefac = 0;
		NOM_INT = 0;
		DEN_INT = 1;
		dvalue = 0;
		return;
	}

	if (iden == 0) {
		sign_prefac = -7777;
		return;
	}

	sign_prefac = 1;
	if (inom < 0) {
		sign_prefac *= -1;
		inom *= -1;
	}
	if (iden < 0) {
		sign_prefac *= -1;
		iden *= -1;
	}

	if ((inom > 1000 || iden > 1000) &&    // catch the const initialisations
			(inom > MAXPRIMSQUARED ||      // in the header file
					iden > MAXPRIMSQUARED)) {

		if (inom > MAXPRIMSQUARED || iden > MAXPRIMSQUARED) {
			cerr << "MAXPRIMSQUARED reached!!! NOM=" << inom << ", DEN=" << iden
					<< endl;
		}
		maxPrimNom = -inom;
		maxPrimDen = -iden;
		NOM = 0;
		DEN = 0;
		NOM_INT = inom;
		DEN_INT = iden;
		dvalue = sign_prefac * double(inom) / double(iden);
	} else {
		maxPrimNom = 0;
		long rest_nom = inom;
		long prim_vec_nom[NPRIMFIELD];
		for (long i = 0; i < NPRIMFIELD; i++)
			prim_vec_nom[i] = 0;
		while (rest_nom != 1 && maxPrimNom < NPRIMFIELD) {
			while (rest_nom % PRIMES[maxPrimNom] == 0) {
				prim_vec_nom[maxPrimNom]++;
				rest_nom /= PRIMES[maxPrimNom];
			}
			maxPrimNom++;
		}
		if (rest_nom != 1) {
			maxPrimNom = -inom;
			maxPrimDen = -iden;
			NOM = 0;
			DEN = 0;
		} else {
			maxPrimDen = 0;
			long rest_den = iden;
			long prim_vec_den[NPRIMFIELD];
			for (long i = 0; i < NPRIMFIELD; i++)
				prim_vec_den[i] = 0;
			while (rest_den != 1 && maxPrimDen < NPRIMFIELD) {
				while (rest_den % PRIMES[maxPrimDen] == 0) {
					prim_vec_den[maxPrimDen]++;
					rest_den /= PRIMES[maxPrimDen];
				}
				maxPrimDen++;
			}
			if (rest_den != 1) {
				maxPrimNom = -inom;
				maxPrimDen = -iden;
				NOM = 0;
				DEN = 0;
			} else {
				long maxPrim = maxPrimNom;
				if (maxPrimDen > maxPrim)
					maxPrim = maxPrimDen;

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

				maxPrimNom = 0;
				maxPrimDen = 0;
				for (long ip = 0; ip < NPRIMFIELD; ip++) {
					if (prim_vec_nom[ip] != 0)
						maxPrimNom = ip + 1;
					if (prim_vec_den[ip] != 0)
						maxPrimDen = ip + 1;
				}

				if (maxPrimNom) {
					NOM = new long[maxPrimNom];
					for (long jp = 0; jp < maxPrimNom; jp++)
						NOM[jp] = prim_vec_nom[jp];
				} else {
					NOM = 0;
				}
				if (maxPrimDen) {
					DEN = new long[maxPrimDen];
					for (long jp = 0; jp < maxPrimDen; jp++)
						DEN[jp] = prim_vec_den[jp];
				} else {
					DEN = 0;
				}
			}
		}
		this->SetINTs();
	}
}

bool TFracNum::SetINTs() {
	if (sign_prefac == 0) {
		NOM_INT = 0;
		DEN_INT = 1;
		dvalue = 0;
	} else {
		NOM_INT = 1;
		DEN_INT = 1;
		if (maxPrimNom < 0 && maxPrimDen < 0) {
			NOM_INT = -maxPrimNom;
			DEN_INT = -maxPrimDen;
		} else {

			long ip = maxPrimNom - 1;
			while (ip >= 0 && NOM[ip] == 0) {
				maxPrimNom = ip;
				ip--;
			}
			//if (maxPrimNom==0) NOM=0;

			ip = maxPrimDen - 1;
			while (ip >= 0 && DEN[ip] == 0) {
				maxPrimDen = ip;
				ip--;
			}

			long ipn = 0;
			while (ipn < maxPrimNom) {
				for (long jj = 0; jj < NOM[ipn]; jj++)
					NOM_INT *= PRIMES[ipn];
				ipn++;
			}
			long ipd = 0;
			while (ipd < maxPrimDen) {
				for (long jj = 0; jj < DEN[ipd]; jj++)
					DEN_INT *= PRIMES[ipd];
				ipd++;
			}
		}
		dvalue = double(sign_prefac) * double(NOM_INT) / double(DEN_INT);
	}
	return true;
}

long TFracNum::DenomCommonDivisor(const TFracNum &b) const {
	long maxPD = maxPrimDen;
	if (maxPD > b.maxPrimDen)
		maxPD = b.maxPrimDen;
	long comdiv = 1;
	for (long i = 0; i < maxPD; i++) {
		long ppot = DEN[i];
		if (b.DEN[i] < ppot)
			ppot = b.DEN[i];
		while (ppot-- > 0)
			comdiv *= PRIMES[i];
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
		if (flipsign)
			res->FlipSign();
		return res;
	}
	cerr << "Error in TFracNum::SumSignedRoots()" << endl << "this:"
			<< this->Dval() << endl;
	this->PrintToErr();
	cerr << "b:" << b->Dval() << endl;
	b->PrintToErr();
	return 0;
}

bool TFracNum::Sqrt() {
	if (sign_prefac == 0 || NOM_INT == 0)
		return true;
	if (debugFracNum) {
		long sqrt_ok = 1;
		for (long i = 0; i < maxPrimNom; i++)
			if (NOM[i] % 2) {
				sqrt_ok = 0;
				break;
			}
		if (sqrt_ok == 1)
			for (long i = 0; i < maxPrimDen; i++)
				if (DEN[i] % 2) {
					sqrt_ok = 0;
					break;
				}
		if (sqrt_ok == 0) {
			cout << "square root not possible for this fracnum :(" << endl;
			this->Print();
		}
	}
	for (long i = 0; i < maxPrimNom; i++)
		if (NOM[i] % 2)
			return false;
	for (long i = 0; i < maxPrimDen; i++)
		if (DEN[i] % 2)
			return false;
	for (long i = 0; i < maxPrimNom; i++)
		NOM[i] /= 2;
	for (long i = 0; i < maxPrimDen; i++)
		DEN[i] /= 2;
	this->SetINTs();
	return true;
}
;

bool TFracNum::FlipSign() {
	sign_prefac *= -1;
	dvalue *= -1.0;
	return true;
}

bool TFracNum::Abs() {
	if (sign_prefac == 0)
		return true;
	sign_prefac = 1;
	if (NOM_INT < 0)
		NOM_INT *= -1;
	if (dvalue < 0)
		dvalue *= -1.0;
	return true;
}

bool TFracNum::Invert() {
	if (sign_prefac == -7777) {
		maxPrimNom = 0;
		maxPrimDen = 0;
		NOM = 0;
		DEN = 0;
		sign_prefac = -6666;
		this->SetINTs();
		return false;
	}
	if (NOM_INT == 0) {
		maxPrimNom = 0;
		maxPrimDen = 0;
		NOM = 0;
		DEN = 0;
		sign_prefac = -7777;
		this->SetINTs();
		return false;
	}
	long MPN = maxPrimNom;
	maxPrimNom = maxPrimDen;
	maxPrimDen = MPN;
	long* NOMPTR = NOM;
	NOM = DEN;
	DEN = NOMPTR;
	this->SetINTs();
	return true;
}

long TFracNum::GetSign() {
	if (dvalue < 0)
		return -1;
	else
		return 1;
}

bool TFracNum::operator==(const TFracNum &b) const {
	if (sign_prefac == 0 && b.sign_prefac == 0)
		return true;
	if (sign_prefac != b.sign_prefac)
		return false;
	if (maxPrimNom != b.maxPrimNom)
		return false;
	if (maxPrimDen != b.maxPrimDen)
		return false;
	for (long i = 0; i < maxPrimNom; i++)
		if (NOM[i] != b.NOM[i])
			return false;
	for (long i = 0; i < maxPrimDen; i++)
		if (DEN[i] != b.DEN[i])
			return false;
	return true;
}
;

bool TFracNum::PrintDifference(const TFracNum &b) const {
	if (sign_prefac == 0 && b.sign_prefac == 0) {
		cout << "Both zero, they are equal." << endl;
		return true;
	}
	if (sign_prefac != b.sign_prefac) {
		cout << "Different sign: " << sign_prefac << "!=" << b.sign_prefac
				<< endl;
		return false;
	}
	if (maxPrimNom != b.maxPrimNom) {
		cout << "Different maxPrimNom: " << maxPrimNom << "!=" << b.maxPrimNom
				<< endl;
		return false;
	}
	if (maxPrimDen != b.maxPrimDen) {
		cout << "Different maxPrimDen: " << maxPrimDen << "!=" << b.maxPrimDen
				<< endl;
		return false;
	}
	for (long i = 0; i < maxPrimNom; i++)
		if (NOM[i] != b.NOM[i]) {
			cout << "Different numerator contribution at prime " << i << ": "
					<< NOM[i] << "!=" << b.NOM[i] << endl;
			return false;
		}
	for (long i = 0; i < maxPrimDen; i++)
		if (DEN[i] != b.DEN[i]) {
			cout << "Different denominator contribution at prime " << i << ": "
					<< DEN[i] << "!=" << b.DEN[i] << endl;
			return false;
		}
	cout << "Well, they're simply equal!" << endl;
	return true;
}
;

char*
TFracNum::HeaderString() {
	char* hstr = new char[30];
	if (sign_prefac == 0) {
		sprintf(hstr, "{0,1}");
		return hstr;
	}
	if (sign_prefac == -7777) {
		sprintf(hstr, "{1,0}");
		return hstr;
	}
	if (sign_prefac == -6666) {
		sprintf(hstr, "{0,0}");
		return hstr;
	}
	if (sign_prefac == 1)
		sprintf(hstr, "{%ld,%ld}", NOM_INT, DEN_INT);
	else
		sprintf(hstr, "{%ld,%ld}", -NOM_INT, DEN_INT);
	return hstr;
}
;

bool TFracNum::operator>(const TFracNum &b) const {
	if (dvalue > b.dvalue)
		return true;
	return false;
}
;

TFracNum TFracNum::operator+(const TFracNum &b) const {
	long den_cdiv = DenomCommonDivisor(b);
	long bdc = b.DEN_INT / den_cdiv;
	long adc = DEN_INT / den_cdiv;
	return TFracNum(
			sign_prefac * NOM_INT * bdc + b.sign_prefac * b.NOM_INT * adc,
			DEN_INT * bdc);
}

TFracNum TFracNum::operator*(const TFracNum &b) const {
	// if one of the two numbers is undetermined,
	// the product is also undetermined
	if (sign_prefac == -6666 || b.sign_prefac == -6666)
		return TFracNum(0, 0, 0, 0, -6666);

	// if one of the two numbers contains division by zero,
	// and the other nominator is zero, the product is undetermined
	if ((sign_prefac == -7777 && b.sign_prefac == 0)
			|| (sign_prefac == 0 && b.sign_prefac == -7777))
		return TFracNum(0, 0, 0, 0, -6666);

	// other cases with division by zero; product is also infinity
	if ((sign_prefac == -7777 || b.sign_prefac == -7777))
		return TFracNum(0, 0, 0, 0, -7777);

	if (sign_prefac * b.sign_prefac == 0)
		return TFracNum(0, 0, 0, 0, 0);

	long maxPrimNom_ = maxPrimNom;
	if (b.maxPrimNom > maxPrimNom_)
		maxPrimNom_ = b.maxPrimNom;
	long maxPrimDen_ = maxPrimDen;
	if (b.maxPrimDen > maxPrimDen_)
		maxPrimDen_ = b.maxPrimDen;
	long maxPrim = maxPrimNom_;
	if (maxPrimDen_ > maxPrim)
		maxPrim = maxPrimDen_;

	long prim_vec_nom[maxPrim];
	long prim_vec_den[maxPrim];

	for (long ip = 0; ip < maxPrim; ip++) {
		prim_vec_nom[ip] = 0;
		prim_vec_den[ip] = 0;
		if (maxPrimNom > ip)
			prim_vec_nom[ip] += NOM[ip];
		if (b.maxPrimNom > ip)
			prim_vec_nom[ip] += b.NOM[ip];
		if (maxPrimDen > ip)
			prim_vec_den[ip] += DEN[ip];
		if (b.maxPrimDen > ip)
			prim_vec_den[ip] += b.DEN[ip];
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
		if (prim_vec_nom[ip] != 0)
			maxPrimNom_ = ip + 1;
		if (prim_vec_den[ip] != 0)
			maxPrimDen_ = ip + 1;
	}

	long* NOM_ = new long[maxPrimNom_];
	for (long jp = 0; jp < maxPrimNom_; jp++)
		NOM_[jp] = prim_vec_nom[jp];

	long* DEN_ = new long[maxPrimDen_];
	for (long jp = 0; jp < maxPrimDen_; jp++)
		DEN_[jp] = prim_vec_den[jp];

	if (debugFracNum) {
		cout << " Initial with maxN=" << maxPrimNom_ << ", maxD=" << maxPrimDen_
				<< ", " << NOM_ << ", " << DEN_ << ", "
				<< sign_prefac * b.sign_prefac << endl;
		cout << "NOM:";
		for (long i = 0; i < maxPrimNom_; i++)
			cout << NOM_[i] << ",";
		cout << endl;
		cout << "DEN:";
		for (long i = 0; i < maxPrimDen_; i++)
			cout << DEN_[i] << ",";
		cout << endl;
	}
	return TFracNum(maxPrimNom_, maxPrimDen_, NOM_, DEN_,
			sign_prefac * b.sign_prefac);
}

double TFracNum::Print() const {
	if (debugFracNum) {
		cout << "nom prime list: " << maxPrimNom << ",pointer " << NOM << endl;
		cout << "den prime list: " << maxPrimDen << ",pointer " << DEN << endl;
		cout << "NOM:";
		for (long i = 0; i < maxPrimNom; i++)
			cout << NOM[i] << ",";
		cout << endl;
		cout << "DEN:";
		for (long i = 0; i < maxPrimDen; i++)
			cout << DEN[i] << ",";
		cout << endl;
	}
	cout << "sign_prefac=" << sign_prefac << endl;
	if (maxPrimNom < 0) {
		if (sign_prefac < 0)
			cout << "-";
		cout << -maxPrimNom << "/" << -maxPrimDen << endl;
		return double(sign_prefac) * double(-maxPrimNom) / double(-maxPrimDen);
	}
	long integrity = 1;
	for (long i = 0; i < maxPrimNom; i++)
		if (NOM[i] < 0 || NOM[i] > 1000)
			integrity = 0;
	for (long i = 0; i < maxPrimDen; i++)
		if (DEN[i] < 0 || DEN[i] > 1000)
			integrity = 0;
	if (integrity == 0)
		return -1;

	long nom = 1;
	long ipn = 0;
	if (sign_prefac < 0)
		cout << "-NOM = ";
	else
		cout << " NOM = ";
	long FirstTerm = 1;
	while (ipn < maxPrimNom) {
		if (NOM[ipn] != 0) {
			cout << PRIMES[ipn] << "^" << NOM[ipn];
			FirstTerm = 0;
		}
		for (long jj = 0; jj < NOM[ipn]; jj++)
			nom *= PRIMES[ipn];
		ipn++;
		if (!FirstTerm && ipn < maxPrimNom && NOM[ipn] != 0)
			cout << " * ";
	}
	cout << " = " << nom << endl;

	long den = 1;
	long ipd = 0;
	cout << " DEN = ";
	FirstTerm = 1;
	while (ipd < maxPrimDen) {
		if (DEN[ipd] != 0) {
			cout << PRIMES[ipd] << "^" << DEN[ipd];
			FirstTerm = 0;
		}
		for (long jj = 0; jj < DEN[ipd]; jj++)
			den *= PRIMES[ipd];
		ipd++;
		if (!FirstTerm && ipd < maxPrimDen && DEN[ipd] != 0)
			cout << " * ";
	}
	cout << " = " << den << endl;
	cout << "NOM_INT=" << NOM_INT << endl;
	cout << "DEN_INT=" << DEN_INT << endl;
	cout << "dvalue=" << dvalue << endl;
	return double(sign_prefac) * double(nom) / double(den);
}

double TFracNum::PrintToErr() const {
	cerr << "nom prime list: " << maxPrimNom << ",pointer " << NOM << endl;
	cerr << "den prime list: " << maxPrimDen << ",pointer " << DEN << endl;
	cerr << "NOM:";
	for (long i = 0; i < maxPrimNom; i++)
		cerr << NOM[i] << ",";
	cerr << endl;
	cerr << "DEN:";
	for (long i = 0; i < maxPrimDen; i++)
		cerr << DEN[i] << ",";
	cerr << endl;
	cerr << "sign_prefac=" << sign_prefac << endl;
	if (maxPrimNom < 0) {
		if (sign_prefac < 0)
			cerr << "-";
		cerr << -maxPrimNom << "/" << -maxPrimDen << endl;
		return double(sign_prefac) * double(-maxPrimNom) / double(-maxPrimDen);
	}
	long integrity = 1;
	for (long i = 0; i < maxPrimNom; i++)
		if (NOM[i] < 0 || NOM[i] > 1000)
			integrity = 0;
	for (long i = 0; i < maxPrimDen; i++)
		if (DEN[i] < 0 || DEN[i] > 1000)
			integrity = 0;
	if (integrity == 0)
		return -1;

	long nom = 1;
	long ipn = 0;
	if (sign_prefac < 0)
		cerr << "-NOM = ";
	else
		cerr << " NOM = ";
	long FirstTerm = 1;
	while (ipn < maxPrimNom) {
		if (NOM[ipn] != 0) {
			cerr << PRIMES[ipn] << "^" << NOM[ipn];
			FirstTerm = 0;
		}
		for (long jj = 0; jj < NOM[ipn]; jj++)
			nom *= PRIMES[ipn];
		ipn++;
		if (!FirstTerm && ipn < maxPrimNom && NOM[ipn] != 0)
			cerr << " * ";
	}
	cerr << " = " << nom << endl;

	long den = 1;
	long ipd = 0;
	cerr << " DEN = ";
	FirstTerm = 1;
	while (ipd < maxPrimDen) {
		if (DEN[ipd] != 0) {
			cerr << PRIMES[ipd] << "^" << DEN[ipd];
			FirstTerm = 0;
		}
		for (long jj = 0; jj < DEN[ipd]; jj++)
			den *= PRIMES[ipd];
		ipd++;
		if (!FirstTerm && ipd < maxPrimDen && DEN[ipd] != 0)
			cerr << " * ";
	}
	cerr << " = " << den << endl;
	cerr << "NOM_INT=" << NOM_INT << endl;
	cerr << "DEN_INT=" << DEN_INT << endl;
	cerr << "dvalue=" << dvalue << endl;
	return double(sign_prefac) * double(nom) / double(den);
}

const char*
TFracNum::FracString() {
	char *formstr = new char[50];
	char *fstr = new char[100];
	if (NOM_INT == 0) {
		sprintf(fstr, "0");
	} else if (DEN_INT == 1) {
		sprintf(formstr, "%%c%s", IOUTSTRING);
		sprintf(fstr, formstr, sign_prefac < 0 ? '-' : '+', NOM_INT);
	} else {
		sprintf(formstr, "%%c%s/%s", IOUTSTRING, IOUTSTRING);
		sprintf(fstr, formstr, sign_prefac < 0 ? '-' : '+', NOM_INT, DEN_INT);
	}
	return fstr;
}
;

const char NULLSTRING[1] = ""; // workaround CINT warning when sprintf(s,"");

const char*
TFracNum::FracStringSqrt() const {
	char *formstr = new char[50];
	char *fstr = new char[200];
	if (NOM_INT == 0) {
		sprintf(fstr, "0");
		return fstr;
	}
	long ipn = 0;
	long SQRT_NOM_INT = 1;
	long NOM_INT_REST = 1;
	while (ipn < maxPrimNom) {
		for (long jj = 0; jj < NOM[ipn] / 2; jj++)
			SQRT_NOM_INT *= PRIMES[ipn];
		if (NOM[ipn] % 2)
			NOM_INT_REST *= PRIMES[ipn];
		ipn++;
	}
	long ipd = 0;
	long SQRT_DEN_INT = 1;
	long DEN_INT_REST = 1;
	while (ipd < maxPrimDen) {
		for (long jj = 0; jj < DEN[ipd] / 2; jj++)
			SQRT_DEN_INT *= PRIMES[ipd];
		if (DEN[ipd] % 2)
			DEN_INT_REST *= PRIMES[ipd];
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

	if (one1 && one2)
		sprintf(sqrtstr, "1");
	sprintf(fstr, "%c%s%s", sign_prefac < 0 ? '-' : '+', sqrtstr, reststr);
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
	if (ell == 0)
		return TFracNum(0, 0, 0, 0, 1);
	long nom_ptr[1] = { ell };
	TFracNum two_to_ell(1, 0, nom_ptr, 0, 1);
	TFracNum fac1(ell, 1, "factorial");
	TFracNum fac2(ell, 2 * ell, "factorial");
	return two_to_ell * fac1 * fac2;
}

TFracNum cm0_sub_ell(long ell, long m0) {
	if (ell == 0)
		return TFracNum(0, 0, 0, 0, 1);
	long nom_ptr[1] = { (ell + m0) / 2 };
	TFracNum two_to_ell(1, 0, nom_ptr, 0, 1);
	TFracNum fac1(ell, 1, "factorial");
	TFracNum fac2(ell, 2 * ell, "factorial");
	return two_to_ell * fac1 * fac2;
}

TFracNum cm0_sub_ell_2(long ell, long m0) {
	//return  am0_to_J(ell, 0, m0);
	if (ell == 0)
		return TFracNum(0, 0, 0, 0, 1);
	long nom_ptr[1] = { (ell + m0) };
	TFracNum two_to_ell(1, 0, nom_ptr, 0, 1);
	TFracNum fac1a(ell, 1, "factorial");
	TFracNum fac1b(ell, 1, "factorial");
	TFracNum fac2a(ell, 2 * ell, "factorial");
	TFracNum fac2b(ell, 2 * ell, "factorial");
	return two_to_ell * fac1a * fac1b * fac2a * fac2b;
}

