#include "TJwfTensor.h"

#include <iostream>
#include <string>

#include <reportingUtils.hpp>

using namespace std;

TTensorTerm::TTensorTerm(char name,
                         long RJ,
                         long* pzm_field,
                         const TFracNum& prefac)
	: _Rome(0),
	  _ome_pzm(0),
	  _Reps(0),
	  _eps_pzm(0),
	  _Rchi(0),
	  _chi_pzm(0),
	  _Rphi(0),
	  _phi_pzm(0),
	  _gam_s_pot(0),
	  _gam_sig_pot(0),
	  _prefac(prefac)
{
	switch(name) {
		case 'o':
			_Rome = RJ;
			_ome_pzm = pzm_field;
			break;
		case 'e':
			_Reps = RJ;
			_eps_pzm = pzm_field;
			break;
		case 'c':
			_Rchi = RJ;
			_chi_pzm = pzm_field;
			break;
		case 'p':
			_Rphi = RJ;
			_phi_pzm = pzm_field;
			break;
	}
}

TTensorTerm::TTensorTerm(const TTensorTerm& S,
                         const TTensorTerm& L,
                         long contractions,
                         long o_share,
                         long e_share,
                         char con_type)
	: _Rome(S._Rome),
	  _Reps(S._Reps),
	  _Rchi(S._Rchi + L._Rchi),
	  _Rphi(S._Rphi + L._Rphi),
	  _gam_s_pot(S._gam_s_pot + L._gam_s_pot),
	  _gam_sig_pot(S._gam_sig_pot + L._gam_sig_pot),
	  _prefac(S._prefac * L._prefac)
{

	if (L._Rome or L._Reps) {
		return;
	}

	long ocount = o_share;
	long ecount = e_share;
	for (long con = 0; con < contractions; con++) {
		long cp = 0;
		if (con_type == 'c') {
			cp = L._chi_pzm[_Rchi - 1];
			_Rchi--;
		} else {
			cp = L._phi_pzm[_Rphi - 1];
			_Rphi--;
		}
		long oe = 0;
		if (ocount) {              // o is to be contracted
			oe = S._ome_pzm[_Rome - 1];
			_Rome--;
		} else {                     // e is to be contracted
			oe = S._eps_pzm[_Reps - 1];
			ecount--;
			_Reps--;
		}
		if (con_type == 'c' and ((oe == 1 and cp == -1) or (oe == -1 and cp == 1)))
		{
			_prefac.FlipSign();
		} else if (con_type == 'p' and ((oe == 1 and cp == 1) or (oe == -1 and cp == -1))) {

		} else if (oe == 0 and cp == 0) {
			if (ocount) {
				_gam_s_pot++;
			} else {
				_gam_sig_pot++;
			}
		} else {
			_prefac = TFracNum::Zero;
			return;
		}
		if (ocount) {
			ocount--;
		}
	}

	_ome_pzm = new long[_Rome];
	for (long i = 0; i < _Rome; i++) {
		_ome_pzm[i] = S._ome_pzm[i];
	}
	_eps_pzm = new long[_Reps];
	for (long i = 0; i < _Reps; i++) {
		_eps_pzm[i] = S._eps_pzm[i];
	}
	_chi_pzm = new long[_Rchi];
	for (long i = 0; i < _Rchi; i++) {
		if (i < L._Rchi) {
			_chi_pzm[i] = L._chi_pzm[i];
		} else {
			_chi_pzm[i] = S._chi_pzm[i - (L._Rchi)];
		}
	}
	_phi_pzm = new long[_Rphi];
	for (long i = 0; i < _Rphi; i++) {
		if (i < L. _Rphi) {
			_phi_pzm[i] = L._phi_pzm[i];
		} else {
			_phi_pzm[i] = S._phi_pzm[i - (L._Rphi)];
		}
	}
}

long TTensorTerm::LJContraction(long ncon, long even) {

	for (long con = 0; con < ncon; con++) {
		long c = _chi_pzm[_Rchi - 1];
		long p = _phi_pzm[_Rphi - 1];
		if (c != p) {
			_prefac = TFracNum::Zero;
			return 0;
		}
		_Rchi--;
		_Rphi--;
	}

	bool error = false;

	if (even == 0) {
		if (_Rome + _Reps + _Rchi + _Rphi != 3) {
			printErr << "TTensorTerm::LJContraction:"
			         << " Contraction ended with wrong number of indices!!"
			         << endl;
			throw;
		} else { // eq. (5.21) - (5.23)
			long os = -100;
			if (_Rome == 1) {
				os = _ome_pzm[0];
			}
			long es = -100;
			if (_Reps == 1) {
				es = _eps_pzm[0];
			}
			long cs = -100;
			if (_Rchi == 1) {
				cs = _chi_pzm[0];
			}
			long ps = -100;
			if (_Rphi == 1) {
				ps = _phi_pzm[0];
			}
			if (_Rome == 0) {
				_Reps--;
				_Rchi--;
				_Rphi--;
				if (es == 1 && cs == -1 && ps == 0) {

				} else if (es == -1 && cs == 1 && ps == 0) {
					_prefac.FlipSign();
				} else if (es == 1 && cs == 0 && ps == 1) {

				} else if (es == -1 && cs == 0 && ps == -1) {
					_prefac.FlipSign();
				} else if (es == 0 && cs == 1 && ps == 1) {
					_prefac.FlipSign();
					_gam_sig_pot++;
				} else if (es == 0 && cs == -1 && ps == -1) {
					_gam_sig_pot++;
				} else {
					_prefac = TFracNum::Zero;
					return 1;
				}
			} else if (_Reps == 0) {
				_Rome--;
				_Rchi--;
				_Rphi--;
				if (os == 1 && cs == -1 && ps == 0) {

				} else if (os == -1 && cs == 1 && ps == 0) {
					_prefac.FlipSign();
				} else if (os == 1 && cs == 0 && ps == 1) {

				} else if (os == -1 && cs == 0 && ps == -1) {
					_prefac.FlipSign();
				} else if (os == 0 && cs == 1 && ps == 1) {
					_prefac.FlipSign();
					_gam_s_pot++;
				} else if (os == 0 && cs == -1 && ps == -1) {
					_gam_s_pot++;
				} else {
					_prefac = TFracNum::Zero;
					return 1;
				}
			} else if (_Rchi == 0) {
				_Rome--;
				_Reps--;
				_Rphi--;
				if (os == 1 && es == -1 && ps == 0) {

				} else if (os == -1 && es == 1 && ps == 0) {
					_prefac.FlipSign();
				} else if (os == 1 && es == 0 && ps == 1) {
					_gam_sig_pot++;
				} else if (os == -1 && es == 0 && ps == -1) {
					_prefac.FlipSign();
					_gam_sig_pot++;
				} else if (os == 0 && es == 1 && ps == 1) {
					_prefac.FlipSign();
					_gam_s_pot++;
				} else if (os == 0 && es == -1 && ps == -1) {
					_gam_s_pot++;
				} else {
					_prefac = TFracNum::Zero;
					return 1;
				}
			} else if (_Rphi == 0) {
				_Rome--;
				_Reps--;
				_Rchi--;
				if (os == 1 && es == -1 && cs == 0) {

				} else if (os == -1 && es == 1 && cs == 0) {
					_prefac.FlipSign();
				} else if (os == 1 && es == 0 && cs == -1) {
					_prefac.FlipSign();
					_gam_sig_pot++;
				} else if (os == -1 && es == 0 && cs == 1) {
					_gam_sig_pot++;
				} else if (os == 0 && es == 1 && cs == -1) {
					_gam_s_pot++;
				} else if (os == 0 && es == -1 && cs == 1) {
					_prefac.FlipSign();
					_gam_s_pot++;
				} else {
					_prefac = TFracNum::Zero;
					return 1;
				}
			} else {
				printWarn << "troule == 5, whatever that means..." << endl;
				error = true;
			}
		}
	}

	else {
		if (_Rome != 0 and _Reps != 0 and _Rchi != 0 and _Rphi != 0) {
			printWarn << "troule == 1, whatever that means..." << endl;
			error = true;
		}
	}
	if (error) {
		printErr << "TTensorTerm::LJContraction: "
		         << "Invalid espilon-contraction occurred " << endl;
		throw;
	}
	return 1;
}

long TTensorTerm::Multiply(char name,
                           long RJ,
                           long* pzm_field,
                           const TFracNum& prefac)
{

	_prefac = _prefac * prefac;

	long Merr = 0;
	if (name == 'o') {
		if (_Rome) {
			Merr = 1;
		} else {
			_Rome = RJ;
			_ome_pzm = pzm_field;
		}
	}
	if (name == 'e') {
		if (_Reps) {
			Merr = 1;
		} else {
			_Reps = RJ;
			_eps_pzm = pzm_field;
		}
	}
	if (name == 'c') {
		if (_Rchi) {
			Merr = 1;
		} else {
			_Rchi = RJ;
			_chi_pzm = pzm_field;
		}
	}
	if (name == 'p') {
		if (_Rphi) {
			Merr = 1;
		} else {
			_Rphi = RJ;
			_phi_pzm = pzm_field;
		}
	}
	if (Merr) {
		printErr << "TTensorTerm::Multiply: Each type can be multiplied only once!"
		         << endl;
		return 0;
	}
	return 1;
}

long TTensorTerm::SpinInnerContraction(long cPsiInt) {
	long res = 0;
	for (long ic = 0; ic < cPsiInt; ic++) {
		res = 0;
		long o = _ome_pzm[_Rome - 1];
		long e = _eps_pzm[_Reps - 1];
		if ((o == 1 and e == -1) or (o == -1 and e == 1)) {
			_prefac.FlipSign();
			res = 1;
		}
		if (o == 0 and e == 0) {
			_gam_s_pot += 1;
			_gam_sig_pot += 1;
			res = 1;
		}
		_Rome--;
		_Reps--;
	}
	return res;
}

bool TTensorTerm::SameStructure(const TTensorTerm& other) const
{
	if (_Rome == other._Rome           and
	    _Reps == other._Reps           and
	    _Rchi == other._Rchi           and
	    _Rphi == other._Rphi           and
	    _gam_s_pot == other._gam_s_pot and
	    _gam_sig_pot == other._gam_sig_pot)
	{
		return true;
	}
	return false;
}

bool TTensorTerm::AddTwoTerms(const TTensorTerm& other) {
	if (!SameStructure(other)) {
		printErr << "NO NO NO these terms cannot be added!" << endl;
		return false;
	} else {
		const TFracNum* sum = _prefac.SumSignedRoots(other._prefac);
		if (sum) {
			_prefac = *sum;
			return true;
		}
	}
	return false;
}

long TTensorTerm::Print(char flag) const
{
	if (flag == 's') {
		cout << _prefac.FracStringSqrt() << " ";
	} else {
		cout << _prefac.FracString() << " ";
	}
	if (_Rome) {
		cout << "o(";
		for (long i = 0; i < _Rome; i++) {
			if (_ome_pzm[i] == 1) {
				cout << "+";
			}
			if (_ome_pzm[i] == 0) {
				cout << "0";
			}
			if (_ome_pzm[i] == -1) {
				cout << "-";
			}
		}
		cout << ")";
	}
	if (_Reps) {
		cout << "e(";
		for (long i = 0; i < _Reps; i++) {
			if (_eps_pzm[i] == 1) {
				cout << "+";
			}
			if (_eps_pzm[i] == 0) {
				cout << "0";
			}
			if (_eps_pzm[i] == -1) {
				cout << "-";
			}
		}
		cout << ")";
	}
	if (_Rchi) {
		cout << "c(";
		for (long i = 0; i < _Rchi; i++) {
			if (_chi_pzm[i] == 1) {
				cout << "+";
			}
			if (_chi_pzm[i] == 0) {
				cout << "0";
			}
			if (_chi_pzm[i] == -1) {
				cout << "-";
			}
		}
		cout << ")";
	}
	if (_Rphi) {
		cout << "p(";
		for (long i = 0; i < _Rphi; i++) {
			if (_phi_pzm[i] == 1) {
				cout << "+";
			}
			if (_phi_pzm[i] == 0) {
				cout << "0";
			}
			if (_phi_pzm[i] == -1) {
				cout << "-";
			}
		}
		cout << ")";
	}
	if (_gam_s_pot) {
		if (_gam_s_pot == 1) {
			cout << " gs";
		} else {
			cout << " gs^" << _gam_s_pot;
		}
	}
	if (_gam_sig_pot) {
		if (_gam_sig_pot == 1) {
			cout << " gsig";
		} else {
			cout << " gsig^" << _gam_sig_pot;
		}
	}
	return 0;
}

//long
//TTensorSum::Print(char flag='n'){ // CINT limitation for overloading
long TTensorSum::Print(char flag) {
	for (long i = 0; i < Nterms; i++) {
		terms[i].Print(flag);
		if (i < Nterms - 1) {
			cout << " ";
		}
	}
	cout << endl;
	return 0;
}
;

long TTensorSum::AddTerm(TTensorTerm* addt) {
	TTensorTerm* nt = new TTensorTerm[Nterms + 1];
	for (long i = 0; i < Nterms; i++) {
		nt[i] = terms[i];
	}
	nt[Nterms] = *addt;
	delete[] terms;
	terms = nt;
	Nterms++;
	return 0;
}

long TTensorSum::SpinInnerContraction(long cPsiInt) {
	long non_zero_terms = 0;
	long val[Nterms];
	for (long i = 0; i < Nterms; i++) {
		val[i] = terms[i].SpinInnerContraction(cPsiInt);
		if (val[i] != 0) {
			non_zero_terms++;
		}
	}
	if (non_zero_terms < Nterms) {
		TTensorTerm* nt = new TTensorTerm[non_zero_terms];
		long j = 0;
		for (long i = 0; i < Nterms; i++) {
			if (val[i]) {
				nt[j] = terms[i];
				j++;
			}
		}
		delete[] terms;
		terms = nt;
		Nterms = non_zero_terms;
	}
	return Nterms;
}

TTensorSum*
TTensorSum::LSContraction(TTensorSum *L, long contr,
long co, long ce, char con_type) {

	TTensorSum *tls = new TTensorSum();

	for (long i = 0; i < Nterms; i++) {
		for (long j = 0; j < L->Nterms; j++) {
			TTensorTerm *nt = new TTensorTerm(terms[i], L->terms[j], contr, co, ce, con_type);
			if (nt->IsNonZero()) {
				tls->AddTerm(nt);
			}
		}
	}

	return tls;
}

TTensorSum*
TTensorSum::LJContraction(long cChiPhi, long even) {

	for (long i = 0; i < Nterms; i++) {
		terms[i].LJContraction(cChiPhi, even);
	}

	TTensorSum *tls = new TTensorSum();

	for (long i = 0; i < Nterms; i++) {
		long found_same = 0;
		for (long j = 0; j < tls->Nterms; j++) {
			if (terms[i].SameStructure(tls->terms[j])) {
				if ((tls->terms[j]).AddTwoTerms(terms[i])) {
					found_same = 1;
					break;
				}
			}
		}
		if (!found_same) {
			TTensorTerm *nt = new TTensorTerm(terms[i]);
			tls->AddTerm(nt);
		}
	}

	TTensorSum* tls_nonzero = new TTensorSum();

	for (long i = 0; i < tls->Nterms; i++) {
		if (tls->terms[i].IsNonZero()) {
			TTensorTerm *nt = new TTensorTerm(tls->terms[i]);
			tls_nonzero->AddTerm(nt);
		}
	}

	return tls_nonzero;
}
