#include "TJwfTensor.h"

#include <string>

#include "ClebschGordanBox.h"

#include <reportingUtils.hpp>

using namespace std;
using rpwa::operator<<;


TTensorTerm::TTensorTerm(const char& name,
                         const vector<long>& pzm_field,
                         const TFracNum& prefac)
	: _ome_pzm(),
	  _eps_pzm(),
	  _chi_pzm(),
	  _phi_pzm(),
	  _gam_s_pot(0),
	  _gam_sig_pot(0),
	  _prefac(prefac)
{
	switch(name) {
		case 'o':
			_ome_pzm = pzm_field;
			break;
		case 'e':
			_eps_pzm = pzm_field;
			break;
		case 'c':
			_chi_pzm = pzm_field;
			break;
		case 'p':
			_phi_pzm = pzm_field;
			break;
	}
}

TTensorTerm::TTensorTerm(const TTensorTerm& S,
                         const TTensorTerm& L,
                         const long& contractions,
                         long o_share,
                         long e_share,
                         const char& conType)
	: _ome_pzm(),
	  _eps_pzm(),
	  _chi_pzm(),
	  _phi_pzm(),
	  _gam_s_pot(S._gam_s_pot + L._gam_s_pot),
	  _gam_sig_pot(S._gam_sig_pot + L._gam_sig_pot),
	  _prefac(S._prefac * L._prefac)
{

	if (not L._ome_pzm.empty() or not L._eps_pzm.empty()) {
		return;
	}

	size_t rOme = S._ome_pzm.size();
	size_t rEps = S._eps_pzm.size();
	size_t rChi = L._chi_pzm.size() + S._chi_pzm.size();
	size_t rPhi = L._phi_pzm.size() + S._phi_pzm.size();

	for (long con = 0; con < contractions; con++) {
		long cp = 0;
		if (conType == 'c') {
			cp = L._chi_pzm[rChi - 1];
			rChi--;
		} else {
			cp = L._phi_pzm[rPhi - 1];
			rPhi--;
		}
		long oe = 0;
		if (o_share) {              // o is to be contracted
			oe = S._ome_pzm[rOme - 1];
			rOme--;
		} else {                     // e is to be contracted
			oe = S._eps_pzm[rEps - 1];
			e_share--;
			rEps--;
		}
		if (conType == 'c' and ((oe == 1 and cp == -1) or (oe == -1 and cp == 1)))
		{
			_prefac.FlipSign();
		} else if (conType == 'p' and ((oe == 1 and cp == 1) or (oe == -1 and cp == -1))) {

		} else if (oe == 0 and cp == 0) {
			if (o_share) {
				_gam_s_pot++;
			} else {
				_gam_sig_pot++;
			}
		} else {
			_prefac = TFracNum::Zero;
			return;
		}
		if (o_share) {
			o_share--;
		}
	}

	_ome_pzm.resize(rOme);
	for (size_t i = 0; i < _ome_pzm.size(); i++) {
		_ome_pzm[i] = S._ome_pzm[i];
	}
	_eps_pzm.resize(rEps);
	for (size_t i = 0; i < _eps_pzm.size(); i++) {
		_eps_pzm[i] = S._eps_pzm[i];
	}
	_chi_pzm.resize(rChi);
	for (size_t i = 0; i < _chi_pzm.size(); i++) {
		if (i < L._chi_pzm.size()) {
			_chi_pzm[i] = L._chi_pzm[i];
		} else {
			_chi_pzm[i] = S._chi_pzm[i - L._chi_pzm.size()];
		}
	}
	_phi_pzm.resize(rPhi);
	for (size_t i = 0; i < _phi_pzm.size(); i++) {
		if (i < L._phi_pzm.size()) {
			_phi_pzm[i] = L._phi_pzm[i];
		} else {
			_phi_pzm[i] = S._phi_pzm[i - L._phi_pzm.size()];
		}
	}

}

long TTensorTerm::LJContraction(const long& nCon, const bool& even) {

	size_t rOme = _ome_pzm.size();
	size_t rEps = _eps_pzm.size();
	size_t rChi = _chi_pzm.size();
	size_t rPhi = _phi_pzm.size();

	for (long con = 0; con < nCon; con++) {
		const long& c = _chi_pzm[rChi - 1];
		const long& p = _phi_pzm[rPhi - 1];
		if (c != p) {
			_prefac = TFracNum::Zero;
			shrinkVectors(rOme, rEps, rChi, rPhi);
			return 0;
		}
		rChi--;
		rPhi--;
	}

	bool error = false;

	if (not even) {
		if (rOme + rEps + rChi + rPhi != 3) {
			printErr << "TTensorTerm::LJContraction:"
			         << " Contraction ended with wrong number of indices!!"
			         << endl;
			throw;
		} else { // eq. (5.21) - (5.23)
			long os = -100;
			if (rOme == 1) {
				os = _ome_pzm[0];
			}
			long es = -100;
			if (rEps == 1) {
				es = _eps_pzm[0];
			}
			long cs = -100;
			if (rChi == 1) {
				cs = _chi_pzm[0];
			}
			long ps = -100;
			if (rPhi == 1) {
				ps = _phi_pzm[0];
			}
			if (rOme == 0) {
				rEps--;
				rChi--;
				rPhi--;
				if (es == 1 and cs == -1 and ps == 0) {

				} else if (es == -1 and cs == 1 and ps == 0) {
					_prefac.FlipSign();
				} else if (es == 1 and cs == 0 and ps == 1) {

				} else if (es == -1 and cs == 0 and ps == -1) {
					_prefac.FlipSign();
				} else if (es == 0 and cs == 1 and ps == 1) {
					_prefac.FlipSign();
					_gam_sig_pot++;
				} else if (es == 0 and cs == -1 and ps == -1) {
					_gam_sig_pot++;
				} else {
					_prefac = TFracNum::Zero;
					shrinkVectors(rOme, rEps, rChi, rPhi);
					return 1;
				}
			} else if (rEps == 0) {
				rOme--;
				rChi--;
				rPhi--;
				if (os == 1 and cs == -1 and ps == 0) {

				} else if (os == -1 and cs == 1 and ps == 0) {
					_prefac.FlipSign();
				} else if (os == 1 and cs == 0 and ps == 1) {

				} else if (os == -1 and cs == 0 and ps == -1) {
					_prefac.FlipSign();
				} else if (os == 0 and cs == 1 and ps == 1) {
					_prefac.FlipSign();
					_gam_s_pot++;
				} else if (os == 0 and cs == -1 and ps == -1) {
					_gam_s_pot++;
				} else {
					_prefac = TFracNum::Zero;
					shrinkVectors(rOme, rEps, rChi, rPhi);
					return 1;
				}
			} else if (rChi == 0) {
				rOme--;
				rEps--;
				rPhi--;
				if (os == 1 and es == -1 and ps == 0) {

				} else if (os == -1 and es == 1 and ps == 0) {
					_prefac.FlipSign();
				} else if (os == 1 and es == 0 and ps == 1) {
					_gam_sig_pot++;
				} else if (os == -1 and es == 0 and ps == -1) {
					_prefac.FlipSign();
					_gam_sig_pot++;
				} else if (os == 0 and es == 1 and ps == 1) {
					_prefac.FlipSign();
					_gam_s_pot++;
				} else if (os == 0 and es == -1 and ps == -1) {
					_gam_s_pot++;
				} else {
					_prefac = TFracNum::Zero;
					shrinkVectors(rOme, rEps, rChi, rPhi);
					return 1;
				}
			} else if (rPhi == 0) {
				rOme--;
				rEps--;
				rChi--;
				if (os == 1 and es == -1 and cs == 0) {

				} else if (os == -1 and es == 1 and cs == 0) {
					_prefac.FlipSign();
				} else if (os == 1 and es == 0 and cs == -1) {
					_prefac.FlipSign();
					_gam_sig_pot++;
				} else if (os == -1 and es == 0 and cs == 1) {
					_gam_sig_pot++;
				} else if (os == 0 and es == 1 and cs == -1) {
					_gam_s_pot++;
				} else if (os == 0 and es == -1 and cs == 1) {
					_prefac.FlipSign();
					_gam_s_pot++;
				} else {
					_prefac = TFracNum::Zero;
					shrinkVectors(rOme, rEps, rChi, rPhi);
					return 1;
				}
			} else {
				printWarn << "troule == 5, whatever that means..." << endl;
				error = true;
			}
		}
	} else {
		if (rOme != 0 and rEps != 0 and rChi != 0 and rPhi != 0) {
			printWarn << "troule == 1, whatever that means..." << endl;
			error = true;
		}
	}
	if (error) {
		printErr << "TTensorTerm::LJContraction: "
		         << "Invalid espilon-contraction occurred " << endl;
		throw;
	}
	shrinkVectors(rOme, rEps, rChi, rPhi);
	return 1;
}

void TTensorTerm::Multiply(const char& name,
                           const vector<long>& pzm_field,
                           const TFracNum& prefac)
{

	_prefac = _prefac * prefac;

	bool error = false;
	if (name == 'o') {
		if (not _ome_pzm.empty()) {
			error = true;
		} else {
			_ome_pzm = pzm_field;
		}
	}
	if (name == 'e') {
		if (not _eps_pzm.empty()) {
			error = true;
		} else {
			_eps_pzm = pzm_field;
		}
	}
	if (name == 'c') {
		if (not _chi_pzm.empty()) {
			error = true;
		} else {
			_chi_pzm = pzm_field;
		}
	}
	if (name == 'p') {
		if (not _phi_pzm.empty()) {
			error = true;
		} else {
			_phi_pzm = pzm_field;
		}
	}
	if (error) {
		printErr << "TTensorTerm::Multiply: Each type can be multiplied only once!" << endl;
		throw;
	}
}

long TTensorTerm::SpinInnerContraction(const long& cPsiInt) {
	long res = 0;
	size_t rOme = _ome_pzm.size();
	size_t rEps = _eps_pzm.size();
	for (long ic = 0; ic < cPsiInt; ic++) {
		res = 0;
		long o = _ome_pzm[rOme - 1];
		long e = _eps_pzm[rEps - 1];
		if ((o == 1 and e == -1) or (o == -1 and e == 1)) {
			_prefac.FlipSign();
			res = 1;
		}
		if (o == 0 and e == 0) {
			_gam_s_pot += 1;
			_gam_sig_pot += 1;
			res = 1;
		}
		rOme--;
		rEps--;
	}
	_ome_pzm.resize(rOme);
	_eps_pzm.resize(rEps);
	return res;
}

bool TTensorTerm::SameStructure(const TTensorTerm& other) const
{
	if (_ome_pzm.size() == other._ome_pzm.size() and
	    _eps_pzm.size() == other._eps_pzm.size() and
	    _chi_pzm.size() == other._chi_pzm.size() and
	    _phi_pzm.size() == other._phi_pzm.size() and
	    _gam_s_pot == other._gam_s_pot             and
	    _gam_sig_pot == other._gam_sig_pot)
	{
		return true;
	}
	return false;
}

void TTensorTerm::AddTwoTerms(const TTensorTerm& other) {
	if (not SameStructure(other)) {
		printErr << "NO NO NO these terms cannot be added!" << endl;
		throw;
	}
	_prefac = _prefac.SumSignedRoots(other._prefac);
}

std::ostream& TTensorTerm::Print(const char& flag, std::ostream& out) const
{
	if (flag == 's') {
		out << _prefac.FracStringSqrt() << " ";
	} else {
		out << _prefac.FracString() << " ";
	}
	if (not _ome_pzm.empty()) {
		out << "o(";
		for (size_t i = 0; i < _ome_pzm.size(); i++) {
			if (_ome_pzm[i] == 1) {
				out << "+";
			}
			if (_ome_pzm[i] == 0) {
				out << "0";
			}
			if (_ome_pzm[i] == -1) {
				out << "-";
			}
		}
		out << ")";
	}
	if (not _eps_pzm.empty()) {
		out << "e(";
		for (size_t i = 0; i < _eps_pzm.size(); i++) {
			if (_eps_pzm[i] == 1) {
				out << "+";
			}
			if (_eps_pzm[i] == 0) {
				out << "0";
			}
			if (_eps_pzm[i] == -1) {
				out << "-";
			}
		}
		out << ")";
	}
	if (not _chi_pzm.empty()) {
		out << "c(";
		for (size_t i = 0; i < _chi_pzm.size(); i++) {
			if (_chi_pzm[i] == 1) {
				out << "+";
			}
			if (_chi_pzm[i] == 0) {
				out << "0";
			}
			if (_chi_pzm[i] == -1) {
				out << "-";
			}
		}
		out << ")";
	}
	if (not _phi_pzm.empty()) {
		out << "p(";
		for (size_t i = 0; i < _phi_pzm.size(); i++) {
			if (_phi_pzm[i] == 1) {
				out << "+";
			}
			if (_phi_pzm[i] == 0) {
				out << "0";
			}
			if (_phi_pzm[i] == -1) {
				out << "-";
			}
		}
		out << ")";
	}
	if (_gam_s_pot) {
		if (_gam_s_pot == 1) {
			out << " gs";
		} else {
			out << " gs^" << _gam_s_pot;
		}
	}
	if (_gam_sig_pot) {
		if (_gam_sig_pot == 1) {
			out << " gsig";
		} else {
			out << " gsig^" << _gam_sig_pot;
		}
	}
	return out;
}


void TTensorTerm::shrinkVectors(const size_t& rOme,
                                const size_t& rEps,
                                const size_t& rChi,
                                const size_t& rPhi)
{
	if(rOme > _ome_pzm.size() or
	   rEps > _eps_pzm.size() or
	   rChi > _chi_pzm.size() or
	   rPhi > _phi_pzm.size())
	{
		printErr << "could not shrink arrays, they are already too small. Aborting..." << endl;
		throw;
	}
	_ome_pzm.resize(rOme);
	_eps_pzm.resize(rEps);
	_chi_pzm.resize(rChi);
	_phi_pzm.resize(rPhi);
}


void TTensorSum::Print(char flag) const
{
	for (size_t i = 0; i < _terms.size(); i++) {
		_terms[i].Print(flag);
		if (i < _terms.size() - 1) {
			cout << " ";
		}
	}
	cout << endl;
}


void TTensorSum::AddTerm(const TTensorTerm& addt) {
	_terms.push_back(addt);
}

size_t TTensorSum::SpinInnerContraction(const long& cPsiInt) {
	vector<TTensorTerm> contracted;
	for(size_t i = 0; i < _terms.size(); ++i) {
		if(_terms[i].SpinInnerContraction(cPsiInt) != 0) {
			contracted.push_back(_terms[i]);
		}
	}
	_terms = contracted;
	return _terms.size();
}

TTensorSum
TTensorSum::LSContraction(const TTensorSum& L,
                          const long& contr,
                          const long& co,
                          const long& ce,
                          const char& conType) const
{

	TTensorSum tls;
	for (size_t i = 0; i < _terms.size(); i++) {
		for (size_t j = 0; j < L.GetNterms(); j++) {
			TTensorTerm nt(_terms[i], L._terms[j], contr, co, ce, conType);
			if (nt.IsNonZero()) {
				tls.AddTerm(nt);
			}
		}
	}
	return tls;
}

TTensorSum
TTensorSum::LJContraction(const long& cChiPhi, const bool& even) {
	for (size_t i = 0; i < _terms.size(); i++) {
		_terms[i].LJContraction(cChiPhi, even);
	}
	TTensorSum tls;
	for (size_t i = 0; i < _terms.size(); i++) {
		bool foundSame = false;
		for (size_t j = 0; j < tls.GetNterms(); j++) {
			if (_terms[i].SameStructure(tls._terms[j])) {
				(tls._terms[j]).AddTwoTerms(_terms[i]);
				foundSame = true;
				break;
			}
		}
		if (not foundSame) {
			TTensorTerm nt(_terms[i]);
			tls.AddTerm(nt);
		}
	}
	tls.removeZeroTerms();
	return tls;
}

void TTensorSum::removeZeroTerms()
{
	for(size_t i= _terms.size(); i > 0; --i) {
		if(not _terms[i-1].IsNonZero()) {
			_terms.erase(_terms.begin() + (i-1));
		}
	}
}
