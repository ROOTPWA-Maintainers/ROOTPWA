#ifndef TJWFTENSOR_HH
#define TJWFTENSOR_HH

#include <iostream>
#include <vector>

#include "TFracNum.h"

class TTensorTerm {

  public:

	TTensorTerm()
	: _ome_pzm(),
	  _eps_pzm(),
	  _chi_pzm(),
	  _phi_pzm(),
	  _gam_s_pot(0),
	  _gam_sig_pot(0),
	  _prefac(TFracNum::Zero) {
	}

	// TODO: optimize this call
	TTensorTerm(char name,
	            const std::vector<long>& pzm_field,
	            const TFracNum& prefac);

	TTensorTerm(const TTensorTerm& S,
	            const TTensorTerm& L,
	            long contractions,
	            long o_share,
	            long e_share,
	            char con_type);

	long LJContraction(long ncon, long even);
	// TODO: optimize this call
	void Multiply(char name, const std::vector<long>& pzm_field, const TFracNum& prefac);
	long SpinInnerContraction(const long& cPsiInt);
	bool SameStructure(const TTensorTerm& rhs) const;
	bool AddTwoTerms(const TTensorTerm& rhs);
	long IsNonZero() const { return (_prefac == TFracNum::Zero) ? 0 : 1; }

	std::ostream& Print(const char& flag = 'n', std::ostream& out = std::cout) const;

	const TFracNum& GetPreFac() const { return _prefac; }
	const long&     GetGamS()   const { return _gam_s_pot; }
	const long&     GetGamSig() const { return _gam_sig_pot; }

private:

	void shrinkVectors(const size_t& rOme,
	                   const size_t& rEps,
	                   const size_t& rChi,
	                   const size_t& rPhi);

	std::vector<long> _ome_pzm;
	std::vector<long> _eps_pzm;
	std::vector<long> _chi_pzm;
	std::vector<long> _phi_pzm;

	long _gam_s_pot;
	long _gam_sig_pot;

	TFracNum _prefac;

};

inline
std::ostream&
operator <<(std::ostream&            out,
            const TTensorTerm&       tensorTerm)
{
	return tensorTerm.Print('n', out);
}


class TTensorSum {

private:
	long Nterms;
	TTensorTerm *terms;

public:
	TTensorSum() {
		Nterms = 0;
		terms = 0;
	}
	;
	long AddTerm(TTensorTerm*);
	long SpinInnerContraction(long);
	TTensorSum* LSContraction(TTensorSum*, long, long, long, char);
	TTensorSum* LJContraction(long, long);
	long GetNterms() {
		return Nterms;
	}
	;
	long Print(char = 'n');

	TTensorTerm* GetTerm(long i) {
		return &terms[i];
	}

};

#endif
