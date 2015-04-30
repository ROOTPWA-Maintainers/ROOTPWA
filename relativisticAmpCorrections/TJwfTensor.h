#ifndef TJWFTENSOR_HH
#define TJWFTENSOR_HH

#include "ClebschGordanBox.h"

class TTensorTerm {

  public:

	TTensorTerm()
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
	  _prefac(TFracNum::Zero) {
	}

	TTensorTerm(char name,
	            long RJ,
	            long* pzm_field,
	            const TFracNum& prefac);

	TTensorTerm(const TTensorTerm& S,
	            const TTensorTerm& L,
	            long contractions,
	            long o_share,
	            long e_share,
	            char con_type);

	long LJContraction(long ncon, long even);
	long Multiply(char name, long RJ, long* pzm_field, const TFracNum& prefac);
	long SpinInnerContraction(long cPsiInt);
	bool SameStructure(const TTensorTerm& rhs) const;
	bool AddTwoTerms(const TTensorTerm& rhs);
	long IsNonZero() const { return (_prefac == TFracNum::Zero) ? 0 : 1; }

	long Print(char flag) const;

	const TFracNum& GetPreFac() const { return _prefac; }
	const long&     GetGamS()   const { return _gam_s_pot; }
	const long&     GetGamSig() const { return _gam_sig_pot; }

private:

	long _Rome;
	long* _ome_pzm;
	long _Reps;
	long* _eps_pzm;
	long _Rchi;
	long* _chi_pzm;
	long _Rphi;
	long* _phi_pzm;

	long _gam_s_pot;
	long _gam_sig_pot;

	TFracNum _prefac;

};

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
	long Print(char);
	long Print() {
		return Print('n');
	}
	; // CINT limitation for overloading

	TTensorTerm* GetTerm(long i) {
		return &terms[i];
	}

};

#endif
