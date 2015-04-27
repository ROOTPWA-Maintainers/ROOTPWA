#ifndef TFhh_h
#define TFhh_h
/*!
  \class TFhh
  \brief Relativistic helicity-coupling structure

  The class describes relativistic helicity-coupling amplitudes, as
  described in the paper [PRD78(2008)], equation (6.6). It collects a field
  of terms TLSContrib with different LS. It also provides the relation to the
  non-relativistic Zeemach amplitudes.

  \author Jan.Friedrich@ph.tum.de
  */

#include <string>
#include <vector>

#include "TLSAmpl.h"

class TFhh {

  public:

	TFhh()
	: _J(0),
	  _lambda(0),
	  _nu(0),
	  _nTerms(0),
	  _LSt(0) { }

	TFhh(const long& J_,
	     const long& S1,
	     const long& S2,
	     const long& lambda_,
	     const long& nu_,
	     const std::vector<TLSAmpl*>& LSampl,
	     const bool& even_contr_);

	TFhh(TFhh*, char);
	TFhh(TFhh*, TFhh*);

	const long&        GetNterms()          const { return _nTerms; }
	bool               IsNuNu()             const { return (_lambda ==  _nu); }
	bool               IsNuMinusNu()        const { return (_lambda == -_nu); }
	const long&        GetLambda()          const { return _lambda; }
	const long&        GetNu()              const { return _nu;     }
	const long&        GetJ()               const { return _J;      }
	const bool&        GetEvenContraction() const { return _evenContraction;}
	TLSContrib**       GetLStPtr()                { return _LSt; }
	const std::string& GetName()            const { return _name_str; }

	void NonRelLimit();
	void PrintNRG() const;
	void Print()    const;

  private:

	std::string _name_str;
	long _J;
	long _lambda;
	long _nu;
	bool _evenContraction;
	long _nTerms;
	TLSContrib* *_LSt;
	long _NNRterms;
	TLSNonRel* *_NRLSt;

	static unsigned int _debugLevel;

};


#endif
