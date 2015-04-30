#ifndef TLSAmpl_h
#define TLSAmpl_h

#include "TSpinWaveFunction.h"

/*!
 \class TLSAmpl
 \brief Relativistic LS-coupling amplitudes


 \author Jan.Friedrich@ph.tum.de
 */
class TLSAmpl {

  public:

	TLSAmpl(const long& RankS1,
	        const long& RankS2,
	        const long& RankL,
	        const long& RankJ,
	        const long& delta,
	        const long& S_L,
	        const long& cPsiInt,
	        const long& cPsiChi,
	        const long& cChiPhi,
	        const long& cPsiPhi,
	        const long& cPhiOme,
	        const long& cChiOme,
	        const long& cPhiEps,
	        const long& cChiEps,
	        const long& contractionNumber);

	const long& GetNterms()      const { return _Nterms; }
	const long& GetJ()           const { return _J; }
	const long& GetL()           const { return _L; }
	const long& GetS()           const { return _S; }
	const long& Getdelta()       const { return _delta; }
	const long& GetContraction() const { return _contractionNumber; }

	bool CheckContraction(const long& L,
	                      const long& S,
	                      const long& cPsI,
	                      const long& cPsC,
	                      const long& cCP,
	                      const long& cPsP,
	                      const long& cPO,
	                      const long& cCO,
	                      const long& cPE,
	                      const long& cCE) const;

	const TTensorTerm& GetTerm(const long& i) const { return _TSScalar.GetTerm(i); }

  private:

	long _J;
	long _L;
	long _S;
	long _delta;

	long _contractionNumber;
	long _cPsiInt;
	long _cPsiChi;
	long _cChiPhi;
	long _cPsiPhi;
	long _cPhiOme;
	long _cChiOme;
	long _cPhiEps;
	long _cChiEps;

	long _Nterms;
	TTensorSum _TSScalar;

	static unsigned int _debugLevel;

};

#endif
