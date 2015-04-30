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

	TLSAmpl(long RankS1,
	        long RankS2,
	        long RankL,
	        long RankJ,
	        long delta,
	        long S_L,
	        long cPsiInt,
	        long cPsiChi,
	        long cChiPhi,
	        long cPsiPhi,
	        long cPhiOme,
	        long cChiOme,
	        long cPhiEps,
	        long cChiEps,
	        long contractionNumber);

	const long& GetNterms()      const { return _Nterms; }
	const long& GetJ()           const { return _J; }
	const long& GetL()           const { return _L; }
	const long& GetS()           const { return _S; }
	const long& Getdelta()       const { return _delta; }
	const long& GetContraction() const { return _contractionNumber; }

	bool CheckContraction(long L,
	                      long S,
	                      long cPsI,
	                      long cPsC,
	                      long cCP,
	                      long cPsP,
	                      long cPO,
	                      long cCO,
	                      long cPE,
	                      long cCE) const;

	const TTensorTerm& GetTerm(long i) const { return _TSScalar->GetTerm(i); }

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
	TTensorSum* _TSScalar;

	static unsigned int _debugLevel;

};

#endif
