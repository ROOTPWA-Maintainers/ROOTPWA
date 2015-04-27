#ifndef TJSS_HH
#define TJSS_HH

#include <string>

#include "TLSAmpl.h"
#include "TFhh.h"


class TJSS {


  public:

	TJSS()
		: _JMother(0),
		  _etaJ(1),
		  _SDecay1(0),
		  _eta1(1),
		  _SDecay2(0),
		  _eta2(1),
		  _NLSAmpl(0),
		  _LSAmplitudes(0),
		  _NFhhAmpl(0),
		  _FhhAmpl(0),
		  _NFhhIdAmpl(0),
		  _FhhIdAmpl(0) { }

	TJSS(long J,  long eJ,
	     long S1, long e1,
	     long S2, long e2)
		: _JMother(J),
		  _etaJ(eJ),
		  _SDecay1(S1),
		  _eta1(e1),
		  _SDecay2(S2),
		  _eta2(e2),
		  _NLSAmpl(0),
		  _LSAmplitudes(0),
		  _NFhhAmpl(0),
		  _FhhAmpl(0),
		  _NFhhIdAmpl(0),
		  _FhhIdAmpl(0) { }

	void CalcAmpl();

  private:

	std::string getContractionPattern(const long& cPsiPhi,
	                                  const long& cPhiOmega,
	                                  const long& cPhiEps,
	                                  const long& PsiInternal,
	                                  const long& cPsiChi,
	                                  const long& cChiOmega,
	                                  const long& cChiEps,
	                                  const long& cChiPhi,
	                                  const long& L) const;

	long _JMother;
	long _etaJ;
	long _SDecay1;
	long _eta1;
	long _SDecay2;
	long _eta2;

	long _NLSAmpl;
	TLSAmpl* *_LSAmplitudes;

	long _NFhhAmpl;
	TFhh* *_FhhAmpl;

	long _NFhhIdAmpl;
	TFhh* *_FhhIdAmpl;

	static unsigned int _debugLevel;

};

#endif
