#ifndef TJSS_HH
#define TJSS_HH

#include "TLSAmpl.h"
#include "TFhh.h"


class TJSS {

  private:

	long JMother;
	long etaJ;
	long SDecay1;
	long eta1;
	long SDecay2;
	long eta2;

	long NLSAmpl;
	TLSAmpl* *LSAmplitudes;

	long NFhhAmpl;
	TFhh* *FhhAmpl;

	long NFhhIdAmpl;
	TFhh* *FhhIdAmpl;

	static unsigned int _debugLevel;

  public:

	TJSS()
		: JMother(0),
		  etaJ(1),
		  SDecay1(0),
		  eta1(1),
		  SDecay2(0),
		  eta2(1),
		  NLSAmpl(0),
		  LSAmplitudes(0),
		  NFhhAmpl(0),
		  FhhAmpl(0),
		  NFhhIdAmpl(0),
		  FhhIdAmpl(0) { }

	TJSS(long J,  long eJ,
	     long S1, long e1,
	     long S2, long e2)
		: JMother(J),
		  etaJ(eJ),
		  SDecay1(S1),
		  eta1(e1),
		  SDecay2(S2),
		  eta2(e2),
		  NLSAmpl(0),
		  LSAmplitudes(0),
		  NFhhAmpl(0),
		  FhhAmpl(0),
		  NFhhIdAmpl(0),
		  FhhIdAmpl(0) { }

	long CalcAmpl();

};

#endif
