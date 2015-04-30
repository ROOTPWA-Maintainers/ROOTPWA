#ifndef TJSS_HH
#define TJSS_HH

#include <string>
#include <vector>
#include <utility>

#include "TLSAmpl.h"
#include "TFhh.h"


class TJSS {

  public:

	TJSS(long J,
	     long eJ,
	     long S1,
	     long e1,
	     long S2,
	     long e2)
		: _JMother(J),
		  _etaJ(eJ),
		  _SDecay1(S1),
		  _eta1(e1),
		  _SDecay2(S2),
		  _eta2(e2),
		  _LSAmplitudes(),
		  _FhhAmpl(),
		  _FhhIdAmpl() { }

	void CalcAmpl();

  private:

	bool checkSelectionRules(const long& PsiInternal,
	                         const long& cPsiChi,
	                         const long& cChiPhi,
	                         const long& cPsiPhi,
	                         const long& IndexContractions,
	                         const long& even_contraction,
	                         const long& L) const;

	bool getContractionRank(const long& PsiInternal,
	                        const long& cPsiChi,
	                        const long& cPsiPhi,
	                        const bool& evenContraction,
	                        std::pair<long, long>& rPair) const;

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

	std::vector<TLSAmpl*> _LSAmplitudes;
	std::vector<TFhh*>    _FhhAmpl;
	std::vector<TFhh*>    _FhhIdAmpl;

	static unsigned int _debugLevel;

};

#endif
