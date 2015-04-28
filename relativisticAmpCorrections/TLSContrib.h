#ifndef TLSCONTRIB_HH
#define TLSCONTRIB_HH

#include <vector>

#include "TLSAmpl.h"
#include "TSpinWaveFunction.h"

//TODO: ask Jan how to name that
struct factors {

	factors()
		: squareOfPrefactor(),
		  exponentOfGammaS(),
		  exponentOfGammaSigma() { }

	factors(const factors& factor)
	// TODO: ask jan about "// force new power fields by multiplication " *1 ""
		: squareOfPrefactor(TFracNum::One * factor.squareOfPrefactor),
		  exponentOfGammaS(factor.exponentOfGammaS),
		  exponentOfGammaSigma(factor.exponentOfGammaSigma)
	{ }

	void swapExponents() {
		const long oldExpOfGammaS = exponentOfGammaS;
		exponentOfGammaS          = exponentOfGammaSigma;
		exponentOfGammaSigma      = oldExpOfGammaS;
	}

	TFracNum squareOfPrefactor;
	long exponentOfGammaS;
	long exponentOfGammaSigma;
};

/*!
 \class TLSContrib
 \brief Relativistic LS-coupling contributions


 \author Jan.Friedrich@ph.tum.de
 */
class TLSContrib {

  public:

	TLSContrib(const TLSContrib* b, const bool& particleExchange);
	TLSContrib(/*const*/ TLSAmpl* A, const long& delta, const TFracNum& scfac);


	void Add(TLSContrib*, bool);


	bool        SameParameter(TLSContrib* b) const { return (_J == b->_J and _L == b->_L and _S == b->_S and _cNum == b->_cNum); }
	size_t      GetNterms()                  const { return _factors.size(); }
	const bool& IsPureRelativistic()         const { return _pureRelativistic; }
	const long& GetJ()                       const { return _J; }
	const long& GetL()                       const { return _L; }
	const long& GetS()                       const { return _S; }
	const long& GetDelta()                   const { return _delta; }
	const long& GetRunningNumber()           const { return _cNum; }

	TFracNum* GetSpinCG()     { return &_SpinCG; }
	TFracNum* GetNormFactor() { return &_NormFactor; }

	const std::vector<factors>& getFactors() const { return _factors; }

	void Print()            const;
	void PrintNR()          const;
	void PrintNRG(TFracNum) const;

  private:
	long     _J;
	long     _L;
	long     _S;
	long     _cNum;
	long     _delta;
	TFracNum _SpinCG;

	TFracNum  _NormFactor;     // Square  of normalisation factor

	// TODO: check with jan how this should be named
	std::vector<factors> _factors;

	bool _pureRelativistic;

	static bool _debug;

};

#endif
