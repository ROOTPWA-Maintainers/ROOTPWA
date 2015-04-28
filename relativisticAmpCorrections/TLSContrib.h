#ifndef TLSCONTRIB_HH
#define TLSCONTRIB_HH

#include "TLSAmpl.h"
#include "TSpinWaveFunction.h"

/*!
 \class TLSContrib
 \brief Relativistic LS-coupling contributions


 \author Jan.Friedrich@ph.tum.de
 */
class TLSContrib {

  public:

	TLSContrib()
		: _J(0),
		  _L(0),
		  _S(0),
		  _delta(0),
		  _Nterms(0),
		  _termFracNum(0),
		  _termg1pot(0),
		  _termg2pot(0) { }

	TLSContrib(TLSContrib* b, bool particleExchange);
	TLSContrib(TLSAmpl* A, long delta, TFracNum scfac);


	long Add(TLSContrib*, bool);


	bool SameParameter(TLSContrib* b) { return (_J == b->_J && _L == b->_L && _S == b->_S && _cNum == b->_cNum); }
	const long& GetNterms() { return _Nterms; }
	const bool& IsPureRelativistic() { return _pureRelativistic; }
	const long& GetJ() { return _J; }
	const long& GetL() { return _L; }
	const long& GetS() { return _S; }
	const long& GetDelta() { return _delta; }
	const long& GetRunningNumber() { return _cNum; }

	TFracNum* GetSpinCG() { return &_SpinCG; }
	TFracNum* GetNormFactor() { return &_NormFactor; }

	TFracNum* GetTermFracNum() { return _termFracNum; }
	long* GetTermg1pot() { return _termg1pot; } // exponent of gamma_s
	long* GetTermg2pot() { return _termg2pot; } // exponent of gamma_s

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

	long      _Nterms;
	TFracNum  _NormFactor;     // Square  of normalisation factor

// TODO: struct these?
	TFracNum* _termFracNum;   // Squares of prefactors
	long* _termg1pot;     // exponent of gamma_s
	long* _termg2pot;     // exponent of gamma_sigma

	bool _pureRelativistic;

	static bool _debug;

};

#endif
