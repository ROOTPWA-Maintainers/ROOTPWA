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

	private:
		long J;
		long L;
		long S;
		long cNum;
		long delta;
		TFracNum SpinCG;

		long Nterms;
		TFracNum NormFactor;     // Square  of normalisation factor
		TFracNum *termFracNum;   // Squares of prefactors
		long    *termg1pot;     // exponent of gamma_s
		long    *termg2pot;     // exponent of gamma_sigma

		bool PureRelativistic;

	public:

		TLSContrib() {
			J=0; L=0; S=0; delta=0;
			Nterms=0;
			termFracNum=0;
			termg1pot=0;
			termg2pot=0;
		};

		TLSContrib(TLSContrib *, bool);
		TLSContrib(TLSAmpl*, long, TFracNum);

		bool SameParameter(TLSContrib* b) {
			if (J==b->J && L==b->L && S==b->S && cNum==b->cNum)
				return true;
			return false;
		};
		long GetNterms() {return Nterms;};
		long Add(TLSContrib*, bool);
		long Print();
		long PrintNR();
		long PrintNRG(TFracNum);
		bool IsPureRelativistic() {return PureRelativistic;};
		long GetJ() {return J;};
		long GetL() {return L;};
		long GetS() {return S;};
		long GetDelta() {return delta;};
		long GetRunningNumber() {return cNum;};
		TFracNum* GetSpinCG() {return &SpinCG;};
		TFracNum* GetNormFactor() {return &NormFactor;};

		TFracNum* GetTermFracNum() {return termFracNum;};
		long*    GetTermg1pot()   {return termg1pot;}; // exponent of gamma_s
		long*    GetTermg2pot()   {return termg2pot;}; // exponent of gamma_s

};

#endif
