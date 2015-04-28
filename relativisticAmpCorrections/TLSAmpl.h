#ifndef TLSAmpl_h
#define TLSAmpl_h

#include "TSpinWaveFunction.h"

/*!
  \class TLSAmpl
  \brief Relativistic LS-coupling amplitudes


  \author Jan.Friedrich@ph.tum.de
  */
class TLSAmpl {

	private:
		long J;
		long L;
		long S;
		long delta;

		long ContractionNumber;
		long cPsI; long cPsC;
		long cCP;  long cPsP;
		long cPO; long cCO;
		long cPE; long cCE;

		long Nterms;
		TTensorSum *TSScalar;

	public:

		TLSAmpl() {
			J=0; L=0; S=0; delta=0;
			Nterms=0;
		};

		TLSAmpl(long RankS1, long RankS2,
				long RankL,  long RankJ,
				long delta,  long S_L,
				long cPsiInt, long cPsiChi,
				long cChiPhi, long cPsiPhi,
				long cPhiOme, long cChiOme,
				long cPhiEps, long cChiEps,
				long cNum);

		long GetNterms() {return Nterms;};

		long GetJ()     {return J;};
		long GetL()     {return L;};
		long GetS()     {return S;};
		long Getdelta() {return delta;};
		long GetContraction() {return ContractionNumber;};
		bool CheckContraction(long L_, long S_, long cPsI_, long cPsC_,
				long cCP_, long cPsP_, long cPO_, long cCO_,
				long cPE_, long cCE_) {
			if ( L!=L_ || S!=S_ || cPsI!=cPsI_ || cPsC!=cPsC_ || cCP!=cCP_ ||
					cPsP!=cPsP_ || cPO!=cPO_ || cCO!=cCO_ || cPE!=cPE_ || cCE!=cCE_ )
				return false;
			return true;
		}
		TTensorTerm *GetTerm(long i) {return TSScalar->GetTerm(i);};

};

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


/*!
  \class TLSNonRel
  \brief Non-relativistic LS-coupling contributions


  \author Jan.Friedrich@ph.tum.de
  */
class TLSNonRel {

	private:
		long J;
		long L;
		long S;
		long Nterms;
		TLSContrib* *RelLS;
		TFracNum GnrPrefac;

	public:
		TLSNonRel(TLSContrib *C);

		long CheckJLS(TLSContrib *C) {
			return (C->GetJ()==J && C->GetL()==L && C->GetS()==S) ? 1 : 0;
		};
		long Add(TLSContrib *C);
		long GetL(){return L;};
		long GetS(){return S;};
		long Print();
		long PrintG();

};
#endif
