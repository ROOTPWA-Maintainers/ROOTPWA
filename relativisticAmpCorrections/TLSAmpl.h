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

#endif
