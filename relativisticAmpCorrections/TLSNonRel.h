#ifndef TLSNONREL_HH
#define TLSNONREL_HH

#include "TLSContrib.h"

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
