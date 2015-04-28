#ifndef TLSNONREL_HH
#define TLSNONREL_HH

#include <vector>

#include "TLSContrib.h"

/*!
  \class TLSNonRel
  \brief Non-relativistic LS-coupling contributions


  \author Jan.Friedrich@ph.tum.de
  */
class TLSNonRel {

	public:
		TLSNonRel(TLSContrib *C);

		long CheckJLS(TLSContrib *C) {
			return (C->GetJ()==_J && C->GetL()==_L && C->GetS()==_S) ? 1 : 0;
		};
		void Add(TLSContrib *C);
		long GetL(){return _L;};
		long GetS(){return _S;};
		long Print();
		long PrintG();

	private:
		long _J;
		long _L;
		long _S;
		std::vector<TLSContrib*> _RelLS;
		TFracNum _GnrPrefac;

};

#endif
