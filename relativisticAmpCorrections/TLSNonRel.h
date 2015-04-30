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
		TLSNonRel(TLSContrib* C);

		bool CheckJLS(const TLSContrib* C) const {
			return (C->GetJ() == _J and C->GetL() == _L and C->GetS() == _S);
		}
		void Add(TLSContrib* C);
		const long& GetL()   const { return _L; }
		const long& GetS()   const { return _S; }
		void Print()  const;
		void PrintG() const;

	private:
		long _J;
		long _L;
		long _S;
		std::vector<TLSContrib*> _RelLS;
		TFracNum _GnrPrefac;

};

#endif
