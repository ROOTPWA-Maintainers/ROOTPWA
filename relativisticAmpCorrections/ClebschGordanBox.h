#ifndef ClebschGordanBox_h
#define ClebschGordanBox_h
/*!
  \class ClebschGordanBox
  \brief Handling of Clebsch-Gordan coefficients

  This Class provides the calculation and storage of Clebsch-Gordan
  coefficients. The calculation follows the derivation presented in
  J. J. Sakurai, "Modern Quantum Mechanics" (Addison-Wesley Hawaii 1984 ???),
  pp. ???-???. Since by this method, all coefficients of a multiplet (J,J1,J2)
  are obtained at once, they (actually their squares, keeping the
  original signs) are stored in a TFracNum field accessed by the method GetCG.

  In order to obtain a specific coefficient (J1,m1,J2,m2|Jm), its index in the
  field can be obtained by CGIndex(J1,m1,J2,m2). This method is defined
  globally.

  The GetCG method can be called a often as the coefficients are needed,
  since the calculation is done only once for each multiplet.

  \author Jan.Friedrich@ph.tum.de
  */

#include <iostream>
#include <map>
#include <vector>

#include "TFracNum.h"

struct quantumNumbers {

	quantumNumbers(const long& J_, const long& J1_, const long& J2_)
		: J(J_),
		  J1(J1_),
		  J2(J2_) { }

	long J;
	long J1;
	long J2;

	std::ostream& print(std::ostream& out = std::cout) const {
		out << "quantumNumbers: J=" << J << ", J1=" << J1 << ", J2=" << J2;
		return out;
	}

	quantumNumbers& operator*=(const long& factor) { J *= factor; J1 *= factor; J2 *= factor; return *this; }

};

class ClebschGordanBox {

  public:

	//! Return field with coefficients for coupling (J1 J2|J)
	/*! The calculation is performed only once per coupling
	  as long as the ClebschGordanBox exists, so the method may be
	  called repeatedly as needed. */
	const std::vector<TFracNum>& GetCG(const quantumNumbers& qN);
	const std::vector<TFracNum>& GetCG(const long& J, const long& J1, const long& J2) { return GetCG(quantumNumbers(J, J1, J2)); }

	static ClebschGordanBox* instance();

	//! Return index for a specific (J1,m1,J2,m2|Jm) Clebsch-Gordan coefficient
	/*! By construction of the field indices, the total spin J has not to be
	  specified here. */
	static long CGIndex(const long& J1,
	                    const long& m1,
	                    const long& J2,
	                    const long& m2)
	{
		return (J1 + m1) * (2 * J2 + 1) + J2 + m2;
	}

  private:

	//! Contructor of the container structure
	ClebschGordanBox()
		: _clebschGordans() { }

	const std::vector<TFracNum> ClebschGordan(const quantumNumbers& qN);

	std::map<quantumNumbers, const std::vector<TFracNum> > _clebschGordans;

	static ClebschGordanBox* _instance;
	static unsigned int _debugCG;
	static bool _debugCGBox;

};


inline
std::ostream&
operator <<(std::ostream&            out,
            const quantumNumbers&    qn)
{
	return qn.print(out);
}


inline
bool operator<(const quantumNumbers& lhs, const quantumNumbers& rhs)
{
	if(lhs.J != rhs.J) {
		return lhs.J < rhs.J;
	}
	if(lhs.J1 != rhs.J1) {
		return lhs.J1 < rhs.J1;
	}
	if(lhs.J2 != rhs.J2) {
		return lhs.J2 < rhs.J2;
	}
	return false;
}


inline
quantumNumbers operator*(const long& factor, quantumNumbers lhs)
{
	lhs *= factor;
	return lhs;
}


inline
quantumNumbers operator*(quantumNumbers lhs, const long& factor)
{
	return factor * lhs;
}


#endif
