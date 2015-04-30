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

#include "TFracNum.h"

class ClebschGordanBox {

  public:

	//! Return field with coefficients for coupling (J1 J2|J)
	/*! The calculation is performed only once per coupling
	  as long as the ClebschGordanBox exists, so the method may be
	  called repeatedly as needed. */
	TFracNum* GetCG(const long& J, const long& J1, const long& J2);

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
	ClebschGordanBox() {
		NCG=0;
		CG=0;
	}

	TFracNum* ClebschGordan(long twoJ, long twoJ1, long twoJ2);

	long NCG;
	long* CGJ;
	long* CGJ1;
	long* CGJ2;
	TFracNum **CG;

	static ClebschGordanBox* _instance;
	static unsigned int debugCG;
	static unsigned int debugCGBox;

};

#endif
