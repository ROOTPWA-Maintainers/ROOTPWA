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
//
// Uncomment the following line
// if you want to work in CINT (root.cern.ch)
//
//#define __JCINT__
#ifndef __JCINT__
#define Int_t    long long
#define Double_t double
#define Bool_t   bool
#endif

#include "TFracNum.h"

class ClebschGordanBox {
	private:
		Int_t NCG;
		Int_t *CGJ;
		Int_t *CGJ1;
		Int_t *CGJ2;
		TFracNum **CG;
		TFracNum* ClebschGordan(Int_t twoJ, Int_t twoJ1, Int_t twoJ2);
	public:
		//! Contructor of the container structure
		ClebschGordanBox() {
			NCG=0;
			CG=0;
		}
		//! Return field with coefficients for coupling (J1 J2|J)
		/*! The calculation is performed only once per coupling 
		  as long as the ClebschGordanBox exists, so the method may be
		  called repeatedly as needed. */
		TFracNum* GetCG(Int_t J, Int_t J1, Int_t J2);
		//ClassDef(ClebschGordanBox,1);
};

//! Return index for a specific (J1,m1,J2,m2|Jm) Clebsch-Gordan coefficient
/*! By construction of the field indices, the total spin J has not to be
  specified here. */
Int_t CGIndex(Int_t J1, Int_t m1, Int_t J2, Int_t m2);

#endif
