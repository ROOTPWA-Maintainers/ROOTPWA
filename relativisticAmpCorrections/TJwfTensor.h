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

#include "ClebschGordanBox.h"

class TTensorTerm {

	private:
		Int_t Rome;
		Int_t *ome_pzm;
		Int_t Reps;
		Int_t *eps_pzm;
		Int_t Rchi;
		Int_t *chi_pzm;
		Int_t Rphi;
		Int_t *phi_pzm;

		Int_t gam_s_pot;
		Int_t gam_sig_pot;

		TFracNum prefac;

	public:
		TTensorTerm () {
			Rome=0; ome_pzm=0;
			Reps=0; eps_pzm=0;
			Rchi=0; chi_pzm=0;
			Rphi=0; phi_pzm=0;
			gam_s_pot=0;
			gam_sig_pot=0;
			prefac=TFracNum::Zero;
		};

		TTensorTerm (char, Int_t, Int_t*, TFracNum*);
		TTensorTerm (TTensorTerm*, TTensorTerm*, Int_t, Int_t, Int_t, char);
		Int_t LJContraction(Int_t, Int_t);
		Int_t Multiply(char, Int_t, Int_t*, TFracNum*);
		Int_t SpinInnerContraction(Int_t);
		Int_t SameStructure(TTensorTerm*);
		Int_t AddTwoTerms(TTensorTerm*);
		Int_t IsNonZero() {
			if (prefac==TFracNum::Zero) return 0;
			else                             return 1;
		};
		Int_t Print(char);
		TFracNum GetPreFac() {return prefac;};
		Int_t GetGamS()   {return gam_s_pot;};
		Int_t GetGamSig() {return gam_sig_pot;};
		//ClassDef(TTensorTerm,1);
};

class TTensorSum {

	private:
		Int_t Nterms;
		TTensorTerm *terms;

	public:
		TTensorSum () {
			Nterms=0;
			terms=0;
		};
		Int_t AddTerm(TTensorTerm*);
		Int_t SpinInnerContraction(Int_t);
		TTensorSum* LSContraction(TTensorSum*, Int_t, Int_t, Int_t, char);
		TTensorSum* LJContraction(Int_t, Int_t);
		Int_t GetNterms() {return Nterms;};
		Int_t Print(char);
		Int_t Print() {return Print('n');}; // CINT limitation for overloading

		TTensorTerm* GetTerm(Int_t i) {return &terms[i];}
		//ClassDef(TTensorSum,1);
};
