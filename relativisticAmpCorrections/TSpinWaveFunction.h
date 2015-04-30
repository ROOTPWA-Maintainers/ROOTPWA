#ifndef TSPINWAVEFUNCTION_HH
#define TSPINWAVEFUNCTION_HH

#include "TJwfTensor.h"

class TSpinWaveFunction {

private:
	long J;
	long max_pzm;
	long *pot3;
	long *mi;
	long *M;
	char type;      // 's' Spin, 'l' Orbital Angular Momentum
	// 'c' Spin conjugated
	TFracNum *coeff;

public:
	TSpinWaveFunction() {
		J = 0;
		max_pzm = 1;
		pot3 = 0;
		mi = 0;
		M = 0;
	}

	TSpinWaveFunction(long J, char type);
	TTensorSum* GetTensorSum(char name, long delta);
	TTensorSum* GetSpinCoupledTensorSum(TSpinWaveFunction*, long, long);

	long CheckCGFormula();

};

#endif
