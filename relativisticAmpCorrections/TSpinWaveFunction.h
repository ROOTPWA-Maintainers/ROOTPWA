#ifndef TSPINWAVEFUNCTION_HH
#define TSPINWAVEFUNCTION_HH

#include "TJwfTensor.h"

class TSpinWaveFunction {

  public:

#if(0)
	TSpinWaveFunction()
		: _J(0),
		  _max_pzm(1),
		  _pot3(0),
		  _mi(0),
		  _M(0) { }
#endif

	TSpinWaveFunction(long J, char type);
	TTensorSum GetTensorSum(char name, long delta);
	TTensorSum GetSpinCoupledTensorSum(const TSpinWaveFunction& E,
	                                   const long& delta,
	                                   const long& S);

	long CheckCGFormula();

  private:
	long _J;
	long _max_pzm;
	long* _pot3;
	long* _mi;
	long* _M;
	char _type;      // 's' Spin, 'l' Orbital Angular Momentum
	                 // 'c' Spin conjugated

	TFracNum* _coeff;

	static unsigned int _debugSpinWave;

};

#endif
