#ifndef TSPINWAVEFUNCTION_HH
#define TSPINWAVEFUNCTION_HH

#include <vector>

#include "TJwfTensor.h"

class TSpinWaveFunction {

  public:

	TSpinWaveFunction(const size_t& J, const char& type);
	TTensorSum GetTensorSum(char name, long delta);
	TTensorSum GetSpinCoupledTensorSum(const TSpinWaveFunction& E,
	                                   const long& delta,
	                                   const long& S);

	long CheckCGFormula();

  private:
	size_t _J;

	long _max_pzm;


	long* _mi;
	long* _M;
	char _type;      // 's' Spin, 'l' Orbital Angular Momentum
	                 // 'c' Spin conjugated

	TFracNum* _coeff;

	static unsigned int _debugSpinWave;

};

#endif
