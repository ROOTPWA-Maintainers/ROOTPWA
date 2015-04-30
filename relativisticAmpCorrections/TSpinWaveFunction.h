#ifndef TSPINWAVEFUNCTION_HH
#define TSPINWAVEFUNCTION_HH

#include <vector>

#include "TJwfTensor.h"

class TSpinWaveFunction {

  public:

	TSpinWaveFunction(const size_t& J, const char& type);
	TTensorSum GetTensorSum(const char& name, const long& delta) const;
	TTensorSum GetSpinCoupledTensorSum(const TSpinWaveFunction& E,
	                                   const long& delta,
	                                   const long& S) const;

#if(0)
	long CheckCGFormula() const;
#endif

  private:
	size_t _J;
	std::vector<long> _mi;
	std::vector<std::pair<long, TFracNum> > _M_and_coeff;

	char _type;      // 's' Spin, 'l' Orbital Angular Momentum
	                 // 'c' Spin conjugated

	static unsigned int _debugSpinWave;

};

#endif
