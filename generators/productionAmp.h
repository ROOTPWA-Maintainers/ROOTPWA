//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Parameterization of Production amplitude
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

#ifndef TPRODUCTIONAMP_HH
#define TPRODUCTIONAMP_HH

#include <complex>


namespace rpwa {

	class productionAmp {

	  public:

			productionAmp(){};
			productionAmp(const std::complex<double>& amplitude);
			virtual ~productionAmp() { }

			virtual std::complex<double> amp(double mass) { return _amp; } //> simple one parameter case

	  private:

			std::complex<double> _amp;


	};

}

#endif

