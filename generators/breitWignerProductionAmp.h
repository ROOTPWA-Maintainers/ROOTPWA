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

#ifndef TBWPRODUCTIONAMP_HH
#define TBWPRODUCTIONAMP_HH

#include "productionAmp.h"

#include <complex>

namespace rpwa {

	class breitWignerProductionAmp : public productionAmp {

	  public:

			breitWignerProductionAmp(double mass,
			                         double width,
			                         std::complex<double> coupling=std::complex<double>(1,0));
			virtual ~breitWignerProductionAmp() { }

			virtual std::complex<double> amp(double mass); //> simple one parameter case

	  private:

			double _mass;
			double _m2;
			double _width;
			double _mw;
			std::complex<double> _coupling;

	};

}

#endif

