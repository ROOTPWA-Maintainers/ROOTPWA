//-----------------------------------------------------------
//
// Description:
//      mass dependent fit likelihood rank 1!
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

#ifndef MASSDEPFITLIKELI_HH
#define MASSDEPFITLIKELI_HH

#include <boost/multi_array.hpp>

#include <Math/IFunction.h>

namespace rpwa {


	class pwacompset;


	class massDepFitLikeli : public ROOT::Math::IBaseFunctionMultiDim {

	public:

		massDepFitLikeli() {}
		virtual ~massDepFitLikeli() {}

		virtual massDepFitLikeli* Clone() const;

		virtual unsigned int NDim() const;

		virtual double DoEval(const double* par) const;

		unsigned int NDataPoints() const;

		void init(pwacompset* compset,
		          const std::vector<double>& massBinCenters,
		          const boost::multi_array<std::complex<double>, 3>& spinDensityMatrices,
		          const boost::multi_array<double, 5>& spinDensityCovarianceMatrices,
		          const boost::multi_array<std::pair<size_t, size_t>, 2>& wavePairMassBinLimits,
		          bool useCovariance);

	private:

		pwacompset* _compset;

		size_t _nrMassBins;
		size_t _nrWaves;

		std::vector<double> _massBinCenters;

		boost::multi_array<std::complex<double>, 3> _spinDensityMatrices;
		boost::multi_array<double, 5> _spinDensityCovarianceMatrices;

		boost::multi_array<std::pair<size_t, size_t>, 2> _wavePairMassBinLimits;

		bool _useCovariance;

	};


} // end namespace rpwa

#endif
