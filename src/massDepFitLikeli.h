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

#include <boost/numeric/ublas/matrix.hpp>

#include <Math/IFunction.h>

namespace rpwa {


	typedef boost::numeric::ublas::matrix<std::complex<double> > spinDensityMatrixType;
	typedef boost::numeric::ublas::matrix<double> complexCovarianceMatrixType;
	typedef boost::numeric::ublas::matrix<complexCovarianceMatrixType> spinDensityCovarianceMatrixType;


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
		          const std::vector<spinDensityMatrixType>& spinDensityMatrices,
		          const std::vector<spinDensityCovarianceMatrixType>& spinDensityCovarianceMatrices,
		          const std::vector<std::vector<std::pair<size_t, size_t> > >& wavePairMassBinLimits,
		          bool useCovariance);

	private:

		pwacompset* _compset;

		std::vector<double> _massBinCenters;

		std::vector<spinDensityMatrixType> _spinDensityMatrices;
		std::vector<spinDensityCovarianceMatrixType> _spinDensityCovarianceMatrices;

		std::vector<std::vector<std::pair<size_t, size_t> > > _wavePairMassBinLimits;

		bool _useCovariance;

	};


} // end namespace rpwa

#endif
