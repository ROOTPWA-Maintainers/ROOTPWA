///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010-2012 Sebastian Neubert (TUM)
//    Copyright 2014-2016 Sebastian Uhl (TUM)
//
//    This file is part of ROOTPWA
//
//    ROOTPWA is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    ROOTPWA is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with ROOTPWA.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
//
// Description:
//      function to minimize for the resonance fit
//      - at the moment only rank 1 is handled
//      - two different methods to fit are implemented
//        * a fit to the spin-density matrix
//        * a fit to the production amplitudes
//      - several methods to take errors/the covariance matrix into account
//        * single part of the complex number
//        * covariance for a complex number
//        * full covariance matrix while fitting to production amplitudes
//
//-------------------------------------------------------------------------


#ifndef RESONANCEFIT_FUNCTION_HH
#define RESONANCEFIT_FUNCTION_HH

#include <boost/multi_array.hpp>

#include <TMatrixT.h>
#include <TVectorT.h>

#include "forward.h"

namespace rpwa {

	namespace resonanceFit {

		class cache;
		class parameters;

		class function {

		public:

			enum useCovarianceMatrix {
				useDiagnalElementsOnly,
				useComplexDiagnalElementsOnly,
				useFullCovarianceMatrix,
				useCovarianceMatrixDefault
			};

			function(const bool fitProductionAmplitudes,
			         const rpwa::resonanceFit::function::useCovarianceMatrix useCovariance);
			~function() {}

			bool init(const rpwa::resonanceFit::modelConstPtr& fitModel,
			          const std::vector<size_t>& nrMassBins,
			          const boost::multi_array<double, 2>& massBinCenters,
			          const boost::multi_array<std::complex<double>, 3>& productionAmplitudes,
			          const boost::multi_array<TMatrixT<double>, 2>& productionAmplitudesCovariance,
			          const boost::multi_array<std::complex<double>, 4>& spinDensityMatrices,
			          const boost::multi_array<TMatrixT<double>, 2>& spinDensityCovarianceMatrices,
			          const boost::multi_array<std::pair<size_t, size_t>, 3>& wavePairMassBinLimits);

			size_t getNrParameters() const;
			size_t getNrDataPoints() const;

			double chiSquare(const std::vector<double>& par) const;
			double chiSquare(const double* par) const;
			double chiSquare(const rpwa::resonanceFit::parameters& fitParameters,
			                 rpwa::resonanceFit::cache& cache) const;

			double logLikelihood(const std::vector<double>& par) const;
			double logLikelihood(const double* par) const;
			double logLikelihood(const rpwa::resonanceFit::parameters& fitParameters,
			                     rpwa::resonanceFit::cache& cache) const;

			double logPriorLikelihood(const std::vector<double>& par) const;
			double logPriorLikelihood(const double* par) const;
			double logPriorLikelihood(const rpwa::resonanceFit::parameters& fitParameters) const;

		private:

			double chiSquareProductionAmplitudes(const rpwa::resonanceFit::parameters& fitParameters,
			                                     rpwa::resonanceFit::cache& cache) const;
			double chiSquareSpinDensityMatrix(const rpwa::resonanceFit::parameters& fitParameters,
			                                  rpwa::resonanceFit::cache& cache) const;

			rpwa::resonanceFit::modelConstPtr _fitModel;

			size_t _nrBins;
			size_t _maxMassBins;
			size_t _nrWaves;

			std::vector<size_t> _idxMassMin;
			std::vector<size_t> _idxMassMax;

			std::vector<size_t> _nrMassBins;
			boost::multi_array<double, 2> _massBinCenters;

			boost::multi_array<std::complex<double>, 3> _productionAmplitudes;
			boost::multi_array<TMatrixT<double>, 2> _productionAmplitudesCovMatInv;

			boost::multi_array<std::complex<double>, 4> _spinDensityMatrices;
			boost::multi_array<TMatrixT<double>, 2> _spinDensityMatricesCovMatInv;
			boost::multi_array<double, 6> _spinDensityMatricesCovMatInvArray;

			boost::multi_array<std::pair<size_t, size_t>, 3> _wavePairMassBinLimits;

			size_t _idxAnchorWave;

			const bool _fitProductionAmplitudes;
			const rpwa::resonanceFit::function::useCovarianceMatrix _useCovariance;

		};

	} // end namespace resonanceFit

} // end namespace rpwa

#endif // RESONANCEFIT_FUNCTION_HH
