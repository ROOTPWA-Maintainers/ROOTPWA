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

			function(const rpwa::resonanceFit::dataConstPtr& fitData,
			         const rpwa::resonanceFit::modelConstPtr& fitModel,
			         const bool useProductionAmplitudes);
			~function() {}

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

			const rpwa::resonanceFit::dataConstPtr _fitData;
			const rpwa::resonanceFit::modelConstPtr _fitModel;

			const size_t _nrBins;
			const size_t _maxNrWaves;
			const size_t _maxNrMassBins;

			std::vector<size_t> _idxMassMin;
			std::vector<size_t> _idxMassMax;

			const bool _useProductionAmplitudes;
			const rpwa::resonanceFit::function::useCovarianceMatrix _useCovariance;

		};

	} // end namespace resonanceFit

} // end namespace rpwa

#endif // RESONANCEFIT_FUNCTION_HH
