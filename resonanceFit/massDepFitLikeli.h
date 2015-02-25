///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010-2012 Sebastian Neubert (TUM)
//    Copyright 2014,2015 Sebastian Uhl (TUM)
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
//      likelihood for resonance fit
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


#ifndef MASSDEPFITLIKELI_HH
#define MASSDEPFITLIKELI_HH

#include <boost/multi_array.hpp>

#include <Math/IFunction.h>
#include <TMatrixT.h>

namespace rpwa {

	namespace massDepFit {

		class model;
		class parameters;

		class likelihood : public ROOT::Math::IBaseFunctionMultiDim {

		public:

			enum useCovarianceMatrix {
				useDiagnalElementsOnly,
				useComplexDiagnalElementsOnly,
				useFullCovarianceMatrix,
				useCovarianceMatrixDefault
			};

			likelihood(const bool fitProductionAmplitudes,
			           const rpwa::massDepFit::likelihood::useCovarianceMatrix useCovariance);
			virtual ~likelihood() {}

			virtual likelihood* Clone() const;

			virtual unsigned int NDim() const;

			virtual double DoEval(const double* par) const;
			virtual double DoEval(const rpwa::massDepFit::parameters& fitParameters) const;

			unsigned int NDataPoints() const;

			bool init(rpwa::massDepFit::model* compset,
			          const std::vector<double>& massBinCenters,
			          const boost::multi_array<std::complex<double>, 3>& productionAmplitudes,
			          const boost::multi_array<double, 6>& productionAmplitudesCovariance,
			          const boost::multi_array<std::complex<double>, 4>& spinDensityMatrices,
			          const boost::multi_array<double, 6>& spinDensityCovarianceMatrices,
			          const boost::multi_array<std::pair<size_t, size_t>, 2>& wavePairMassBinLimits);

		private:

			double DoEvalProductionAmplitudes(const rpwa::massDepFit::parameters& fitParameters) const;
			double DoEvalSpinDensityMatrix(const rpwa::massDepFit::parameters& fitParameters) const;

			rpwa::massDepFit::model* _compset;

			size_t _nrBins;
			size_t _nrMassBins;
			size_t _nrWaves;

			size_t _idxMassMin;
			size_t _idxMassMax;

			std::vector<double> _massBinCenters;

			boost::multi_array<std::complex<double>, 3> _productionAmplitudes;
			boost::multi_array<double, 6> _productionAmplitudesCovariance;
			boost::multi_array<TMatrixT<double>, 2> _productionAmplitudesCovMatInv;

			boost::multi_array<std::complex<double>, 4> _spinDensityMatrices;
			boost::multi_array<double, 6> _spinDensityCovarianceMatrices;

			boost::multi_array<std::pair<size_t, size_t>, 2> _wavePairMassBinLimits;

			size_t _idxAnchorWave;

			const bool _fitProductionAmplitudes;
			rpwa::massDepFit::likelihood::useCovarianceMatrix _useCovariance;

		};

	} // end namespace massDepFit

} // end namespace rpwa

#endif
