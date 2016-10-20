///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2016 Sebastian Uhl (TUM)
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
//      storage for the data required by the resonance fit
//      - vectors of production amplitudes
//      - spin-density matrices
//      - the corresponding covariances
//      - phase-space integrals of waves
//
//-------------------------------------------------------------------------


#ifndef RESONANCEFIT_DATA_HH
#define RESONANCEFIT_DATA_HH

#include <boost/multi_array.hpp>

#include <TMatrixT.h>

#include "function.h"

namespace rpwa {

	namespace resonanceFit {

		class data {

		public:

			data(const std::vector<size_t>& nrMassBins,
			     const boost::multi_array<double, 2>& massBinCenters,
			     const boost::multi_array<std::complex<double>, 3>& productionAmplitudes,
			     const boost::multi_array<TMatrixT<double>, 2>& productionAmplitudesCovMatInv,
			     const boost::multi_array<std::complex<double>, 4>& spinDensityMatrixElements,
			     const boost::multi_array<TMatrixT<double>, 2>& spinDensityMatrixElementsCovMatInv,
			     const rpwa::resonanceFit::function::useCovarianceMatrix useCovariance,
			     const boost::multi_array<std::pair<double, double>, 3>& plottingIntensities,
			     const boost::multi_array<std::pair<double, double>, 4>& plottingSpinDensityMatrixElementsReal,
			     const boost::multi_array<std::pair<double, double>, 4>& plottingSpinDensityMatrixElementsImag,
			     const boost::multi_array<std::pair<double, double>, 4>& plottingPhases,
			     const boost::multi_array<std::pair<double, double>, 3>& sysPlottingIntensities,
			     const boost::multi_array<std::pair<double, double>, 4>& sysPlottingSpinDensityMatrixElementsReal,
			     const boost::multi_array<std::pair<double, double>, 4>& sysPlottingSpinDensityMatrixElementsImag,
			     const boost::multi_array<std::pair<double, double>, 4>& sysPlottingPhases);
			~data() {}

			size_t nrBins() const { return _nrMassBins.size(); }
			size_t maxMassBins() const { return *(std::max_element(_nrMassBins.begin(), _nrMassBins.end())); }
			size_t nrWaves() const { return *(_productionAmplitudes.shape()+2); }

			const std::vector<size_t>& nrMassBins() const { return _nrMassBins; }
			const boost::multi_array<double, 2>& massBinCenters() const { return _massBinCenters; }

			const boost::multi_array<std::complex<double>, 3>& productionAmplitudes() const { return _productionAmplitudes; }
			const boost::multi_array<TMatrixT<double>, 2>& productionAmplitudesCovMatInv() const { return _productionAmplitudesCovMatInv; }

			const boost::multi_array<std::complex<double>, 4>& spinDensityMatrixElements() const { return _spinDensityMatrixElements; }
			const boost::multi_array<TMatrixT<double>, 2>& spinDensityMatrixElementsCovMatInv() const { return _spinDensityMatrixElementsCovMatInv; }
			const boost::multi_array<double, 6>& spinDensityMatrixElementsCovMatInvArray() const { return _spinDensityMatrixElementsCovMatInvArray; }

			rpwa::resonanceFit::function::useCovarianceMatrix useCovariance() const { return _useCovariance; }

			const boost::multi_array<std::pair<double, double>, 3>& plottingIntensities() const { return _plottingIntensities; }
			const boost::multi_array<std::pair<double, double>, 4>& plottingSpinDensityMatrixElementsReal() const { return _plottingSpinDensityMatrixElementsReal; }
			const boost::multi_array<std::pair<double, double>, 4>& plottingSpinDensityMatrixElementsImag() const { return _plottingSpinDensityMatrixElementsImag; }
			const boost::multi_array<std::pair<double, double>, 4>& plottingPhases() const { return _plottingPhases; }

			const boost::multi_array<std::pair<double, double>, 3>& sysPlottingIntensities() const { return _sysPlottingIntensities; }
			const boost::multi_array<std::pair<double, double>, 4>& sysPlottingSpinDensityMatrixElementsReal() const { return _sysPlottingSpinDensityMatrixElementsReal; }
			const boost::multi_array<std::pair<double, double>, 4>& sysPlottingSpinDensityMatrixElementsImag() const { return _sysPlottingSpinDensityMatrixElementsImag; }
			const boost::multi_array<std::pair<double, double>, 4>& sysPlottingPhases() const { return _sysPlottingPhases; }

		private:

			std::vector<size_t> _nrMassBins;
			boost::multi_array<double, 2> _massBinCenters;

			// data used for fitting

			boost::multi_array<std::complex<double>, 3> _productionAmplitudes;
			boost::multi_array<TMatrixT<double>, 2> _productionAmplitudesCovMatInv;

			boost::multi_array<std::complex<double>, 4> _spinDensityMatrixElements;
			boost::multi_array<TMatrixT<double>, 2> _spinDensityMatrixElementsCovMatInv;
			boost::multi_array<double, 6> _spinDensityMatrixElementsCovMatInvArray;

			rpwa::resonanceFit::function::useCovarianceMatrix _useCovariance;

			// data used for plotting

			boost::multi_array<std::pair<double, double>, 3> _plottingIntensities;
			boost::multi_array<std::pair<double, double>, 4> _plottingSpinDensityMatrixElementsReal;
			boost::multi_array<std::pair<double, double>, 4> _plottingSpinDensityMatrixElementsImag;
			boost::multi_array<std::pair<double, double>, 4> _plottingPhases;

			boost::multi_array<std::pair<double, double>, 3> _sysPlottingIntensities;
			boost::multi_array<std::pair<double, double>, 4> _sysPlottingSpinDensityMatrixElementsReal;
			boost::multi_array<std::pair<double, double>, 4> _sysPlottingSpinDensityMatrixElementsImag;
			boost::multi_array<std::pair<double, double>, 4> _sysPlottingPhases;

		};

	} // end namespace resonanceFit

} // end namespace rpwa

#endif // RESONANCEFIT_DATA_HH
