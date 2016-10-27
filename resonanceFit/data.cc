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
//      implementation of the storage for the data required by the
//      resonance fit
//
//-------------------------------------------------------------------------


#include "data.h"

#include <reportingUtils.hpp>

#include "resonanceFitHelper.h"


rpwa::resonanceFit::data::data(const std::vector<size_t>& nrMassBins,
                               const boost::multi_array<double, 2>& massBinCenters,
                               const boost::multi_array<std::pair<size_t, size_t>, 3>& wavePairMassBinLimits,
                               const boost::multi_array<double, 3>& phaseSpaceIntegrals,
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
                               const boost::multi_array<std::pair<double, double>, 4>& sysPlottingPhases)
	: _nrMassBins(nrMassBins),
	  _massBinCenters(massBinCenters),
	  _wavePairMassBinLimits(wavePairMassBinLimits),
	  _phaseSpaceIntegrals(phaseSpaceIntegrals),
	  _productionAmplitudes(productionAmplitudes),
	  _spinDensityMatrixElements(spinDensityMatrixElements),
	  _useCovariance(useCovariance),
	  _plottingIntensities(plottingIntensities),
	  _plottingSpinDensityMatrixElementsReal(plottingSpinDensityMatrixElementsReal),
	  _plottingSpinDensityMatrixElementsImag(plottingSpinDensityMatrixElementsImag),
	  _plottingPhases(plottingPhases),
	  _sysPlottingIntensities(sysPlottingIntensities),
	  _sysPlottingSpinDensityMatrixElementsReal(sysPlottingSpinDensityMatrixElementsReal),
	  _sysPlottingSpinDensityMatrixElementsImag(sysPlottingSpinDensityMatrixElementsImag),
	  _sysPlottingPhases(sysPlottingPhases)
{
	// get dimensions from one array and make sure that all other arrays
	// have the same dimensions
	const size_t nrBins = _nrMassBins.size();
	if(nrBins == 0) {
		printErr << "number of bins is zero, cannot perform the fit. Aborting..." << std::endl;
		throw;
	}

	const size_t maxMassBins = *(std::max_element(_nrMassBins.begin(), _nrMassBins.end()));
	if(maxMassBins == 0) {
		printErr << "maximal number of mass bins is zero, cannot perform the fit. Aborting..." << std::endl;
		throw;
	}

	const size_t nrWaves = *(_productionAmplitudes.shape()+2);
	if(nrWaves == 0) {
		printErr << "number of waves is zero, cannot perform the fit. Aborting..." << std::endl;
		throw;
	}

	checkSize(_nrMassBins,
	          nrBins, "number of bins is not correct for number of mass bins.");
	checkSize(_massBinCenters,
	          nrBins, "number of bins is not correct for centers of mass bins.",
	          maxMassBins, "maximal number of mass bins is not correct for centers of mass bins.");

	checkSize(_wavePairMassBinLimits,
	          nrBins, "number of bins is not correct for bin ranges of wave pairs.",
	          nrWaves, "maximal number of mass bins is not correct for bin ranges of wave pairs.",
	          nrWaves, "maximal number of mass bins is not correct for bin ranges of wave pairs.");

	checkSize(_phaseSpaceIntegrals,
	          nrBins, "number of bins is not correct for phase-space integrals.",
	          maxMassBins, "maximal number of mass bins is not correct for phase-space integrals.",
	          nrWaves, "number of waves is not correct for phase-space integrals.");

	checkSize(_productionAmplitudes,
	          nrBins, "number of bins is not correct for production amplitudes.",
	          maxMassBins, "maximal number of mass bins is not correct for production amplitudes.",
	          nrWaves, "number of waves is not correct for production amplitudes.");

	// TMatrixT<double>::operator= does not adjust the size of matrices
	_productionAmplitudesCovMatInv.resize(std::vector<size_t>(productionAmplitudesCovMatInv.shape(), productionAmplitudesCovMatInv.shape()+productionAmplitudesCovMatInv.num_dimensions()));
	checkSize(_productionAmplitudesCovMatInv,
	          nrBins, "number of bins is not correct for inverted covariance matrices of production amplitudes.",
	          maxMassBins, "maximal number of mass bins is not correct for inverted covariance matrices of production amplitudes.");
	for(size_t idxBin = 0; idxBin < *(productionAmplitudesCovMatInv.shape()); ++idxBin) {
		for(size_t idxMass = 0; idxMass < *(productionAmplitudesCovMatInv.shape()+1); ++idxMass) {
			if(productionAmplitudesCovMatInv[idxBin][idxMass].GetNrows() == 0 and productionAmplitudesCovMatInv[idxBin][idxMass].GetNcols() == 0) {
				continue;
			}

			_productionAmplitudesCovMatInv[idxBin][idxMass].ResizeTo(productionAmplitudesCovMatInv[idxBin][idxMass]);
			_productionAmplitudesCovMatInv[idxBin][idxMass] = productionAmplitudesCovMatInv[idxBin][idxMass];

			if((Int_t)(2*nrWaves) != _productionAmplitudesCovMatInv[idxBin][idxMass].GetNrows()) {
				printErr << "number of waves is not correct for inverted covariance matrices of production amplitudes. Aborting..." << std::endl;
				throw;
			}
			if((Int_t)(2*nrWaves) != _productionAmplitudesCovMatInv[idxBin][idxMass].GetNcols()) {
				printErr << "number of waves is not correct for inverted covariance matrices of production amplitudes. Aborting..." << std::endl;
				throw;
			}
		}
	}

	checkSize(_spinDensityMatrixElements,
	          nrBins, "number of bins is not correct for spin-density matrix elements.",
	          maxMassBins, "maximal number of mass bins is not correct for spin-density matrix elements.",
	          nrWaves, "number of waves is not correct for spin-density matrix elements.",
	          nrWaves, "number of waves is not correct for spin-density matrix elements.");

	// TMatrixT<double>::operator= does not adjust the size of matrices
	_spinDensityMatrixElementsCovMatInv.resize(std::vector<size_t>(spinDensityMatrixElementsCovMatInv.shape(), spinDensityMatrixElementsCovMatInv.shape()+spinDensityMatrixElementsCovMatInv.num_dimensions()));
	checkSize(_spinDensityMatrixElementsCovMatInv,
	          nrBins, "number of bins is not correct for inverted covariance matrices of spin-density matrix elements.",
	          maxMassBins, "maximal number of mass bins is not correct for inverted covariance matrices of spin-density matrix elements.");
	for(size_t idxBin = 0; idxBin < *(spinDensityMatrixElementsCovMatInv.shape()); ++idxBin) {
		for(size_t idxMass = 0; idxMass < *(spinDensityMatrixElementsCovMatInv.shape()+1); ++idxMass) {
			if(spinDensityMatrixElementsCovMatInv[idxBin][idxMass].GetNrows() == 0 and spinDensityMatrixElementsCovMatInv[idxBin][idxMass].GetNcols() == 0) {
				continue;
			}

			_spinDensityMatrixElementsCovMatInv[idxBin][idxMass].ResizeTo(spinDensityMatrixElementsCovMatInv[idxBin][idxMass]);
			_spinDensityMatrixElementsCovMatInv[idxBin][idxMass] = spinDensityMatrixElementsCovMatInv[idxBin][idxMass];

			if((Int_t)(nrWaves*(nrWaves+1)) != _spinDensityMatrixElementsCovMatInv[idxBin][idxMass].GetNrows()) {
				printErr << "number of waves is not correct for inverted covariance matrices of spin-density matrix elements. Aborting..." << std::endl;
				throw;
			}
			if((Int_t)(nrWaves*(nrWaves+1)) != _spinDensityMatrixElementsCovMatInv[idxBin][idxMass].GetNcols()) {
				printErr << "number of waves is not correct for inverted covariance matrices of spin-density matrix elements. Aborting..." << std::endl;
				throw;
			}
		}
	}

	// copy the diagonal of the covariance matrices of the spin-density
	// matrix elements for faster calculations if the full covariance
	// matrix is not used
	_spinDensityMatrixElementsCovMatInvArray.resize(boost::extents[nrBins][maxMassBins][nrWaves][nrWaves][2][2]);
	for(size_t idxBin = 0; idxBin < nrBins; ++idxBin) {
		for(size_t idxMass = 0; idxMass < maxMassBins; ++idxMass) {
			if(_spinDensityMatrixElementsCovMatInv[idxBin][idxMass].GetNrows() == 0 and _spinDensityMatrixElementsCovMatInv[idxBin][idxMass].GetNcols() == 0) {
				continue;
			}

			for(size_t iWave = 0; iWave < nrWaves; ++iWave) {
				for(size_t jWave = 0; jWave < nrWaves; ++jWave) {
					const size_t base = nrWaves*(nrWaves+1) - (nrWaves-iWave)*(nrWaves-iWave+1) + 2*(jWave-iWave);
					_spinDensityMatrixElementsCovMatInvArray[idxBin][idxMass][iWave][jWave][0][0] = _spinDensityMatrixElementsCovMatInv[idxBin][idxMass](base,   base  );
					_spinDensityMatrixElementsCovMatInvArray[idxBin][idxMass][iWave][jWave][0][1] = _spinDensityMatrixElementsCovMatInv[idxBin][idxMass](base,   base+1);
					_spinDensityMatrixElementsCovMatInvArray[idxBin][idxMass][iWave][jWave][1][0] = _spinDensityMatrixElementsCovMatInv[idxBin][idxMass](base+1, base  );
					_spinDensityMatrixElementsCovMatInvArray[idxBin][idxMass][iWave][jWave][1][1] = _spinDensityMatrixElementsCovMatInv[idxBin][idxMass](base+1, base+1);
				}
			}
		}
	}

	checkSize(_plottingIntensities,
	          nrBins, "number of bins is not correct for intensities for plotting.",
	          maxMassBins, "maximal number of mass bins is not correct for intensities for plotting.",
	          nrWaves, "number of waves is not correct for intensities for plotting.");
	checkSize(_plottingSpinDensityMatrixElementsReal,
	          nrBins, "number of bins is not correct for real part of spin-density matrix elements for plotting.",
	          maxMassBins, "maximal number of mass bins is not correct for real part of spin-density matrix elements for plotting.",
	          nrWaves, "number of waves is not correct for real part of spin-density matrix elements for plotting.",
	          nrWaves, "number of waves is not correct for real part of spin-density matrix elements for plotting.");
	checkSize(_plottingSpinDensityMatrixElementsImag,
	          nrBins, "number of bins is not correct for imaginary part of spin-density matrix elements for plotting.",
	          maxMassBins, "maximal number of mass bins is not correct for imaginary part of spin-density matrix elements for plotting.",
	          nrWaves, "number of waves is not correct for imaginary part of spin-density matrix elements for plotting.",
	          nrWaves, "number of waves is not correct for imaginary part of spin-density matrix elements for plotting.");
	checkSize(_plottingPhases,
	          nrBins, "number of bins is not correct for phases for plotting.",
	          maxMassBins, "maximal number of mass bins is not correct for phases for plotting.",
	          nrWaves, "number of waves is not correct for phases for plotting.",
	          nrWaves, "number of waves is not correct for phases for plotting.");

	checkSize(_sysPlottingIntensities,
	          nrBins, "number of bins is not correct for intensities for plotting of systematic errors.",
	          maxMassBins, "maximal number of mass bins is not correct for intensities for plotting of systematic errors.",
	          nrWaves, "number of waves is not correct for intensities for plotting of systematic errors.");
	checkSize(_sysPlottingSpinDensityMatrixElementsReal,
	          nrBins, "number of bins is not correct for real part of spin-density matrix elements for plotting of systematic errors.",
	          maxMassBins, "maximal number of mass bins is not correct for real part of spin-density matrix elements for plotting of systematic errors.",
	          nrWaves, "number of waves is not correct for real part of spin-density matrix elements for plotting of systematic errors.",
	          nrWaves, "number of waves is not correct for real part of spin-density matrix elements for plotting of systematic errors.");
	checkSize(_sysPlottingSpinDensityMatrixElementsImag,
	          nrBins, "number of bins is not correct for imaginary part of spin-density matrix elements for plotting of systematic errors.",
	          maxMassBins, "maximal number of mass bins is not correct for imaginary part of spin-density matrix elements for plotting of systematic errors.",
	          nrWaves, "number of waves is not correct for imaginary part of spin-density matrix elements for plotting of systematic errors.",
	          nrWaves, "number of waves is not correct for imaginary part of spin-density matrix elements for plotting of systematic errors.");
	checkSize(_sysPlottingPhases,
	          nrBins, "number of bins is not correct for phases for plotting of systematic errors.",
	          maxMassBins, "maximal number of mass bins is not correct for phases for plotting of systematic errors.",
	          nrWaves, "number of waves is not correct for phases for plotting of systematic errors.",
	          nrWaves, "number of waves is not correct for phases for plotting of systematic errors.");
}


bool
rpwa::resonanceFit::data::hasSameMassBinning() const
{
	for(size_t idxBin = 0; idxBin < _nrMassBins.size(); ++idxBin) {
		if(_nrMassBins[idxBin] != _nrMassBins[0]) {
			return false;
		}
		for(size_t idxMass = 0; idxMass < _nrMassBins[idxBin]; ++idxMass) {
			if(_massBinCenters[idxBin][idxMass] != _massBinCenters[0][idxMass]) {
				return false;
			}
		}
	}

	return true;
}
