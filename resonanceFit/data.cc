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


namespace {


	template<typename T>
	void
	checkSize(const std::vector<T>& obj,
	          const size_t dim1,
	          const std::string& msg1)
	{
		if(dim1 != obj.size()) {
			printErr << msg1 << " expected: " << dim1 << ", found: " << obj.size() << ". Aborting..." << std::endl;
			throw;
		}
	}


	template<typename T, const size_t dim = 0>
	void
	checkSize(const boost::multi_array<T, 1+dim>& obj,
	          const size_t dim1,
	          const std::string& msg1)
	{
		if(dim1 != *(obj.shape()+dim)) {
			printErr << msg1 << " expected: " << dim1 << ", found: " << *(obj.shape()+dim) << ". Aborting..." << std::endl;
			throw;
		}
	}


	template<typename T, const size_t dim = 0>
	void
	checkSize(const boost::multi_array<T, 2+dim>& obj,
	          const size_t dim1,
	          const std::string& msg1,
	          const size_t dim2,
	          const std::string& msg2)
	{
		if(dim1 != *(obj.shape()+dim)) {
			printErr << msg1 << " expected: " << dim1 << ", found: " << *(obj.shape()+dim) << ". Aborting..." << std::endl;
			throw;
		}
		checkSize<T, dim+1>(obj, dim2, msg2);
	}


	template<typename T, const size_t dim = 0>
	void
	checkSize(const boost::multi_array<T, 3+dim>& obj,
	          const size_t dim1,
	          const std::string& msg1,
	          const size_t dim2,
	          const std::string& msg2,
	          const size_t dim3,
	          const std::string& msg3)
	{
		if(dim1 != *(obj.shape()+dim)) {
			printErr << msg1 << " expected: " << dim1 << ", found: " << *(obj.shape()+dim) << ". Aborting..." << std::endl;
			throw;
		}
		checkSize<T, dim+1>(obj, dim2, msg2, dim3, msg3);
	}


	template<typename T, const size_t dim = 0>
	void
	checkSize(const boost::multi_array<T, 4+dim>& obj,
	          const size_t dim1,
	          const std::string& msg1,
	          const size_t dim2,
	          const std::string& msg2,
	          const size_t dim3,
	          const std::string& msg3,
	          const size_t dim4,
	          const std::string& msg4)
	{
		if(dim1 != *(obj.shape()+dim)) {
			printErr << msg1 << " expected: " << dim1 << ", found: " << *(obj.shape()+dim) << ". Aborting..." << std::endl;
			throw;
		}
		checkSize<T, dim+1>(obj, dim2, msg2, dim3, msg3, dim4, msg4);
	}


}


rpwa::resonanceFit::data::data(const std::vector<size_t>& nrMassBins,
                               const boost::multi_array<double, 2>& massBinCenters,
                               const boost::multi_array<std::complex<double>, 3>& productionAmplitudes,
                               const boost::multi_array<TMatrixT<double>, 2>& productionAmplitudesCovMatInv,
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
	  _productionAmplitudes(productionAmplitudes),
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
