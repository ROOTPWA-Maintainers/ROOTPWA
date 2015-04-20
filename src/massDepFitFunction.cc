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
//      implementation of function to minimize for the resonance fit
//
//-------------------------------------------------------------------------


#include "massDepFitFunction.h"

#include "massDepFitCache.h"
#include "massDepFitModel.h"
#include "massDepFitParameters.h"
#include "reportingUtils.hpp"


rpwa::massDepFit::function::function(const bool fitProductionAmplitudes,
                                     const rpwa::massDepFit::function::useCovarianceMatrix useCovariance)
	: _fitProductionAmplitudes(fitProductionAmplitudes),
	  _useCovariance(useCovariance)
{
	// if useCovariance has not been overwritten from the command line set
	// reasonable defaults depending on what to fit to
	if (_useCovariance == useCovarianceMatrixDefault) {
		if (fitProductionAmplitudes) {
			_useCovariance = useFullCovarianceMatrix;
		} else {
			_useCovariance = useComplexDiagnalElementsOnly;
		}
	}
	assert(_useCovariance != useCovarianceMatrixDefault);

	std::ostringstream output;
	output << "created 'function' object for a fit to the ";
	if (fitProductionAmplitudes) {
		output << "production amplitudes";
	} else {
		output << "spin-density matrix";
	}
	output << " using ";
	if (_useCovariance == useDiagnalElementsOnly) {
		output << "only the diagonal elements of the covariance matrix";
	} else if (_useCovariance == useComplexDiagnalElementsOnly) {
		output << "blocks of 2x2 along the diagonal of the covariance matrix corresponding to a complex number";
	} else if (_useCovariance == useFullCovarianceMatrix) {
		output << "the full covariance matrix";
	} else {
		assert(false);
	}
	output << ".";
	printInfo << output.str() << std::endl;
}


bool
rpwa::massDepFit::function::init(rpwa::massDepFit::model* compset,
                                 const std::vector<double>& massBinCenters,
                                 const boost::multi_array<std::complex<double>, 3>& productionAmplitudes,
                                 const boost::multi_array<double, 6>& productionAmplitudesCovariance,
                                 const boost::multi_array<std::complex<double>, 4>& spinDensityMatrices,
                                 const boost::multi_array<double, 6>& spinDensityCovarianceMatrices,
                                 const boost::multi_array<std::pair<size_t, size_t>, 2>& wavePairMassBinLimits)
{
	if(not _fitProductionAmplitudes && _useCovariance == useFullCovarianceMatrix) {
		printErr << "cannot use full covariance matrix while fitting to spin-density matrix." << std::endl;
		return false;
	}

	_compset = compset;

	_massBinCenters = massBinCenters;

	_productionAmplitudes.resize(std::vector<size_t>(productionAmplitudes.shape(), productionAmplitudes.shape()+productionAmplitudes.num_dimensions()));
	_productionAmplitudes = productionAmplitudes;
	_productionAmplitudesCovariance.resize(std::vector<size_t>(productionAmplitudesCovariance.shape(), productionAmplitudesCovariance.shape()+productionAmplitudesCovariance.num_dimensions()));
	_productionAmplitudesCovariance = productionAmplitudesCovariance;

	_spinDensityMatrices.resize(std::vector<size_t>(spinDensityMatrices.shape(), spinDensityMatrices.shape()+spinDensityMatrices.num_dimensions()));
	_spinDensityMatrices = spinDensityMatrices;
	_spinDensityCovarianceMatrices.resize(std::vector<size_t>(spinDensityCovarianceMatrices.shape(), spinDensityCovarianceMatrices.shape()+spinDensityCovarianceMatrices.num_dimensions()));
	_spinDensityCovarianceMatrices = spinDensityCovarianceMatrices;

	_wavePairMassBinLimits.resize(std::vector<size_t>(wavePairMassBinLimits.shape(), wavePairMassBinLimits.shape()+wavePairMassBinLimits.num_dimensions()));
	_wavePairMassBinLimits = wavePairMassBinLimits;

	_idxAnchorWave = _compset->getAnchorWave();

	_nrBins = _spinDensityMatrices.size();
	_nrMassBins = _massBinCenters.size();
	_nrWaves = _wavePairMassBinLimits.size();

	_idxMassMin = _nrMassBins;
	_idxMassMax = 0;
	for(size_t idxBin=0; idxBin<_nrBins; ++idxBin) {
		for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
			_idxMassMin = std::min(_idxMassMin, _wavePairMassBinLimits[idxWave][idxWave].first);
			_idxMassMax = std::max(_idxMassMax, _wavePairMassBinLimits[idxWave][idxWave].second);
		}
	}

	if(_fitProductionAmplitudes) {
		// do some stuff specific to the fit to the production amplitudes

		// get a list of waves that are zero (those have to be excluded
		// form the inversion of the covariance matrix below) and test
		// that the anchor wave is non-zero over the complete fit range
		bool zeroAnchorWave = false;
		boost::multi_array<std::vector<size_t>, 2> zeroWaves(boost::extents[_nrBins][_nrMassBins]);
		for(size_t idxBin=0; idxBin<_nrBins; ++idxBin) {
			for(size_t idxMass=_idxMassMin; idxMass<=_idxMassMax; ++idxMass) {
				for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
					bool zeroThisWave = true;
					zeroThisWave &= (_productionAmplitudes[idxBin][idxMass][idxWave].real() == 0.);
					zeroThisWave &= (_productionAmplitudes[idxBin][idxMass][idxWave].imag() == 0.);
					zeroThisWave &= (_productionAmplitudesCovariance[idxBin][idxMass][idxWave][idxWave][0][0] == 0.);
					zeroThisWave &= (_productionAmplitudesCovariance[idxBin][idxMass][idxWave][idxWave][0][1] == 0.);
					zeroThisWave &= (_productionAmplitudesCovariance[idxBin][idxMass][idxWave][idxWave][1][0] == 0.);
					zeroThisWave &= (_productionAmplitudesCovariance[idxBin][idxMass][idxWave][idxWave][1][1] == 0.);

					if(zeroThisWave || idxMass < _wavePairMassBinLimits[idxWave][idxWave].first || idxMass > _wavePairMassBinLimits[idxWave][idxWave].second) {
						zeroWaves[idxBin][idxMass].push_back(idxWave);
					}

					if(idxWave == _idxAnchorWave) {
						zeroAnchorWave |= zeroThisWave;
					}
				}
			}
		}

		// error if anchor wave is zero in one mass bin
		if(zeroAnchorWave) {
			printErr << "production amplitudes of anchor wave zero in same mass bins (mass limit in mass-independent fit)." << std::endl;
			return false;
		}

		// test if the anchor wave is real valued
		// test that any non-anchor wave is not real valued
		bool realAnchorWave = true;
		bool realOtherWaves = false;
		for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
			bool realThisWave = true;
			for(size_t idxBin=0; idxBin<_nrBins; ++idxBin) {
				for(size_t idxMass=0; idxMass<_nrMassBins; ++idxMass) {
					realThisWave &= (_productionAmplitudes[idxBin][idxMass][idxWave].imag() == 0.);
					realThisWave &= (_productionAmplitudesCovariance[idxBin][idxMass][idxWave][idxWave][0][1] == 0.);
					realThisWave &= (_productionAmplitudesCovariance[idxBin][idxMass][idxWave][idxWave][1][0] == 0.);
					realThisWave &= (_productionAmplitudesCovariance[idxBin][idxMass][idxWave][idxWave][1][1] == 0.);
				}
			}

			if(idxWave == _idxAnchorWave) {
				realAnchorWave &= realThisWave;
			} else {
				realOtherWaves |= realThisWave;
			}
		}

		// error if any non-anchor wave is real
		if(realOtherWaves) {
			printErr << "production amplitudes cannot be fitted if a non-anchor wave is real valued." << std::endl;
			return false;
		}

		_productionAmplitudesCovMatInv.resize(boost::extents[_nrBins][_nrMassBins]);
		for(size_t idxBin=0; idxBin<_nrBins; ++idxBin) {
			for(size_t idxMass=_idxMassMin; idxMass<=_idxMassMax; ++idxMass) {
				// determine whether the anchor wave should be used in the current bin
				bool skipAnchor = false;
				for (std::vector<size_t>::const_iterator it=zeroWaves[idxBin][idxMass].begin(); it!=zeroWaves[idxBin][idxMass].end(); ++it) {
					if (*it == _idxAnchorWave) {
						skipAnchor = true;
					}
				}

				// import covariance matrix of production amplitudes
				const size_t matrixSize = 2 * (_nrWaves - zeroWaves[idxBin][idxMass].size()) - (skipAnchor ? 0 : 1);
				TMatrixT<double> reducedCovMat(matrixSize, matrixSize);

				if(realAnchorWave) {
					for(size_t idxWave=0, idxSkip=0; idxWave<_nrWaves; ++idxWave) {
						if(idxSkip < zeroWaves[idxBin][idxMass].size() && zeroWaves[idxBin][idxMass][idxSkip] == idxWave) {
							++idxSkip;
							continue;
						}

						for(size_t jdxWave=0, jdxSkip=0; jdxWave<_nrWaves; ++jdxWave) {
							if(jdxSkip < zeroWaves[idxBin][idxMass].size() && zeroWaves[idxBin][idxMass][jdxSkip] == jdxWave) {
								++jdxSkip;
								continue;
							}

							const Int_t rowSkip = 2*(idxWave-idxSkip) + ((idxWave>_idxAnchorWave && !skipAnchor) ? -1 : 0);
							const Int_t colSkip = 2*(jdxWave-jdxSkip) + ((jdxWave>_idxAnchorWave && !skipAnchor) ? -1 : 0);

							reducedCovMat(rowSkip + 0, colSkip + 0) = _productionAmplitudesCovariance[idxBin][idxMass][idxWave][jdxWave][0][0];
							if(jdxWave != _idxAnchorWave) {
								reducedCovMat(rowSkip + 0, colSkip + 1) = _productionAmplitudesCovariance[idxBin][idxMass][idxWave][jdxWave][0][1];
							}
							if(idxWave != _idxAnchorWave) {
								reducedCovMat(rowSkip + 1, colSkip + 0) = _productionAmplitudesCovariance[idxBin][idxMass][idxWave][jdxWave][1][0];
							}
							if(idxWave != _idxAnchorWave && jdxWave != _idxAnchorWave) {
								reducedCovMat(rowSkip + 1, colSkip + 1) = _productionAmplitudesCovariance[idxBin][idxMass][idxWave][jdxWave][1][1];
							}
						}
					}
				} else {
					TMatrixT<double> covariance(matrixSize + (skipAnchor ? 0 : 1), matrixSize + (skipAnchor ? 0 : 1));
					TMatrixT<double> jacobian(matrixSize, matrixSize + (skipAnchor ? 0 : 1));

					for(size_t idxWave=0, idxSkip=0; idxWave<_nrWaves; ++idxWave) {
						if(idxSkip < zeroWaves[idxBin][idxMass].size() && zeroWaves[idxBin][idxMass][idxSkip] == idxWave) {
							++idxSkip;
							continue;
						}

						for(size_t jdxWave=0, jdxSkip=0; jdxWave<_nrWaves; ++jdxWave) {
							if(jdxSkip < zeroWaves[idxBin][idxMass].size() && zeroWaves[idxBin][idxMass][jdxSkip] == jdxWave) {
								++jdxSkip;
								continue;
							}

							covariance(2*(idxWave-idxSkip) + 0, 2*(jdxWave-jdxSkip) + 0) = _productionAmplitudesCovariance[idxBin][idxMass][idxWave][jdxWave][0][0];
							covariance(2*(idxWave-idxSkip) + 0, 2*(jdxWave-jdxSkip) + 1) = _productionAmplitudesCovariance[idxBin][idxMass][idxWave][jdxWave][0][1];
							covariance(2*(idxWave-idxSkip) + 1, 2*(jdxWave-jdxSkip) + 0) = _productionAmplitudesCovariance[idxBin][idxMass][idxWave][jdxWave][1][0];
							covariance(2*(idxWave-idxSkip) + 1, 2*(jdxWave-jdxSkip) + 1) = _productionAmplitudesCovariance[idxBin][idxMass][idxWave][jdxWave][1][1];

							const double n = abs(_productionAmplitudes[idxBin][idxMass][_idxAnchorWave]);
							const double n3 = std::pow(n, 3);
							const double xa1 = _productionAmplitudes[idxBin][idxMass][_idxAnchorWave].real();
							const double xa2 = _productionAmplitudes[idxBin][idxMass][_idxAnchorWave].imag();
							const double xi1 = _productionAmplitudes[idxBin][idxMass][idxWave].real();
							const double xi2 = _productionAmplitudes[idxBin][idxMass][idxWave].imag();

							const Int_t rowSkip = 2*(idxWave-idxSkip) + ((idxWave>_idxAnchorWave && !skipAnchor) ? -1 : 0);
							if(idxWave == _idxAnchorWave && jdxWave == _idxAnchorWave) {
								jacobian(rowSkip + 0, 2*(jdxWave-jdxSkip) + 0) = xa1 / n;
								jacobian(rowSkip + 0, 2*(jdxWave-jdxSkip) + 1) = xa2 / n;
							} else if(jdxWave == _idxAnchorWave) {
								jacobian(rowSkip + 0, 2*(jdxWave-jdxSkip) + 0) = xi1 / n - xa1 * (xi1*xa1 + xi2*xa2) / n3;
								jacobian(rowSkip + 0, 2*(jdxWave-jdxSkip) + 1) = xi2 / n - xa2 * (xi1*xa1 + xi2*xa2) / n3;
								if(idxWave != _idxAnchorWave) {
									jacobian(rowSkip + 1, 2*(jdxWave-jdxSkip) + 0) =   xi2 / n - xa1 * (xi2*xa1 - xi1*xa2) / n3;
									jacobian(rowSkip + 1, 2*(jdxWave-jdxSkip) + 1) = - xi1 / n - xa2 * (xi2*xa1 - xi1*xa2) / n3;
								}
							} else if(idxWave == jdxWave) {
								jacobian(rowSkip + 0, 2*(jdxWave-jdxSkip) + 0) =   xa1 / n;
								jacobian(rowSkip + 0, 2*(jdxWave-jdxSkip) + 1) =   xa2 / n;
								jacobian(rowSkip + 1, 2*(jdxWave-jdxSkip) + 0) = - xa2 / n;
								jacobian(rowSkip + 1, 2*(jdxWave-jdxSkip) + 1) =   xa1 / n;
							}
						}
					}

					TMatrixT<double> jacobianT(TMatrixT<double>::kTransposed, jacobian);

					reducedCovMat = jacobian * covariance * jacobianT;
				}

				reducedCovMat.Invert();

				// import covariance matrix of production amplitudes
				_productionAmplitudesCovMatInv[idxBin][idxMass].ResizeTo(2*_nrWaves - 1, 2*_nrWaves - 1);

				for(size_t idxWave=0, idxSkip=0; idxWave<_nrWaves; ++idxWave) {
					if(idxSkip < zeroWaves[idxBin][idxMass].size() && zeroWaves[idxBin][idxMass][idxSkip] == idxWave) {
						++idxSkip;
						continue;
					}

					for(size_t jdxWave=0, jdxSkip=0; jdxWave<_nrWaves; ++jdxWave) {
						if(jdxSkip < zeroWaves[idxBin][idxMass].size() && zeroWaves[idxBin][idxMass][jdxSkip] == jdxWave) {
							++jdxSkip;
							continue;
						}

						if(idxWave != jdxWave && _useCovariance != useFullCovarianceMatrix) {
							continue;
						}

						const Int_t row = 2*idxWave + (idxWave>_idxAnchorWave ? -1 : 0);
						const Int_t col = 2*jdxWave + (jdxWave>_idxAnchorWave ? -1 : 0);
						const Int_t rowSkip = 2*(idxWave-idxSkip) + ((idxWave>_idxAnchorWave && !skipAnchor) ? -1 : 0);
						const Int_t colSkip = 2*(jdxWave-jdxSkip) + ((jdxWave>_idxAnchorWave && !skipAnchor) ? -1 : 0);

						_productionAmplitudesCovMatInv[idxBin][idxMass](row + 0, col + 0) = reducedCovMat(rowSkip + 0, colSkip + 0);
						if(jdxWave != _idxAnchorWave && _useCovariance != useDiagnalElementsOnly) {
							_productionAmplitudesCovMatInv[idxBin][idxMass](row + 0, col + 1) = reducedCovMat(rowSkip + 0, colSkip + 1);
						}
						if(idxWave != _idxAnchorWave && _useCovariance != useDiagnalElementsOnly) {
							_productionAmplitudesCovMatInv[idxBin][idxMass](row + 1, col + 0) = reducedCovMat(rowSkip + 1, colSkip + 0);
						}
						if(idxWave != _idxAnchorWave && jdxWave != _idxAnchorWave) {
							_productionAmplitudesCovMatInv[idxBin][idxMass](row + 1, col + 1) = reducedCovMat(rowSkip + 1, colSkip + 1);
						}
					}
				}
			}
		}

		// modify measured production amplitude such that the anchor wave is always real and positive
		for(size_t idxBin=0; idxBin<_nrBins; ++idxBin) {
			for(size_t idxMass=_idxMassMin; idxMass<=_idxMassMax; ++idxMass) {
				const std::complex<double> anchorPhase = _productionAmplitudes[idxBin][idxMass][_idxAnchorWave] / abs(_productionAmplitudes[idxBin][idxMass][_idxAnchorWave]);
				for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
					_productionAmplitudes[idxBin][idxMass][idxWave] /= anchorPhase;
				}
			}
		}
	} else {
		// do some stuff specific to the fit to the spin-density matrix

		// get a list of waves that are zero (those have to be excluded
		// form the inversion of the covariance matrix below)
		boost::multi_array<std::vector<size_t>, 2> zeroWaves(boost::extents[_nrBins][_nrMassBins]);
		for(size_t idxBin=0; idxBin<_nrBins; ++idxBin) {
			for(size_t idxMass=_idxMassMin; idxMass<=_idxMassMax; ++idxMass) {
				for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
					bool zeroThisWave = true;
					for(size_t jdxWave=0; jdxWave<_nrWaves; ++jdxWave) {
						zeroThisWave &= (_spinDensityMatrices[idxBin][idxMass][idxWave][jdxWave].real() == 0.);
						zeroThisWave &= (_spinDensityMatrices[idxBin][idxMass][idxWave][jdxWave].imag() == 0.);
						zeroThisWave &= (_spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][0][0] == 0.);
						zeroThisWave &= (_spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][0][1] == 0.);
						zeroThisWave &= (_spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][1][0] == 0.);
						zeroThisWave &= (_spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][1][1] == 0.);
					}

					if(zeroThisWave || idxMass < _wavePairMassBinLimits[idxWave][idxWave].first || idxMass > _wavePairMassBinLimits[idxWave][idxWave].second) {
						zeroWaves[idxBin][idxMass].push_back(idxWave);
					}
				}
			}
		}

		_spinDensityMatricesCovMatInv.resize(boost::extents[_nrBins][_nrMassBins]);
		for(size_t idxBin=0; idxBin<_nrBins; ++idxBin) {
			for(size_t idxMass=_idxMassMin; idxMass<=_idxMassMax; ++idxMass) {
				// import covariance matrix of spin-density matrix elements
				const size_t matrixSize((_nrWaves - zeroWaves[idxBin][idxMass].size())*(_nrWaves - zeroWaves[idxBin][idxMass].size()));
				TMatrixT<double> reducedCovMat(matrixSize, matrixSize);

				for(size_t idxWave=0, idxSkip=0; idxWave<_nrWaves; ++idxWave) {
					if(idxSkip < zeroWaves[idxBin][idxMass].size() && zeroWaves[idxBin][idxMass][idxSkip] == idxWave) {
						++idxSkip;
						continue;
					}

					for(size_t jdxWave=idxWave, jdxSkip=idxSkip; jdxWave<_nrWaves; ++jdxWave) {
						if(jdxSkip < zeroWaves[idxBin][idxMass].size() && zeroWaves[idxBin][idxMass][jdxSkip] == jdxWave) {
							++jdxSkip;
							continue;
						}

						const Int_t idxIJskip = (idxWave-idxSkip)*(_nrWaves-zeroWaves[idxBin][idxMass].size()) + (jdxWave-jdxSkip);
						const Int_t idxJIskip = (jdxWave-jdxSkip)*(_nrWaves-zeroWaves[idxBin][idxMass].size()) + (idxWave-idxSkip);

						if(idxWave==jdxWave) {
							reducedCovMat(idxIJskip, idxIJskip) = _spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][0][0];
						} else {
							if (_useCovariance == useDiagnalElementsOnly) {
								reducedCovMat(idxIJskip, idxIJskip) = _spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][0][0];
								reducedCovMat(idxJIskip, idxJIskip) = _spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][1][1];
							} else if(_useCovariance == useComplexDiagnalElementsOnly) {
								reducedCovMat(idxIJskip, idxIJskip) = _spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][0][0];
								reducedCovMat(idxIJskip, idxJIskip) = _spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][0][1];
								reducedCovMat(idxJIskip, idxIJskip) = _spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][1][0];
								reducedCovMat(idxJIskip, idxJIskip) = _spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][1][1];
							} else {
								// this should have returned an error before
								assert(false);
							}
						}
					}
				}

				reducedCovMat.Invert();

				// import covariance matrix of spin-density matrix elements
				_spinDensityMatricesCovMatInv[idxBin][idxMass].ResizeTo(_nrWaves*_nrWaves, _nrWaves*_nrWaves);

				for(size_t idxWave=0, idxSkip=0; idxWave<_nrWaves; ++idxWave) {
					if(idxSkip < zeroWaves[idxBin][idxMass].size() && zeroWaves[idxBin][idxMass][idxSkip] == idxWave) {
						++idxSkip;
						continue;
					}

					for(size_t jdxWave=idxWave, jdxSkip=idxSkip; jdxWave<_nrWaves; ++jdxWave) {
						if(jdxSkip < zeroWaves[idxBin][idxMass].size() && zeroWaves[idxBin][idxMass][jdxSkip] == jdxWave) {
							++jdxSkip;
							continue;
						}

						const Int_t idxIJ = idxWave*_nrWaves + jdxWave;
						const Int_t idxJI = jdxWave*_nrWaves + idxWave;
						const Int_t idxIJskip = (idxWave-idxSkip)*(_nrWaves-zeroWaves[idxBin][idxMass].size()) + (jdxWave-jdxSkip);
						const Int_t idxJIskip = (jdxWave-jdxSkip)*(_nrWaves-zeroWaves[idxBin][idxMass].size()) + (idxWave-idxSkip);

						if(idxWave==jdxWave) {
							_spinDensityMatricesCovMatInv[idxBin][idxMass](idxIJ, idxIJ) = reducedCovMat(idxIJskip, idxIJskip);
						} else {
							_spinDensityMatricesCovMatInv[idxBin][idxMass](idxIJ, idxIJ) = reducedCovMat(idxIJskip, idxIJskip);
							_spinDensityMatricesCovMatInv[idxBin][idxMass](idxIJ, idxJI) = reducedCovMat(idxIJskip, idxJIskip);
							_spinDensityMatricesCovMatInv[idxBin][idxMass](idxJI, idxIJ) = reducedCovMat(idxJIskip, idxIJskip);
							_spinDensityMatricesCovMatInv[idxBin][idxMass](idxJI, idxJI) = reducedCovMat(idxJIskip, idxJIskip);
						}
					}
				}
			}
		}
	}

	return true;
}


size_t
rpwa::massDepFit::function::getNrParameters() const
{
	return _compset->getNrParameters();
}


size_t
rpwa::massDepFit::function::getNrDataPoints() const
{
	size_t nrPts(0);

	if(_fitProductionAmplitudes) {
		// calculate data points:
		// * production amplitudes in general are complex numbers
		// * for the anchor wave it might be real
		// * remember (Re,Im) => factor 2
		for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
			nrPts += _wavePairMassBinLimits[idxWave][idxWave].second - _wavePairMassBinLimits[idxWave][idxWave].first + 1;
			if(idxWave != _idxAnchorWave) {
				nrPts += _wavePairMassBinLimits[idxWave][idxWave].second - _wavePairMassBinLimits[idxWave][idxWave].first + 1;
			}
		}
	} else {
		// calculate data points:
		// * diagonal elements are real numbers
		// * off-diagonal elements are complex numbers
		// * remember (Re,Im) => factor 2
		// * diagonal elements are only checked once, off-diagonal elements with
		//   the two different combinations (i,j) and (j,i)
		for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
			for(size_t jdxWave=0; jdxWave<_nrWaves; ++jdxWave) {
				nrPts += _wavePairMassBinLimits[idxWave][jdxWave].second - _wavePairMassBinLimits[idxWave][jdxWave].first + 1;
			}
		}
	}

	nrPts *= _nrBins;

	return nrPts;
}


double
rpwa::massDepFit::function::chiSquare(const std::vector<double>& par) const
{
#if __cplusplus >= 201103L
	return chiSquare(par.data());
#else
	return chiSquare(&par[0]);
#endif
}


double
rpwa::massDepFit::function::chiSquare(const double* par) const
{
	rpwa::massDepFit::parameters fitParameters(_compset->getNrComponents()+1,           // nr components + final-state mass-dependence
	                                           _compset->getMaxChannelsInComponent(),
	                                           _compset->getMaxParametersInComponent(),
	                                           _nrBins);
	rpwa::massDepFit::cache cache(_nrWaves,
	                              _compset->getNrComponents()+1,           // nr components + final-state mass-dependence
	                              _compset->getMaxChannelsInComponent(),
	                              _nrBins,
	                              _nrMassBins);

	// import parameters (couplings, branchings, resonance parameters, ...)
	_compset->importParameters(par, fitParameters, cache);

	return chiSquare(fitParameters, cache);
}


double
rpwa::massDepFit::function::chiSquare(const rpwa::massDepFit::parameters& fitParameters,
                                      rpwa::massDepFit::cache& cache) const
{
	if(_fitProductionAmplitudes) {
		return chiSquareProductionAmplitudes(fitParameters, cache);
	} else {
		return chiSquareSpinDensityMatrix(fitParameters, cache);
	}
}


double
rpwa::massDepFit::function::logLikelihood(const std::vector<double>& par) const
{
	return -0.5 * chiSquare(par);
}


double
rpwa::massDepFit::function::logLikelihood(const double* par) const
{
	return -0.5 * chiSquare(par);
}


double
rpwa::massDepFit::function::logLikelihood(const rpwa::massDepFit::parameters& fitParameters,
                                          rpwa::massDepFit::cache& cache) const
{
	return -0.5 * chiSquare(fitParameters, cache);
}


double
rpwa::massDepFit::function::chiSquareProductionAmplitudes(const rpwa::massDepFit::parameters& fitParameters,
                                                          rpwa::massDepFit::cache& cache) const
{
	double chi2=0;

	// loop over bins
	for(unsigned idxBin=0; idxBin<_nrBins; ++idxBin) {
		// loop over mass-bins
		for(unsigned idxMass=_idxMassMin; idxMass<=_idxMassMax; ++idxMass) {
			const double mass = _massBinCenters[idxMass];

			// phase of fit in anchor wave
			const std::complex<double> anchorFit = _compset->productionAmplitude(fitParameters, cache, _idxAnchorWave, idxBin, mass, idxMass);
			const std::complex<double> anchorFitPhase = anchorFit / abs(anchorFit);

			TVectorT<double> prodAmpDiffVect(2*_nrWaves - 1);

			// sum over the contributions to chi2
			for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
				// check that this mass bin should be taken into account for this
				// combination of waves
				if(idxMass < _wavePairMassBinLimits[idxWave][idxWave].first || idxMass > _wavePairMassBinLimits[idxWave][idxWave].second) {
					continue;
				}

				// calculate target spin density matrix element
				const std::complex<double> prodAmpFit = _compset->productionAmplitude(fitParameters, cache, idxWave, idxBin, mass, idxMass) / anchorFitPhase;

				const std::complex<double> prodAmpDiff = prodAmpFit - _productionAmplitudes[idxBin][idxMass][idxWave];

				const Int_t row = 2*idxWave + (idxWave>_idxAnchorWave ? -1 : 0);
				prodAmpDiffVect(row + 0) = prodAmpDiff.real();
				if(idxWave != _idxAnchorWave) {
					prodAmpDiffVect(row + 1) = prodAmpDiff.imag();
				}
			} // end loop over idxWave

			chi2 += _productionAmplitudesCovMatInv[idxBin][idxMass].Similarity(prodAmpDiffVect);
		} // end loop over mass-bins
	} // end loop over bins

	return chi2;
}


double
rpwa::massDepFit::function::chiSquareSpinDensityMatrix(const rpwa::massDepFit::parameters& fitParameters,
                                                       rpwa::massDepFit::cache& cache) const
{
	double chi2=0;

	// loop over bins
	for(unsigned idxBin=0; idxBin<_nrBins; ++idxBin) {
		// loop over mass-bins
		for(unsigned idxMass=_idxMassMin; idxMass<=_idxMassMax; ++idxMass) {
			const double mass = _massBinCenters[idxMass];

			// sum over the contributions to chi2 -> rho_ij
			for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
				for(size_t jdxWave=idxWave; jdxWave<_nrWaves; ++jdxWave) {
					// check that this mass bin should be taken into account for this
					// combination of waves
					if(idxMass < _wavePairMassBinLimits[idxWave][jdxWave].first || idxMass > _wavePairMassBinLimits[idxWave][jdxWave].second) {
						continue;
					}

					// calculate target spin density matrix element
					const std::complex<double> rhoFit = _compset->spinDensityMatrix(fitParameters, cache, idxWave, jdxWave, idxBin, mass, idxMass);

					const std::complex<double> rhoDiff = rhoFit - _spinDensityMatrices[idxBin][idxMass][idxWave][jdxWave];

					double dchi;
					if(idxWave==jdxWave) {
						dchi = norm(rhoDiff) / _spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][0][0];
					} else {
						if (_useCovariance == useDiagnalElementsOnly) {
							dchi  = rhoDiff.real()*rhoDiff.real() / _spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][0][0];
							dchi += rhoDiff.imag()*rhoDiff.imag() / _spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][1][1];
						} else if(_useCovariance == useComplexDiagnalElementsOnly) {
							dchi  = rhoDiff.real()*rhoDiff.real() * _spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][1][1];
							dchi -= rhoDiff.real()*rhoDiff.imag() * _spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][0][1];
							dchi -= rhoDiff.real()*rhoDiff.imag() * _spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][1][0];
							dchi += rhoDiff.imag()*rhoDiff.imag() * _spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][0][0];

							dchi /= _spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][0][0]*_spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][1][1]
							        - _spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][0][1]*_spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][1][0];
						} else {
							// this should have returned an error during the call to init()
							assert(false);
						}
					}
					chi2 += dchi;
				} // end loop over jdxWave
			} // end loop over idxWave
		} // end loop over mass-bins
	} // end loop over bins

	return chi2;
}
