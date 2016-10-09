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
//      implementation of function to minimize for the resonance fit
//
//-------------------------------------------------------------------------


#include "function.h"

#include "cache.h"
#include "components.h"
#include "model.h"
#include "parameters.h"
#include "reportingUtils.hpp"


rpwa::resonanceFit::function::function(const bool fitProductionAmplitudes,
                                       const rpwa::resonanceFit::function::useCovarianceMatrix useCovariance)
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
rpwa::resonanceFit::function::init(const rpwa::resonanceFit::modelConstPtr& fitModel,
                                   const std::vector<size_t>& nrMassBins,
                                   const boost::multi_array<double, 2>& massBinCenters,
                                   const boost::multi_array<std::complex<double>, 3>& productionAmplitudes,
                                   const boost::multi_array<TMatrixT<double>, 2>& productionAmplitudesCovariance,
                                   const boost::multi_array<std::complex<double>, 4>& spinDensityMatrices,
                                   const boost::multi_array<TMatrixT<double>, 2>& spinDensityCovarianceMatrices,
                                   const boost::multi_array<std::pair<size_t, size_t>, 3>& wavePairMassBinLimits)
{
	if(not _fitProductionAmplitudes && _useCovariance == useFullCovarianceMatrix) {
		printErr << "cannot use full covariance matrix while fitting to spin-density matrix." << std::endl;
		return false;
	}

	_fitModel = fitModel;

	_nrMassBins = nrMassBins;

	_massBinCenters.resize(std::vector<size_t>(massBinCenters.shape(), massBinCenters.shape()+massBinCenters.num_dimensions()));
	_massBinCenters = massBinCenters;

	_productionAmplitudes.resize(std::vector<size_t>(productionAmplitudes.shape(), productionAmplitudes.shape()+productionAmplitudes.num_dimensions()));
	_productionAmplitudes = productionAmplitudes;

	_spinDensityMatrices.resize(std::vector<size_t>(spinDensityMatrices.shape(), spinDensityMatrices.shape()+spinDensityMatrices.num_dimensions()));
	_spinDensityMatrices = spinDensityMatrices;

	_wavePairMassBinLimits.resize(std::vector<size_t>(wavePairMassBinLimits.shape(), wavePairMassBinLimits.shape()+wavePairMassBinLimits.num_dimensions()));
	_wavePairMassBinLimits = wavePairMassBinLimits;

	_idxAnchorWave = _fitModel->getAnchorWave();

	_nrBins = _spinDensityMatrices.shape()[0];
	_maxMassBins = _spinDensityMatrices.shape()[1];
	_nrWaves = _spinDensityMatrices.shape()[2];

	_idxMassMax.resize(_nrBins);
	_idxMassMin.resize(_nrBins);
	for(size_t idxBin=0; idxBin<_nrBins; ++idxBin) {
		_idxMassMin[idxBin] = _nrMassBins[idxBin];
		_idxMassMax[idxBin] = 0;
		for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
			_idxMassMin[idxBin] = std::min(_idxMassMin[idxBin], _wavePairMassBinLimits[idxBin][idxWave][idxWave].first);
			_idxMassMax[idxBin] = std::max(_idxMassMax[idxBin], _wavePairMassBinLimits[idxBin][idxWave][idxWave].second);
		}
	}

	if(_fitProductionAmplitudes) {
		// do some stuff specific to the fit to the production amplitudes

		// get a list of waves that are zero (those have to be excluded
		// from the inversion of the covariance matrix below) and test
		// that the anchor wave is non-zero over the complete fit range
		bool zeroAnchorWave = false;
		boost::multi_array<std::vector<size_t>, 2> zeroWaves(boost::extents[_nrBins][_maxMassBins]);
		for(size_t idxBin=0; idxBin<_nrBins; ++idxBin) {
			for(size_t idxMass = _idxMassMin[idxBin]; idxMass <= _idxMassMax[idxBin]; ++idxMass) {
				for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
					bool zeroThisWave = true;
					zeroThisWave &= (_productionAmplitudes[idxBin][idxMass][idxWave].real() == 0.);
					zeroThisWave &= (_productionAmplitudes[idxBin][idxMass][idxWave].imag() == 0.);
					zeroThisWave &= (productionAmplitudesCovariance[idxBin][idxMass](2*idxWave,   2*idxWave  ) == 0.);
					zeroThisWave &= (productionAmplitudesCovariance[idxBin][idxMass](2*idxWave,   2*idxWave+1) == 0.);
					zeroThisWave &= (productionAmplitudesCovariance[idxBin][idxMass](2*idxWave+1, 2*idxWave  ) == 0.);
					zeroThisWave &= (productionAmplitudesCovariance[idxBin][idxMass](2*idxWave+1, 2*idxWave+1) == 0.);

					if(zeroThisWave or idxMass < _wavePairMassBinLimits[idxBin][idxWave][idxWave].first or idxMass > _wavePairMassBinLimits[idxBin][idxWave][idxWave].second) {
						zeroWaves[idxBin][idxMass].push_back(idxWave);
					}

					if(idxWave == _idxAnchorWave) {
						zeroAnchorWave |= zeroThisWave;
					}

					// check that a wave is not zero in its fit range
					if(zeroThisWave and idxMass >= _wavePairMassBinLimits[idxBin][idxWave][idxWave].first and idxMass <= _wavePairMassBinLimits[idxBin][idxWave][idxWave].second) {
						printErr << "production amplitudes of wave " << idxWave << " zero in its fit range (e.g. mass limit in mass-independent fit)." << std::endl;
						return false;
					}
				}
			}
		}

		// error if anchor wave is zero in one mass bin
		if(zeroAnchorWave) {
			printErr << "production amplitudes of anchor wave zero in some mass bins (mass limit in mass-independent fit)." << std::endl;
			return false;
		}

		// test if the anchor wave is real valued
		// test that any non-anchor wave is not real valued
		bool realAnchorWave = true;
		bool realOtherWaves = false;
		for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
			bool realThisWave = true;
			for(size_t idxBin=0; idxBin<_nrBins; ++idxBin) {
				for(size_t idxMass = _idxMassMin[idxBin]; idxMass <= _idxMassMax[idxBin]; ++idxMass) {
					realThisWave &= (_productionAmplitudes[idxBin][idxMass][idxWave].imag() == 0.);
					realThisWave &= (productionAmplitudesCovariance[idxBin][idxMass](2*idxWave,   2*idxWave+1) == 0.);
					realThisWave &= (productionAmplitudesCovariance[idxBin][idxMass](2*idxWave+1, 2*idxWave  ) == 0.);
					realThisWave &= (productionAmplitudesCovariance[idxBin][idxMass](2*idxWave+1, 2*idxWave+1) == 0.);
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

		_productionAmplitudesCovMatInv.resize(boost::extents[_nrBins][_maxMassBins]);
		for(size_t idxBin=0; idxBin<_nrBins; ++idxBin) {
			for(size_t idxMass = _idxMassMin[idxBin]; idxMass <= _idxMassMax[idxBin]; ++idxMass) {
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

							reducedCovMat(rowSkip, colSkip) = productionAmplitudesCovariance[idxBin][idxMass](2*idxWave, 2*jdxWave);
							if(jdxWave != _idxAnchorWave) {
								reducedCovMat(rowSkip, colSkip+1) = productionAmplitudesCovariance[idxBin][idxMass](2*idxWave, 2*jdxWave+1);
							}
							if(idxWave != _idxAnchorWave) {
								reducedCovMat(rowSkip+1, colSkip) = productionAmplitudesCovariance[idxBin][idxMass](2*idxWave+1, 2*jdxWave);
							}
							if(idxWave != _idxAnchorWave && jdxWave != _idxAnchorWave) {
								reducedCovMat(rowSkip+1, colSkip+1) = productionAmplitudesCovariance[idxBin][idxMass](2*idxWave+1, 2*jdxWave+1);
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

							covariance(2*(idxWave-idxSkip),   2*(jdxWave-jdxSkip)  ) = productionAmplitudesCovariance[idxBin][idxMass](2*idxWave,   2*jdxWave  );
							covariance(2*(idxWave-idxSkip),   2*(jdxWave-jdxSkip)+1) = productionAmplitudesCovariance[idxBin][idxMass](2*idxWave,   2*jdxWave+1);
							covariance(2*(idxWave-idxSkip)+1, 2*(jdxWave-jdxSkip)  ) = productionAmplitudesCovariance[idxBin][idxMass](2*idxWave+1, 2*jdxWave  );
							covariance(2*(idxWave-idxSkip)+1, 2*(jdxWave-jdxSkip)+1) = productionAmplitudesCovariance[idxBin][idxMass](2*idxWave+1, 2*jdxWave+1);

							const double n = abs(_productionAmplitudes[idxBin][idxMass][_idxAnchorWave]);
							const double n3 = std::pow(n, 3);
							const double xa1 = _productionAmplitudes[idxBin][idxMass][_idxAnchorWave].real();
							const double xa2 = _productionAmplitudes[idxBin][idxMass][_idxAnchorWave].imag();
							const double xi1 = _productionAmplitudes[idxBin][idxMass][idxWave].real();
							const double xi2 = _productionAmplitudes[idxBin][idxMass][idxWave].imag();

							const Int_t rowSkip = 2*(idxWave-idxSkip) + ((idxWave>_idxAnchorWave && !skipAnchor) ? -1 : 0);
							if(idxWave == _idxAnchorWave && jdxWave == _idxAnchorWave) {
								jacobian(rowSkip, 2*(jdxWave-jdxSkip)  ) = xa1 / n;
								jacobian(rowSkip, 2*(jdxWave-jdxSkip)+1) = xa2 / n;
							} else if(jdxWave == _idxAnchorWave) {
								jacobian(rowSkip,   2*(jdxWave-jdxSkip)  ) =   xi1 / n - xa1 * (xi1*xa1 + xi2*xa2) / n3;
								jacobian(rowSkip,   2*(jdxWave-jdxSkip)+1) =   xi2 / n - xa2 * (xi2*xa2 + xi1*xa1) / n3;
								jacobian(rowSkip+1, 2*(jdxWave-jdxSkip)  ) =   xi2 / n - xa1 * (xi2*xa1 - xi1*xa2) / n3;
								jacobian(rowSkip+1, 2*(jdxWave-jdxSkip)+1) = - xi1 / n - xa2 * (xi2*xa1 - xi1*xa2) / n3;
							} else if(idxWave == jdxWave) {
								jacobian(rowSkip  , 2*(jdxWave-jdxSkip)  ) =   xa1 / n;
								jacobian(rowSkip  , 2*(jdxWave-jdxSkip)+1) =   xa2 / n;
								jacobian(rowSkip+1, 2*(jdxWave-jdxSkip)  ) = - xa2 / n;
								jacobian(rowSkip+1, 2*(jdxWave-jdxSkip)+1) =   xa1 / n;
							}
						}
					}

					TMatrixT<double> jacobianT(TMatrixT<double>::kTransposed, jacobian);

					reducedCovMat = jacobian * covariance * jacobianT;
				}

				// set entries in covariance matrix to zero according to which parts are to be used
				if(_useCovariance != useFullCovarianceMatrix) {
					for(size_t idxWave = 0, idxSkip = 0; idxWave < _nrWaves; ++idxWave) {
						if(idxSkip < zeroWaves[idxBin][idxMass].size() and zeroWaves[idxBin][idxMass][idxSkip] == idxWave) {
							++idxSkip;
							continue;
						}

						for(size_t jdxWave = 0, jdxSkip = 0; jdxWave < _nrWaves; ++jdxWave) {
							if(jdxSkip < zeroWaves[idxBin][idxMass].size() and zeroWaves[idxBin][idxMass][jdxSkip] == jdxWave) {
								++jdxSkip;
								continue;
							}

							const Int_t rowSkip = 2*(idxWave-idxSkip) + ((idxWave>_idxAnchorWave and not skipAnchor) ? -1 : 0);
							const Int_t colSkip = 2*(jdxWave-jdxSkip) + ((jdxWave>_idxAnchorWave and not skipAnchor) ? -1 : 0);

							if(idxWave == jdxWave) {
								if(_useCovariance == useDiagnalElementsOnly) {
									if(jdxWave != _idxAnchorWave) {
										reducedCovMat(rowSkip, colSkip+1) = 0;
									}
									if(idxWave != _idxAnchorWave) {
										reducedCovMat(rowSkip+1, colSkip) = 0;
									}
								}
							} else {
								reducedCovMat(rowSkip, colSkip) = 0;
								if(jdxWave != _idxAnchorWave) {
									reducedCovMat(rowSkip, colSkip+1) = 0;
								}
								if(idxWave != _idxAnchorWave) {
									reducedCovMat(rowSkip+1, colSkip) = 0;
								}
								if(idxWave != _idxAnchorWave and jdxWave != _idxAnchorWave) {
									reducedCovMat(rowSkip+1, colSkip+1) = 0;
								}
							}
						}
					}
				}

				reducedCovMat.Invert();

				// import covariance matrix of production amplitudes
				_productionAmplitudesCovMatInv[idxBin][idxMass].ResizeTo(2*_nrWaves, 2*_nrWaves);

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

						_productionAmplitudesCovMatInv[idxBin][idxMass](2*idxWave, 2*jdxWave) = reducedCovMat(rowSkip, colSkip);
						if(jdxWave != _idxAnchorWave) {
							_productionAmplitudesCovMatInv[idxBin][idxMass](2*idxWave, 2*jdxWave+1) = reducedCovMat(rowSkip, colSkip+1);
						}
						if(idxWave != _idxAnchorWave) {
							_productionAmplitudesCovMatInv[idxBin][idxMass](2*idxWave+1, 2*jdxWave) = reducedCovMat(rowSkip+1, colSkip);
						}
						if(idxWave != _idxAnchorWave && jdxWave != _idxAnchorWave) {
							_productionAmplitudesCovMatInv[idxBin][idxMass](2*idxWave+1, 2*jdxWave+1) = reducedCovMat(rowSkip+1, colSkip+1);
						}
					}
				}
			}
		}

		// modify measured production amplitude such that the anchor wave is always real and positive
		for(size_t idxBin=0; idxBin<_nrBins; ++idxBin) {
			for(size_t idxMass = _idxMassMin[idxBin]; idxMass <= _idxMassMax[idxBin]; ++idxMass) {
				const std::complex<double> anchorPhase = _productionAmplitudes[idxBin][idxMass][_idxAnchorWave] / abs(_productionAmplitudes[idxBin][idxMass][_idxAnchorWave]);
				for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
					_productionAmplitudes[idxBin][idxMass][idxWave] /= anchorPhase;
				}
			}
		}
	} else {
		// do some stuff specific to the fit to the spin-density matrix

		// get a list of waves that are zero (those have to be excluded
		// from the inversion of the covariance matrix below)
		boost::multi_array<std::vector<size_t>, 2> zeroWaves(boost::extents[_nrBins][_maxMassBins]);
		for(size_t idxBin=0; idxBin<_nrBins; ++idxBin) {
			for(size_t idxMass = _idxMassMin[idxBin]; idxMass <= _idxMassMax[idxBin]; ++idxMass) {
				for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
					bool zeroThisWave = true;
					for(size_t jdxWave=0; jdxWave<_nrWaves; ++jdxWave) {
						const size_t idx = _nrWaves*(_nrWaves+1) - ((jdxWave >= idxWave) ? ((_nrWaves-idxWave)*(_nrWaves-idxWave+1) - 2*(jdxWave-idxWave)) : ((_nrWaves-jdxWave)*(_nrWaves-jdxWave+1) - 2*(idxWave-jdxWave)));

						zeroThisWave &= (_spinDensityMatrices[idxBin][idxMass][idxWave][jdxWave].real() == 0.);
						zeroThisWave &= (_spinDensityMatrices[idxBin][idxMass][idxWave][jdxWave].imag() == 0.);
						zeroThisWave &= (spinDensityCovarianceMatrices[idxBin][idxMass](idx,   idx  ) == 0.);
						zeroThisWave &= (spinDensityCovarianceMatrices[idxBin][idxMass](idx,   idx+1) == 0.);
						zeroThisWave &= (spinDensityCovarianceMatrices[idxBin][idxMass](idx+1, idx  ) == 0.);
						zeroThisWave &= (spinDensityCovarianceMatrices[idxBin][idxMass](idx+1, idx+1) == 0.);
					}

					if(zeroThisWave or idxMass < _wavePairMassBinLimits[idxBin][idxWave][idxWave].first or idxMass > _wavePairMassBinLimits[idxBin][idxWave][idxWave].second) {
						zeroWaves[idxBin][idxMass].push_back(idxWave);
					}

					// check that a wave is not zero in its fit range
					if(zeroThisWave and idxMass >= _wavePairMassBinLimits[idxBin][idxWave][idxWave].first and idxMass <= _wavePairMassBinLimits[idxBin][idxWave][idxWave].second) {
						printErr << "spin-density matrix element of wave " << idxWave << " zero in its fit range (e.g. mass limit in mass-independent fit)." << std::endl;
						return false;
					}
				}
			}
		}

		if(_useCovariance == useFullCovarianceMatrix) {
			_spinDensityMatricesCovMatInv.resize(boost::extents[_nrBins][_maxMassBins]);
		} else {
			_spinDensityMatricesCovMatInvArray.resize(boost::extents[_nrBins][_maxMassBins][_nrWaves][_nrWaves][2][2]);
		}
		for(size_t idxBin=0; idxBin<_nrBins; ++idxBin) {
			for(size_t idxMass = _idxMassMin[idxBin]; idxMass <= _idxMassMax[idxBin]; ++idxMass) {
				// import covariance matrix of spin-density matrix elements
				const size_t matrixSize(_nrWaves * (_nrWaves+1));
				const size_t reducedMatrixSize((_nrWaves - zeroWaves[idxBin][idxMass].size())*(_nrWaves - zeroWaves[idxBin][idxMass].size()));
				TMatrixT<double> reducedCovMat(reducedMatrixSize, reducedMatrixSize);

				{
					// i is for loop over rows
					size_t redIdx = 0;
					for(size_t iWave1 = 0, iSkip1 = 0; iWave1 < _nrWaves; ++iWave1) {
						if(iSkip1 < zeroWaves[idxBin][idxMass].size() and zeroWaves[idxBin][idxMass][iSkip1] == iWave1) {
							++iSkip1;
							continue;
						}
						for(size_t iWave2 = iWave1, iSkip2 = iSkip1; iWave2 < _nrWaves; ++iWave2) {
							if(iSkip2 < zeroWaves[idxBin][idxMass].size() and zeroWaves[idxBin][idxMass][iSkip2] == iWave2) {
								++iSkip2;
								continue;
							}
							const size_t idx = _nrWaves*(_nrWaves+1) - (_nrWaves-iWave1)*(_nrWaves-iWave1+1) + 2*(iWave2-iWave1);

							// j is for loop over columns
							size_t redJdx = 0;
							for(size_t jWave1 = 0, jSkip1 = 0; jWave1 < _nrWaves; ++jWave1) {
								if(jSkip1 < zeroWaves[idxBin][idxMass].size() and zeroWaves[idxBin][idxMass][jSkip1] == jWave1) {
									++jSkip1;
									continue;
								}
								for(size_t jWave2 = jWave1, jSkip2 = jSkip1; jWave2 < _nrWaves; ++jWave2) {
									if(jSkip2 < zeroWaves[idxBin][idxMass].size() and zeroWaves[idxBin][idxMass][jSkip2] == jWave2) {
										++jSkip2;
										continue;
									}
									const size_t jdx = _nrWaves*(_nrWaves+1) - (_nrWaves-jWave1)*(_nrWaves-jWave1+1) + 2*(jWave2-jWave1);

									if(iWave1 == iWave2) { // one row
										if(jWave1 == jWave2) { // one column
											if((iWave1 == jWave1 and iWave2 == jWave2) or _useCovariance == useFullCovarianceMatrix) {
												reducedCovMat(redIdx,   redJdx  ) = spinDensityCovarianceMatrices[idxBin][idxMass](idx,   jdx  );
											}
										} else { // two columns
											if(_useCovariance == useFullCovarianceMatrix) {
												reducedCovMat(redIdx,   redJdx  ) = spinDensityCovarianceMatrices[idxBin][idxMass](idx,   jdx  );
												reducedCovMat(redIdx,   redJdx+1) = spinDensityCovarianceMatrices[idxBin][idxMass](idx,   jdx+1);
											}
										}
									} else { // two rows
										if(jWave1 == jWave2) { // one column
											if(_useCovariance == useFullCovarianceMatrix) {
												reducedCovMat(redIdx,   redJdx  ) = spinDensityCovarianceMatrices[idxBin][idxMass](idx,   jdx  );
												reducedCovMat(redIdx+1, redJdx  ) = spinDensityCovarianceMatrices[idxBin][idxMass](idx+1, jdx  );
											}
										} else { // two columns
											if((iWave1 == jWave1 and iWave2 == jWave2) or _useCovariance == useFullCovarianceMatrix) {
												reducedCovMat(redIdx,   redJdx  ) = spinDensityCovarianceMatrices[idxBin][idxMass](idx,   jdx  );
												reducedCovMat(redIdx+1, redJdx+1) = spinDensityCovarianceMatrices[idxBin][idxMass](idx+1, jdx+1);
												if(_useCovariance != useDiagnalElementsOnly) {
													reducedCovMat(redIdx,   redJdx+1) = spinDensityCovarianceMatrices[idxBin][idxMass](idx,   jdx+1);
													reducedCovMat(redIdx+1, redJdx  ) = spinDensityCovarianceMatrices[idxBin][idxMass](idx+1, jdx  );
												}
											}
										}
									}

									if(jWave1 == jWave2) {
										redJdx += 1;
									} else {
										redJdx += 2;
									}
								}
							}

							if(iWave1 == iWave2) {
								redIdx += 1;
							} else {
								redIdx += 2;
							}
						}
					}
				}

				reducedCovMat.Invert();

				// import covariance matrix of spin-density matrix elements
				if(_useCovariance == useFullCovarianceMatrix) {
					_spinDensityMatricesCovMatInv[idxBin][idxMass].ResizeTo(matrixSize, matrixSize);

					// i is for loop over rows
					size_t redIdx = 0;
					for(size_t iWave1 = 0, iSkip1 = 0; iWave1 < _nrWaves; ++iWave1) {
						if(iSkip1 < zeroWaves[idxBin][idxMass].size() and zeroWaves[idxBin][idxMass][iSkip1] == iWave1) {
							++iSkip1;
							continue;
						}
						for(size_t iWave2 = iWave1, iSkip2 = iSkip1; iWave2 < _nrWaves; ++iWave2) {
							if(iSkip2 < zeroWaves[idxBin][idxMass].size() and zeroWaves[idxBin][idxMass][iSkip2] == iWave2) {
								++iSkip2;
								continue;
							}
							const size_t idx = _nrWaves*(_nrWaves+1) - (_nrWaves-iWave1)*(_nrWaves-iWave1+1) + 2*(iWave2-iWave1);

							// j is for loop over columns
							size_t redJdx = 0;
							for(size_t jWave1 = 0, jSkip1 = 0; jWave1 < _nrWaves; ++jWave1) {
								if(jSkip1 < zeroWaves[idxBin][idxMass].size() and zeroWaves[idxBin][idxMass][jSkip1] == jWave1) {
									++jSkip1;
									continue;
								}
								for(size_t jWave2 = jWave1, jSkip2 = jSkip1; jWave2 < _nrWaves; ++jWave2) {
									if(jSkip2 < zeroWaves[idxBin][idxMass].size() and zeroWaves[idxBin][idxMass][jSkip2] == jWave2) {
										++jSkip2;
										continue;
									}
									const size_t jdx = _nrWaves*(_nrWaves+1) - (_nrWaves-jWave1)*(_nrWaves-jWave1+1) + 2*(jWave2-jWave1);

									if(iWave1 == iWave2) { // one row
										if(jWave1 == jWave2) { // one column
											_spinDensityMatricesCovMatInv[idxBin][idxMass](idx,   jdx  ) = reducedCovMat(redIdx,   redJdx  );
										} else { // two columns
											_spinDensityMatricesCovMatInv[idxBin][idxMass](idx,   jdx  ) = reducedCovMat(redIdx,   redJdx  );
											_spinDensityMatricesCovMatInv[idxBin][idxMass](idx,   jdx+1) = reducedCovMat(redIdx,   redJdx+1);
										}
									} else { // two rows
										if(jWave1 == jWave2) { // one column
											_spinDensityMatricesCovMatInv[idxBin][idxMass](idx,   jdx  ) = reducedCovMat(redIdx,   redJdx  );
											_spinDensityMatricesCovMatInv[idxBin][idxMass](idx+1, jdx  ) = reducedCovMat(redIdx+1, redJdx  );
										} else { // two columns
											_spinDensityMatricesCovMatInv[idxBin][idxMass](idx,   jdx  ) = reducedCovMat(redIdx,   redJdx  );
											_spinDensityMatricesCovMatInv[idxBin][idxMass](idx,   jdx+1) = reducedCovMat(redIdx,   redJdx+1);
											_spinDensityMatricesCovMatInv[idxBin][idxMass](idx+1, jdx  ) = reducedCovMat(redIdx+1, redJdx  );
											_spinDensityMatricesCovMatInv[idxBin][idxMass](idx+1, jdx+1) = reducedCovMat(redIdx+1, redJdx+1);
										}
									}

									if(jWave1 == jWave2) {
										redJdx += 1;
									} else {
										redJdx += 2;
									}
								}
							}

							if(iWave1 == iWave2) {
								redIdx += 1;
							} else {
								redIdx += 2;
							}
						}
					}
				}else {
					// i is for loop over rows
					size_t redIdx = 0;
					for(size_t iWave = 0, iSkip = 0; iWave < _nrWaves; ++iWave) {
						if(iSkip < zeroWaves[idxBin][idxMass].size() and zeroWaves[idxBin][idxMass][iSkip] == iWave) {
							++iSkip;
							continue;
						}
						for(size_t jWave = iWave, jSkip = iSkip; jWave < _nrWaves; ++jWave) {
							if(jSkip < zeroWaves[idxBin][idxMass].size() and zeroWaves[idxBin][idxMass][jSkip] == jWave) {
								++jSkip;
								continue;
							}

							if(iWave == jWave) {
								_spinDensityMatricesCovMatInvArray[idxBin][idxMass][iWave][jWave][0][0] = reducedCovMat(redIdx,   redIdx  );
							} else {
								_spinDensityMatricesCovMatInvArray[idxBin][idxMass][iWave][jWave][0][0] = reducedCovMat(redIdx,   redIdx  );
								_spinDensityMatricesCovMatInvArray[idxBin][idxMass][iWave][jWave][0][1] = reducedCovMat(redIdx,   redIdx+1);
								_spinDensityMatricesCovMatInvArray[idxBin][idxMass][iWave][jWave][1][0] = reducedCovMat(redIdx+1, redIdx  );
								_spinDensityMatricesCovMatInvArray[idxBin][idxMass][iWave][jWave][1][1] = reducedCovMat(redIdx+1, redIdx+1);
							}

							if(iWave == jWave) {
								redIdx += 1;
							} else {
								redIdx += 2;
							}
						}
					}
				}
			}
		}
	}

	return true;
}


size_t
rpwa::resonanceFit::function::getNrParameters() const
{
	return _fitModel->getNrParameters();
}


size_t
rpwa::resonanceFit::function::getNrDataPoints() const
{
	size_t nrPts(0);

	if(_fitProductionAmplitudes) {
		// calculate data points:
		// * production amplitudes in general are complex numbers
		// * for the anchor wave it might be real
		// * remember (Re,Im) => factor 2
		for(size_t idxBin = 0; idxBin < _nrBins; ++idxBin) {
			for(size_t idxWave = 0; idxWave < _nrWaves; ++idxWave) {
				nrPts += _wavePairMassBinLimits[idxBin][idxWave][idxWave].second - _wavePairMassBinLimits[idxBin][idxWave][idxWave].first + 1;
				if(idxWave != _idxAnchorWave) {
					nrPts += _wavePairMassBinLimits[idxBin][idxWave][idxWave].second - _wavePairMassBinLimits[idxBin][idxWave][idxWave].first + 1;
				}
			}
		}
	} else {
		// calculate data points:
		// * diagonal elements are real numbers
		// * off-diagonal elements are complex numbers
		// * remember (Re,Im) => factor 2
		// * diagonal elements are only checked once, off-diagonal elements with
		//   the two different combinations (i,j) and (j,i)
		for(size_t idxBin = 0; idxBin < _nrBins; ++idxBin) {
			for(size_t idxWave = 0; idxWave < _nrWaves; ++idxWave) {
				for(size_t jdxWave = 0; jdxWave < _nrWaves; ++jdxWave) {
					nrPts += _wavePairMassBinLimits[idxBin][idxWave][jdxWave].second - _wavePairMassBinLimits[idxBin][idxWave][jdxWave].first + 1;
				}
			}
		}
	}

	return nrPts;
}


double
rpwa::resonanceFit::function::chiSquare(const std::vector<double>& par) const
{
	return chiSquare(par.data());
}


double
rpwa::resonanceFit::function::chiSquare(const double* par) const
{
	// in C++11 we can use a static variable per thread so that the
	// parameters are kept over function calls and we can implement some
	// caching
	thread_local rpwa::resonanceFit::parameters fitParameters(_fitModel->getNrComponents()+1,            // nr components + final-state mass-dependence
	                                                          _fitModel->getMaxChannelsInComponent(),
	                                                          _fitModel->getMaxParametersInComponent(),
	                                                          _nrBins);
	thread_local rpwa::resonanceFit::cache cache(_nrWaves,
	                                             _fitModel->getNrComponents()+1,          // nr components + final-state mass-dependence
	                                             _fitModel->getMaxChannelsInComponent(),
	                                             _nrBins,
	                                             _maxMassBins);

	// import parameters (couplings, branchings, resonance parameters, ...)
	_fitModel->importParameters(par, fitParameters, cache);

	return chiSquare(fitParameters, cache);
}


double
rpwa::resonanceFit::function::chiSquare(const rpwa::resonanceFit::parameters& fitParameters,
                                        rpwa::resonanceFit::cache& cache) const
{
	if(_fitProductionAmplitudes) {
		return chiSquareProductionAmplitudes(fitParameters, cache);
	} else {
		return chiSquareSpinDensityMatrix(fitParameters, cache);
	}
}


double
rpwa::resonanceFit::function::logLikelihood(const std::vector<double>& par) const
{
	return -0.5 * chiSquare(par);
}


double
rpwa::resonanceFit::function::logLikelihood(const double* par) const
{
	return -0.5 * chiSquare(par);
}


double
rpwa::resonanceFit::function::logLikelihood(const rpwa::resonanceFit::parameters& fitParameters,
                                            rpwa::resonanceFit::cache& cache) const
{
	return -0.5 * chiSquare(fitParameters, cache);
}


double
rpwa::resonanceFit::function::logPriorLikelihood(const std::vector<double>& par) const
{
	return logPriorLikelihood(par.data());
}


double
rpwa::resonanceFit::function::logPriorLikelihood(const double* par) const
{
	rpwa::resonanceFit::parameters fitParameters(_fitModel->getNrComponents()+1,            // nr components + final-state mass-dependence
	                                             _fitModel->getMaxChannelsInComponent(),
	                                             _fitModel->getMaxParametersInComponent(),
	                                             _nrBins);
	rpwa::resonanceFit::cache cache(_nrWaves,
	                                _fitModel->getNrComponents()+1,          // nr components + final-state mass-dependence
	                                _fitModel->getMaxChannelsInComponent(),
	                                _nrBins,
	                                _maxMassBins);

	// import parameters (couplings, branchings, resonance parameters, ...)
	_fitModel->importParameters(par, fitParameters, cache);

	return logPriorLikelihood(fitParameters);
}


double
rpwa::resonanceFit::function::logPriorLikelihood(const rpwa::resonanceFit::parameters& fitParameters) const
{
	double logPrior = 0;

	const size_t nrComponents = _fitModel->getNrComponents();
	for (size_t idxComponent = 0; idxComponent < nrComponents; ++idxComponent) {
		const rpwa::resonanceFit::componentConstPtr& component = _fitModel->getComponent(idxComponent);

		const size_t nrParameters = component->getNrParameters();
		for (size_t idxParameter = 0; idxParameter < nrParameters; ++idxParameter) {
			// fixed parameters to no contribute to prior likelihood
			if (component->getParameterFixed(idxParameter))
				continue;

			// parameters with 0 error are assumed to have a flat prior
			if (component->getParameterError(idxParameter) == 0.0)
				continue;

			logPrior += -0.5 * std::pow((fitParameters.getParameter(idxComponent, idxParameter) - component->getParameterStart(idxParameter)) / component->getParameterError(idxParameter), 2.);
		}
	}

	return logPrior;
}


double
rpwa::resonanceFit::function::chiSquareProductionAmplitudes(const rpwa::resonanceFit::parameters& fitParameters,
                                                            rpwa::resonanceFit::cache& cache) const
{
	double chi2=0;

	// loop over bins
	for(unsigned idxBin=0; idxBin<_nrBins; ++idxBin) {
		// loop over mass-bins
		for(unsigned idxMass = _idxMassMin[idxBin]; idxMass <= _idxMassMax[idxBin]; ++idxMass) {
			const double mass = _massBinCenters[idxBin][idxMass];

			// phase of fit in anchor wave
			const std::complex<double> anchorFit = _fitModel->productionAmplitude(fitParameters, cache, _idxAnchorWave, idxBin, mass, idxMass);
			const std::complex<double> anchorFitPhase = anchorFit / abs(anchorFit);

			TVectorT<double> prodAmpDiffVect(2*_nrWaves);

			// sum over the contributions to chi2
			for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
				// check that this mass bin should be taken into account for this
				// combination of waves
				if(idxMass < _wavePairMassBinLimits[idxBin][idxWave][idxWave].first or idxMass > _wavePairMassBinLimits[idxBin][idxWave][idxWave].second) {
					continue;
				}

				// calculate target spin density matrix element
				const std::complex<double> prodAmpFit = _fitModel->productionAmplitude(fitParameters, cache, idxWave, idxBin, mass, idxMass) / anchorFitPhase;

				const std::complex<double> prodAmpDiff = prodAmpFit - _productionAmplitudes[idxBin][idxMass][idxWave];

				const Int_t row = 2*idxWave;
				prodAmpDiffVect(row) = prodAmpDiff.real();
				if(idxWave != _idxAnchorWave) {
					prodAmpDiffVect(row+1) = prodAmpDiff.imag();
				}
			} // end loop over idxWave

			chi2 += _productionAmplitudesCovMatInv[idxBin][idxMass].Similarity(prodAmpDiffVect);
		} // end loop over mass-bins
	} // end loop over bins

	return chi2;
}


double
rpwa::resonanceFit::function::chiSquareSpinDensityMatrix(const rpwa::resonanceFit::parameters& fitParameters,
                                                         rpwa::resonanceFit::cache& cache) const
{
	double chi2=0;

	// loop over bins
	for(unsigned idxBin=0; idxBin<_nrBins; ++idxBin) {
		// loop over mass-bins
		for(unsigned idxMass = _idxMassMin[idxBin]; idxMass <= _idxMassMax[idxBin]; ++idxMass) {
			const double mass = _massBinCenters[idxBin][idxMass];

			// sum over the contributions to chi2 -> rho_ij
			for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
				for(size_t jdxWave=idxWave; jdxWave<_nrWaves; ++jdxWave) {
					// check that this mass bin should be taken into account for this
					// combination of waves
					if(idxMass < _wavePairMassBinLimits[idxBin][idxWave][jdxWave].first or idxMass > _wavePairMassBinLimits[idxBin][idxWave][jdxWave].second) {
						continue;
					}

					// calculate target spin density matrix element
					const std::complex<double> rhoFit = _fitModel->spinDensityMatrix(fitParameters, cache, idxWave, jdxWave, idxBin, mass, idxMass);

					const std::complex<double> rhoDiff = rhoFit - _spinDensityMatrices[idxBin][idxMass][idxWave][jdxWave];

					if(idxWave==jdxWave) {
						chi2 += rhoDiff.real() * _spinDensityMatricesCovMatInvArray[idxBin][idxMass][idxWave][jdxWave][0][0] * rhoDiff.real();
					} else {
						if (_useCovariance == useDiagnalElementsOnly) {
							chi2 += rhoDiff.real() * _spinDensityMatricesCovMatInvArray[idxBin][idxMass][idxWave][jdxWave][0][0] * rhoDiff.real();
							chi2 += rhoDiff.imag() * _spinDensityMatricesCovMatInvArray[idxBin][idxMass][idxWave][jdxWave][1][1] * rhoDiff.imag();
						} else if(_useCovariance == useComplexDiagnalElementsOnly) {
							chi2 += rhoDiff.real() * _spinDensityMatricesCovMatInvArray[idxBin][idxMass][idxWave][jdxWave][0][0] * rhoDiff.real();
							chi2 += rhoDiff.real() * _spinDensityMatricesCovMatInvArray[idxBin][idxMass][idxWave][jdxWave][0][1] * rhoDiff.imag();
							chi2 += rhoDiff.imag() * _spinDensityMatricesCovMatInvArray[idxBin][idxMass][idxWave][jdxWave][1][0] * rhoDiff.real();
							chi2 += rhoDiff.imag() * _spinDensityMatricesCovMatInvArray[idxBin][idxMass][idxWave][jdxWave][1][1] * rhoDiff.imag();
						} else {
							// this should have returned an error during the call to init()
							assert(false);
						}
					}
				} // end loop over jdxWave
			} // end loop over idxWave
		} // end loop over mass-bins
	} // end loop over bins

	return chi2;
}
