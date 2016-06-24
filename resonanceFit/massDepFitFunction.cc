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


#include "massDepFitFunction.h"

#include "massDepFitCache.h"
#include "massDepFitComponents.h"
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
                                 const boost::multi_array<TMatrixT<double>, 2>& productionAmplitudesCovariance,
                                 const boost::multi_array<std::complex<double>, 4>& spinDensityMatrices,
                                 const boost::multi_array<TMatrixT<double>, 2>& spinDensityCovarianceMatrices,
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

	_spinDensityMatrices.resize(std::vector<size_t>(spinDensityMatrices.shape(), spinDensityMatrices.shape()+spinDensityMatrices.num_dimensions()));
	_spinDensityMatrices = spinDensityMatrices;

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
					zeroThisWave &= (productionAmplitudesCovariance[idxBin][idxMass](2*idxWave,   2*idxWave  ) == 0.);
					zeroThisWave &= (productionAmplitudesCovariance[idxBin][idxMass](2*idxWave,   2*idxWave+1) == 0.);
					zeroThisWave &= (productionAmplitudesCovariance[idxBin][idxMass](2*idxWave+1, 2*idxWave  ) == 0.);
					zeroThisWave &= (productionAmplitudesCovariance[idxBin][idxMass](2*idxWave+1, 2*idxWave+1) == 0.);

					if(zeroThisWave || idxMass < _wavePairMassBinLimits[idxWave][idxWave].first || idxMass > _wavePairMassBinLimits[idxWave][idxWave].second) {
						zeroWaves[idxBin][idxMass].push_back(idxWave);
					}

					if(idxWave == _idxAnchorWave) {
						zeroAnchorWave |= zeroThisWave;
					}

					// check that a wave is not zero in its fit range
					if(zeroThisWave && idxMass >= _wavePairMassBinLimits[idxWave][idxWave].first && idxMass <= _wavePairMassBinLimits[idxWave][idxWave].second) {
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
				for(size_t idxMass=0; idxMass<_nrMassBins; ++idxMass) {
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

							reducedCovMat(rowSkip + 0, colSkip + 0) = productionAmplitudesCovariance[idxBin][idxMass](2*idxWave, 2*jdxWave);
							if(jdxWave != _idxAnchorWave) {
								reducedCovMat(rowSkip + 0, colSkip + 1) = productionAmplitudesCovariance[idxBin][idxMass](2*idxWave, 2*jdxWave+1);
							}
							if(idxWave != _idxAnchorWave) {
								reducedCovMat(rowSkip + 1, colSkip + 0) = productionAmplitudesCovariance[idxBin][idxMass](2*idxWave+1, 2*jdxWave);
							}
							if(idxWave != _idxAnchorWave && jdxWave != _idxAnchorWave) {
								reducedCovMat(rowSkip + 1, colSkip + 1) = productionAmplitudesCovariance[idxBin][idxMass](2*idxWave+1, 2*jdxWave+1);
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

							covariance(2*(idxWave-idxSkip) + 0, 2*(jdxWave-jdxSkip) + 0) = productionAmplitudesCovariance[idxBin][idxMass](2*idxWave,   2*jdxWave  );
							covariance(2*(idxWave-idxSkip) + 0, 2*(jdxWave-jdxSkip) + 1) = productionAmplitudesCovariance[idxBin][idxMass](2*idxWave,   2*jdxWave+1);
							covariance(2*(idxWave-idxSkip) + 1, 2*(jdxWave-jdxSkip) + 0) = productionAmplitudesCovariance[idxBin][idxMass](2*idxWave+1, 2*jdxWave  );
							covariance(2*(idxWave-idxSkip) + 1, 2*(jdxWave-jdxSkip) + 1) = productionAmplitudesCovariance[idxBin][idxMass](2*idxWave+1, 2*jdxWave+1);

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
						const size_t idx = _nrWaves*(_nrWaves+1) - ((jdxWave >= idxWave) ? ((_nrWaves-idxWave)*(_nrWaves-idxWave+1) - 2*(jdxWave-idxWave)) : ((_nrWaves-jdxWave)*(_nrWaves-jdxWave+1) - 2*(idxWave-jdxWave)));

						zeroThisWave &= (_spinDensityMatrices[idxBin][idxMass][idxWave][jdxWave].real() == 0.);
						zeroThisWave &= (_spinDensityMatrices[idxBin][idxMass][idxWave][jdxWave].imag() == 0.);
						zeroThisWave &= (spinDensityCovarianceMatrices[idxBin][idxMass](idx,   idx  ) == 0.);
						zeroThisWave &= (spinDensityCovarianceMatrices[idxBin][idxMass](idx,   idx+1) == 0.);
						zeroThisWave &= (spinDensityCovarianceMatrices[idxBin][idxMass](idx+1, idx  ) == 0.);
						zeroThisWave &= (spinDensityCovarianceMatrices[idxBin][idxMass](idx+1, idx+1) == 0.);
					}

					if(zeroThisWave || idxMass < _wavePairMassBinLimits[idxWave][idxWave].first || idxMass > _wavePairMassBinLimits[idxWave][idxWave].second) {
						zeroWaves[idxBin][idxMass].push_back(idxWave);
					}

					// check that a wave is not zero in its fit range
					if(zeroThisWave && idxMass >= _wavePairMassBinLimits[idxWave][idxWave].first && idxMass <= _wavePairMassBinLimits[idxWave][idxWave].second) {
						printErr << "spin-density matrix element of wave " << idxWave << " zero in its fit range (e.g. mass limit in mass-independent fit)." << std::endl;
						return false;
					}
				}
			}
		}

		_spinDensityMatricesCovMatInv.resize(boost::extents[_nrBins][_nrMassBins]);
		for(size_t idxBin=0; idxBin<_nrBins; ++idxBin) {
			for(size_t idxMass=_idxMassMin; idxMass<=_idxMassMax; ++idxMass) {
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
				_spinDensityMatricesCovMatInv[idxBin][idxMass].ResizeTo(matrixSize, matrixSize);
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
	return chiSquare(par.data());
}


double
rpwa::massDepFit::function::chiSquare(const double* par) const
{
	// in C++11 we can use a static variable per thread so that the
	// parameters are kept over function calls and we can implement some
	// caching
	thread_local rpwa::massDepFit::parameters fitParameters(_compset->getNrComponents()+1,           // nr components + final-state mass-dependence
	                                                        _compset->getMaxChannelsInComponent(),
	                                                        _compset->getMaxParametersInComponent(),
	                                                        _nrBins);
	thread_local rpwa::massDepFit::cache cache(_nrWaves,
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
rpwa::massDepFit::function::logPriorLikelihood(const std::vector<double>& par) const
{
	return logPriorLikelihood(par.data());
}


double
rpwa::massDepFit::function::logPriorLikelihood(const double* par) const
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

	return logPriorLikelihood(fitParameters);
}


double
rpwa::massDepFit::function::logPriorLikelihood(const rpwa::massDepFit::parameters& fitParameters) const
{
	double logPrior = 0;

	const size_t nrComponents = _compset->getNrComponents();
	for (size_t idxComponent = 0; idxComponent < nrComponents; ++idxComponent) {
		const rpwa::massDepFit::component* component = _compset->getComponent(idxComponent);

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

					const size_t idx = _nrWaves*(_nrWaves+1) - (_nrWaves-idxWave)*(_nrWaves-idxWave+1) + 2*(jdxWave-idxWave);
					if(idxWave==jdxWave) {
						chi2 += norm(rhoDiff) * _spinDensityMatricesCovMatInv[idxBin][idxMass](idx, idx);
					} else {
						if (_useCovariance == useDiagnalElementsOnly) {
							chi2 += rhoDiff.real()*rhoDiff.real() * _spinDensityMatricesCovMatInv[idxBin][idxMass](idx,   idx  );
							chi2 += rhoDiff.imag()*rhoDiff.imag() * _spinDensityMatricesCovMatInv[idxBin][idxMass](idx+1, idx+1);
						} else if(_useCovariance == useComplexDiagnalElementsOnly) {
							chi2 += rhoDiff.real()*rhoDiff.real() * _spinDensityMatricesCovMatInv[idxBin][idxMass](idx,   idx  );
							chi2 += rhoDiff.real()*rhoDiff.imag() * _spinDensityMatricesCovMatInv[idxBin][idxMass](idx,   idx+1);
							chi2 += rhoDiff.real()*rhoDiff.imag() * _spinDensityMatricesCovMatInv[idxBin][idxMass](idx+1, idx  );
							chi2 += rhoDiff.imag()*rhoDiff.imag() * _spinDensityMatricesCovMatInv[idxBin][idxMass](idx+1, idx+1);
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
