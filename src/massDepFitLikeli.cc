#include "massDepFitLikeli.h"

#include "massDepFitModel.h"
#include "massDepFitParameters.h"
#include "reportingUtils.hpp"


rpwa::massDepFit::likelihood::likelihood(const bool fitProductionAmplitudes,
                                         const rpwa::massDepFit::likelihood::useCovarianceMatrix useCovariance)
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
	output << "created 'likelihood' object for a fit to the ";
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


rpwa::massDepFit::likelihood*
rpwa::massDepFit::likelihood::Clone() const
{
	return new likelihood(*this);
}


unsigned int
rpwa::massDepFit::likelihood::NDim() const
{
	return _compset->getNrParameters();
}


unsigned int
rpwa::massDepFit::likelihood::NDataPoints() const
{
	unsigned int nrPts(0);

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


bool
rpwa::massDepFit::likelihood::init(rpwa::massDepFit::model* compset,
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

	// do some stuff specific to the fit to the production amplitudes
	if(_fitProductionAmplitudes) {
		// test that the anchor wave is non-zero over the complete fit range
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

					if(zeroThisWave) {
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
				// import covariance matrix of production amplitudes
				const size_t matrixSize = 2 * (_nrWaves - zeroWaves[idxBin][idxMass].size());
				TMatrixT<double> reducedCovMat(matrixSize - 1, matrixSize - 1);

				if(realAnchorWave) {
					size_t idxSkip=0;
					for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
						if(idxSkip < zeroWaves[idxBin][idxMass].size() && zeroWaves[idxBin][idxMass][idxSkip] == idxWave) {
							++idxSkip;
							continue;
						}

						size_t jdxSkip=0;
						for(size_t jdxWave=0; jdxWave<_nrWaves; ++jdxWave) {
							if(jdxSkip < zeroWaves[idxBin][idxMass].size() && zeroWaves[idxBin][idxMass][jdxSkip] == jdxWave) {
								++jdxSkip;
								continue;
							}

							const Int_t row = 2*(idxWave-idxSkip) + (idxWave>_idxAnchorWave ? -1 : 0);
							const Int_t col = 2*(jdxWave-jdxSkip) + (jdxWave>_idxAnchorWave ? -1 : 0);

							reducedCovMat(row + 0, col + 0) = _productionAmplitudesCovariance[idxBin][idxMass][idxWave][jdxWave][0][0];
							if(jdxWave != _idxAnchorWave) {
								reducedCovMat(row + 0, col + 1) = _productionAmplitudesCovariance[idxBin][idxMass][idxWave][jdxWave][0][1];
							}
							if(idxWave != _idxAnchorWave) {
								reducedCovMat(row + 1, col + 0) = _productionAmplitudesCovariance[idxBin][idxMass][idxWave][jdxWave][1][0];
							}
							if(idxWave != _idxAnchorWave && jdxWave != _idxAnchorWave) {
								reducedCovMat(row + 1, col + 1) = _productionAmplitudesCovariance[idxBin][idxMass][idxWave][jdxWave][1][1];
							}
						}
					}
				} else {
					TMatrixT<double> covariance(matrixSize, matrixSize);
					TMatrixT<double> jacobian(matrixSize - 1, matrixSize);

					size_t idxSkip=0;
					for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
						if(idxSkip < zeroWaves[idxBin][idxMass].size() && zeroWaves[idxBin][idxMass][idxSkip] == idxWave) {
							++idxSkip;
							continue;
						}

						size_t jdxSkip=0;
						for(size_t jdxWave=0; jdxWave<_nrWaves; ++jdxWave) {
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

							const Int_t row = 2*(idxWave-idxSkip) + (idxWave>_idxAnchorWave ? -1 : 0);
							if(idxWave == _idxAnchorWave && jdxWave == _idxAnchorWave) {
								jacobian(row + 0, 2*(jdxWave-jdxSkip) + 0) = xa1 / n;
								jacobian(row + 0, 2*(jdxWave-jdxSkip) + 1) = xa2 / n;
							} else if(jdxWave == _idxAnchorWave) {
								jacobian(row + 0, 2*(jdxWave-jdxSkip) + 0) = xi1 / n - xa1 * (xi1*xa1 + xi2*xa2) / n3;
								jacobian(row + 0, 2*(jdxWave-jdxSkip) + 1) = xi2 / n - xa2 * (xi1*xa1 + xi2*xa2) / n3;
								if(idxWave != _idxAnchorWave) {
									jacobian(row + 1, 2*(jdxWave-jdxSkip) + 0) =   xi2 / n - xa1 * (xi2*xa1 - xi1*xa2) / n3;
									jacobian(row + 1, 2*(jdxWave-jdxSkip) + 1) = - xi1 / n - xa2 * (xi2*xa1 - xi1*xa2) / n3;
								}
							} else if(idxWave == jdxWave) {
								jacobian(row + 0, 2*(jdxWave-jdxSkip) + 0) =   xa1 / n;
								jacobian(row + 0, 2*(jdxWave-jdxSkip) + 1) =   xa2 / n;
								jacobian(row + 1, 2*(jdxWave-jdxSkip) + 0) = - xa2 / n;
								jacobian(row + 1, 2*(jdxWave-jdxSkip) + 1) =   xa1 / n;
							}
						}
					}

					TMatrixT<double> jacobianT(TMatrixT<double>::kTransposed, jacobian);

					reducedCovMat = jacobian * covariance * jacobianT;
				}

				reducedCovMat.Invert();

				// import covariance matrix of production amplitudes
				_productionAmplitudesCovMatInv[idxBin][idxMass].ResizeTo(2*_nrWaves - 1, 2*_nrWaves - 1);

				size_t idxSkip=0;
				for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
					if(idxSkip < zeroWaves[idxBin][idxMass].size() && zeroWaves[idxBin][idxMass][idxSkip] == idxWave) {
						++idxSkip;
						continue;
					}

					size_t jdxSkip=0;
					for(size_t jdxWave=0; jdxWave<_nrWaves; ++jdxWave) {
						if(jdxSkip < zeroWaves[idxBin][idxMass].size() && zeroWaves[idxBin][idxMass][jdxSkip] == jdxWave) {
							++jdxSkip;
							continue;
						}

						if(idxWave != jdxWave && _useCovariance != useFullCovarianceMatrix) {
							continue;
						}

						const Int_t row = 2*idxWave + (idxWave>_idxAnchorWave ? -1 : 0);
						const Int_t col = 2*jdxWave + (jdxWave>_idxAnchorWave ? -1 : 0);

						_productionAmplitudesCovMatInv[idxBin][idxMass](row + 0, col + 0) = reducedCovMat(row - 2*idxSkip + 0, col - 2*jdxSkip + 0);
						if(jdxWave != _idxAnchorWave && _useCovariance != useDiagnalElementsOnly) {
							_productionAmplitudesCovMatInv[idxBin][idxMass](row + 0, col + 1) = reducedCovMat(row - 2*idxSkip + 0, col - 2*jdxSkip + 1);
						}
						if(idxWave != _idxAnchorWave && _useCovariance != useDiagnalElementsOnly) {
							_productionAmplitudesCovMatInv[idxBin][idxMass](row + 1, col + 0) = reducedCovMat(row - 2*idxSkip + 1, col - 2*jdxSkip + 0);
						}
						if(idxWave != _idxAnchorWave && jdxWave != _idxAnchorWave) {
							_productionAmplitudesCovMatInv[idxBin][idxMass](row + 1, col + 1) = reducedCovMat(row - 2*idxSkip + 1, col - 2*jdxSkip + 1);
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
	}

	return true;
}


double
rpwa::massDepFit::likelihood::DoEval(const double* par) const
{
	// import parameters (couplings, branchings, resonance parameters, ...)
	rpwa::massDepFit::parameters fitParameters(_compset->getNrComponents()+1,           // nr components + final-state mass-dependence
	                                           _compset->getMaxChannelsInComponent(),
	                                           _compset->getMaxParametersInComponent(),
	                                           _nrBins);
	_compset->importParameters(par, fitParameters);

	return DoEval(fitParameters);
}


double
rpwa::massDepFit::likelihood::DoEval(const rpwa::massDepFit::parameters& fitParameters) const
{
	if(_fitProductionAmplitudes) {
		return DoEvalProductionAmplitudes(fitParameters);
	} else {
		return DoEvalSpinDensityMatrix(fitParameters);
	}
}


double
rpwa::massDepFit::likelihood::DoEvalProductionAmplitudes(const rpwa::massDepFit::parameters& fitParameters) const
{
	double chi2=0;

	// loop over bins
	for(unsigned idxBin=0; idxBin<_nrBins; ++idxBin) {
		// loop over mass-bins
		for(unsigned idxMass=_idxMassMin; idxMass<=_idxMassMax; ++idxMass) {
			const double mass = _massBinCenters[idxMass];

			// phase of fit in anchor wave
			const std::complex<double> anchorFit = _compset->productionAmplitude(fitParameters, _idxAnchorWave, idxBin, mass, idxMass);
			const std::complex<double> anchorFitPhase = anchorFit / abs(anchorFit);

			TMatrixT<double> prodAmpDiffMat(2*_nrWaves - 1, 1);

			// sum over the contributions to chi2
			for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
				// check that this mass bin should be taken into account for this
				// combination of waves
				if(idxMass < _wavePairMassBinLimits[idxWave][idxWave].first || idxMass > _wavePairMassBinLimits[idxWave][idxWave].second) {
					continue;
				}

				// calculate target spin density matrix element
				const std::complex<double> prodAmpFit = _compset->productionAmplitude(fitParameters, idxWave, idxBin, mass, idxMass) / anchorFitPhase;

				const std::complex<double> prodAmpDiff = prodAmpFit - _productionAmplitudes[idxBin][idxMass][idxWave];

				const Int_t row = 2*idxWave + (idxWave>_idxAnchorWave ? -1 : 0);
				prodAmpDiffMat(row + 0, 0) = prodAmpDiff.real();
				if(idxWave != _idxAnchorWave) {
					prodAmpDiffMat(row + 1, 0) = prodAmpDiff.imag();
				}
			} // end loop over idxWave

			const TMatrixT<double> prodAmpDiffMatT(TMatrixT<double>::kTransposed, prodAmpDiffMat);

			const TMatrixT<double> dChi2Mat(prodAmpDiffMatT * _productionAmplitudesCovMatInv[idxBin][idxMass] * prodAmpDiffMat);

			chi2 += dChi2Mat(0, 0);
		} // end loop over mass-bins
	} // end loop over bins

	return chi2;
}


double
rpwa::massDepFit::likelihood::DoEvalSpinDensityMatrix(const rpwa::massDepFit::parameters& fitParameters) const
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
					const std::complex<double> rhoFit = _compset->spinDensityMatrix(fitParameters, idxWave, jdxWave, idxBin, mass, idxMass);

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
