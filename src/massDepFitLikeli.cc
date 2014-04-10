#include "massDepFitLikeli.h"

#include "massDepFitModel.h"
#include "reportingUtils.hpp"


rpwa::massDepFit::likelihood*
rpwa::massDepFit::likelihood::Clone() const {
	return new likelihood(*this);
}


unsigned int
rpwa::massDepFit::likelihood::NDim() const {
	return _compset->getNrParameters();
}


unsigned int
rpwa::massDepFit::likelihood::NDataPoints() const {
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
                                   const boost::multi_array<std::pair<size_t, size_t>, 2>& wavePairMassBinLimits,
                                   bool fitProductionAmplitudes,
                                   bool useCovariance)
{
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

	_fitProductionAmplitudes = fitProductionAmplitudes;
	_useCovariance = useCovariance;

	_nrBins = _spinDensityMatrices.size();
	_nrMassBins = _massBinCenters.size();
	_nrWaves = _wavePairMassBinLimits.size();

	// do some stuff specific to the fit to the production amplitudes
	if(_fitProductionAmplitudes) {
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

		// at the moment production amplitude can only be fitted, if the anchor wave is real valued
		if(not realAnchorWave) {
			printErr << "production amplitudes cannot be fitted if the anchor wave is not real valued." << std::endl;
			return false;
		}

		// error if any non-anchor wave is real
		if(realOtherWaves) {
			printErr << "production amplitudes cannot be fitted if a non-anchor wave is real valued." << std::endl;
			return false;
		}

		_productionAmplitudesCovMatInv.resize(boost::extents[_nrBins][_nrMassBins]);
		for(size_t idxBin=0; idxBin<_nrBins; ++idxBin) {
			for(size_t idxMass=0; idxMass<_nrMassBins; ++idxMass) {
				// import covariance matrix of production amplitudes
				_productionAmplitudesCovMatInv[idxBin][idxMass].ResizeTo(2*_nrWaves - 1, 2*_nrWaves - 1);

				for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
					for(size_t jdxWave=0; jdxWave<_nrWaves; ++jdxWave) {
						if(idxWave != jdxWave && not _useCovariance) {
							continue;
						}

						const Int_t row = 2*idxWave + (idxWave>_idxAnchorWave ? -1 : 0);
						const Int_t col = 2*jdxWave + (jdxWave>_idxAnchorWave ? -1 : 0);

						_productionAmplitudesCovMatInv[idxBin][idxMass](row + 0, col + 0) = _productionAmplitudesCovariance[idxBin][idxMass][idxWave][jdxWave][0][0];
						if(idxWave != _idxAnchorWave || jdxWave != _idxAnchorWave) {
							_productionAmplitudesCovMatInv[idxBin][idxMass](row + 0, col + 1) = _productionAmplitudesCovariance[idxBin][idxMass][idxWave][jdxWave][0][1];
							_productionAmplitudesCovMatInv[idxBin][idxMass](row + 1, col + 0) = _productionAmplitudesCovariance[idxBin][idxMass][idxWave][jdxWave][1][0];
							_productionAmplitudesCovMatInv[idxBin][idxMass](row + 1, col + 1) = _productionAmplitudesCovariance[idxBin][idxMass][idxWave][jdxWave][1][1];
						}
					}
				}

				_productionAmplitudesCovMatInv[idxBin][idxMass].Invert();
			}
		}

		// modify measured production amplitude such that the anchor wave is always real and positive
		for(size_t idxBin=0; idxBin<_nrBins; ++idxBin) {
			for(size_t idxMass=0; idxMass<_nrMassBins; ++idxMass) {
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
rpwa::massDepFit::likelihood::DoEval(const double* par) const {
	// set parameters for resonances, background and phase space
	_compset->setParameters(par);

	if(_fitProductionAmplitudes) {
		return DoEvalProductionAmplitudes();
	} else {
		return DoEvalSpinDensityMatrix();
	}
}


double
rpwa::massDepFit::likelihood::DoEvalProductionAmplitudes() const {
	double chi2=0;

	// loop over bins
	for(unsigned idxBin=0; idxBin<_nrBins; ++idxBin) {
		// loop over mass-bins
		for(unsigned idxMass=0; idxMass<_nrMassBins; ++idxMass) {
			const double mass = _massBinCenters[idxMass];

			TMatrixT<double> prodAmpDiffMat(2*_nrWaves - 1, 1);

			// sum over the contributions to chi2
			for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
				// check that this mass bin should be taken into account for this
				// combination of waves
				if(idxMass < _wavePairMassBinLimits[idxWave][idxWave].first || idxMass > _wavePairMassBinLimits[idxWave][idxWave].second) {
					continue;
				}

				// phase of fit in anchor wave
				const std::complex<double> anchorFit = _compset->productionAmplitude(_idxAnchorWave, idxBin, mass, idxMass);
				const std::complex<double> anchorFitPhase = anchorFit / abs(anchorFit);

				// calculate target spin density matrix element
				const std::complex<double> prodAmpFit = _compset->productionAmplitude(idxWave, idxBin, mass, idxMass) / anchorFitPhase;

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
rpwa::massDepFit::likelihood::DoEvalSpinDensityMatrix() const {
	double chi2=0;

	// loop over bins
	for(unsigned idxBin=0; idxBin<_nrBins; ++idxBin) {
		// loop over mass-bins
		for(unsigned idxMass=0; idxMass<_nrMassBins; ++idxMass) {
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
					const std::complex<double> rhoFit = _compset->spinDensityMatrix(idxWave, jdxWave, idxBin, mass, idxMass);

					const std::complex<double> rhoDiff = rhoFit - _spinDensityMatrices[idxBin][idxMass][idxWave][jdxWave];

					double dchi;
					if(idxWave==jdxWave) {
						dchi = norm(rhoDiff) / _spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][0][0];
					} else {
						if(_useCovariance) {
							dchi  = rhoDiff.real()*rhoDiff.real() * _spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][1][1];
							dchi -= rhoDiff.real()*rhoDiff.imag() * _spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][0][1];
							dchi -= rhoDiff.real()*rhoDiff.imag() * _spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][1][0];
							dchi += rhoDiff.imag()*rhoDiff.imag() * _spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][0][0];

							dchi /= _spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][0][0]*_spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][1][1]
							        - _spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][0][1]*_spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][1][0];
						} else {
							dchi  = rhoDiff.real()*rhoDiff.real() / _spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][0][0];
							dchi += rhoDiff.imag()*rhoDiff.imag() / _spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][1][1];
						}
					}
					chi2 += dchi;
				} // end loop over jdxWave
			} // end loop over idxWave
		} // end loop over mass-bins
	} // end loop over bins

	return chi2;
}
