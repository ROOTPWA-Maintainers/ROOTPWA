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

#include <TVectorT.h>

#include <reportingUtils.hpp>

#include "cache.h"
#include "components.h"
#include "data.h"
#include "model.h"
#include "parameters.h"


rpwa::resonanceFit::function::function(const rpwa::resonanceFit::dataConstPtr& fitData,
                                       const rpwa::resonanceFit::modelConstPtr& fitModel,
                                       const bool useProductionAmplitudes)
	: _fitData(fitData),
	  _fitModel(fitModel),
	  _nrBins(fitData->nrBins()),
	  _maxNrWaves(_fitData->maxNrWaves()),
	  _maxNrMassBins(_fitData->maxNrMassBins()),
	  _useProductionAmplitudes(useProductionAmplitudes),
	  _useCovariance(_fitData->useCovariance())
{
	if(not _useProductionAmplitudes and _useCovariance == useFullCovarianceMatrix) {
		printErr << "cannot use full covariance matrix while fitting to spin-density matrix." << std::endl;
		throw;
	}

	_idxMassMax.resize(_nrBins);
	_idxMassMin.resize(_nrBins);
	for(size_t idxBin = 0; idxBin < _nrBins; ++idxBin) {
		_idxMassMin[idxBin] = _maxNrMassBins;
		_idxMassMax[idxBin] = 0;
		for(size_t idxWave = 0; idxWave < _fitData->nrWaves(idxBin); ++idxWave) {
			_idxMassMin[idxBin] = std::min(_idxMassMin[idxBin], _fitData->wavePairMassBinLimits()[idxBin][idxWave][idxWave].first);
			_idxMassMax[idxBin] = std::max(_idxMassMax[idxBin], _fitData->wavePairMassBinLimits()[idxBin][idxWave][idxWave].second);
		}
	}

	std::ostringstream output;
	output << "created 'function' object for a fit to the ";
	if(_useProductionAmplitudes) {
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


size_t
rpwa::resonanceFit::function::getNrParameters() const
{
	return _fitModel->getNrParameters();
}


size_t
rpwa::resonanceFit::function::getNrDataPoints() const
{
	size_t nrPts(0);

	if(_useProductionAmplitudes) {
		// calculate data points:
		// * production amplitudes in general are complex numbers
		// * for the anchor wave it might be real
		// * remember (Re,Im) => factor 2
		for(size_t idxBin = 0; idxBin < _nrBins; ++idxBin) {
			for(size_t idxWave = 0; idxWave < _fitData->nrWaves(idxBin); ++idxWave) {
				nrPts += _fitData->wavePairMassBinLimits()[idxBin][idxWave][idxWave].second - _fitData->wavePairMassBinLimits()[idxBin][idxWave][idxWave].first + 1;
				if(idxWave != _fitModel->anchorWaveIndex(idxBin)) {
					nrPts += _fitData->wavePairMassBinLimits()[idxBin][idxWave][idxWave].second - _fitData->wavePairMassBinLimits()[idxBin][idxWave][idxWave].first + 1;
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
			for(size_t idxWave = 0; idxWave < _fitData->nrWaves(idxBin); ++idxWave) {
				for(size_t jdxWave = 0; jdxWave < _fitData->nrWaves(idxBin); ++jdxWave) {
					nrPts += _fitData->wavePairMassBinLimits()[idxBin][idxWave][jdxWave].second - _fitData->wavePairMassBinLimits()[idxBin][idxWave][jdxWave].first + 1;
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
	thread_local rpwa::resonanceFit::cache cache(_maxNrWaves,
	                                             _fitModel->getNrComponents()+1,          // nr components + final-state mass-dependence
	                                             _fitModel->getMaxChannelsInComponent(),
	                                             _nrBins,
	                                             _maxNrMassBins);

	// import parameters (couplings, branchings, resonance parameters, ...)
	_fitModel->importParameters(par, fitParameters, cache);

	return chiSquare(fitParameters, cache);
}


double
rpwa::resonanceFit::function::chiSquare(const rpwa::resonanceFit::parameters& fitParameters,
                                        rpwa::resonanceFit::cache& cache) const
{
	if(_useProductionAmplitudes) {
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
	rpwa::resonanceFit::cache cache(_maxNrWaves,
	                                _fitModel->getNrComponents()+1,          // nr components + final-state mass-dependence
	                                _fitModel->getMaxChannelsInComponent(),
	                                _nrBins,
	                                _maxNrMassBins);

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
			if(component->getParameter(idxParameter).fixed())
				continue;

			// parameters with 0 error are assumed to have a flat prior
			if(component->getParameter(idxParameter).startError() == 0.0)
				continue;

			logPrior += -0.5 * std::pow((fitParameters.getParameter(idxComponent, idxParameter) - component->getParameter(idxParameter).startValue()) / component->getParameter(idxParameter).startError(), 2.);
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
			const double mass = _fitData->massBinCenters()[idxBin][idxMass];

			// phase of fit in anchor wave
			const std::complex<double> anchorFit = _fitModel->productionAmplitude(fitParameters, cache, _fitModel->anchorWaveIndex(idxBin), idxBin, mass, idxMass);
			const std::complex<double> anchorFitPhase = anchorFit / abs(anchorFit);

			TVectorT<double> prodAmpDiffVect(2*_fitData->nrWaves(idxBin));

			// sum over the contributions to chi2
			for(size_t idxWave = 0; idxWave < _fitData->nrWaves(idxBin); ++idxWave) {
				// check that this mass bin should be taken into account for this
				// combination of waves
				if(idxMass < _fitData->wavePairMassBinLimits()[idxBin][idxWave][idxWave].first or idxMass > _fitData->wavePairMassBinLimits()[idxBin][idxWave][idxWave].second) {
					continue;
				}

				// calculate target spin density matrix element
				const std::complex<double> prodAmpFit = _fitModel->productionAmplitude(fitParameters, cache, idxWave, idxBin, mass, idxMass) / anchorFitPhase;

				const std::complex<double> prodAmpDiff = prodAmpFit - _fitData->productionAmplitudes()[idxBin][idxMass][idxWave];

				const Int_t row = 2*idxWave;
				prodAmpDiffVect(row) = prodAmpDiff.real();
				if(idxWave != _fitModel->anchorWaveIndex(idxBin)) {
					prodAmpDiffVect(row+1) = prodAmpDiff.imag();
				}
			} // end loop over idxWave

			chi2 += _fitData->productionAmplitudesCovMatInv()[idxBin][idxMass].Similarity(prodAmpDiffVect);
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
			const double mass = _fitData->massBinCenters()[idxBin][idxMass];

			// sum over the contributions to chi2 -> rho_ij
			for(size_t idxWave = 0; idxWave < _fitData->nrWaves(idxBin); ++idxWave) {
				for(size_t jdxWave = idxWave; jdxWave < _fitData->nrWaves(idxBin); ++jdxWave) {
					// check that this mass bin should be taken into account for this
					// combination of waves
					if(idxMass < _fitData->wavePairMassBinLimits()[idxBin][idxWave][jdxWave].first or idxMass > _fitData->wavePairMassBinLimits()[idxBin][idxWave][jdxWave].second) {
						continue;
					}

					// calculate target spin density matrix element
					const std::complex<double> rhoFit = _fitModel->spinDensityMatrix(fitParameters, cache, idxWave, jdxWave, idxBin, mass, idxMass);

					const std::complex<double> rhoDiff = rhoFit - _fitData->spinDensityMatrixElements()[idxBin][idxMass][idxWave][jdxWave];

					if(idxWave==jdxWave) {
						chi2 += rhoDiff.real() * _fitData->spinDensityMatrixElementsCovMatInvArray()[idxBin][idxMass][idxWave][jdxWave][0][0] * rhoDiff.real();
					} else {
						if (_useCovariance == useDiagnalElementsOnly) {
							chi2 += rhoDiff.real() * _fitData->spinDensityMatrixElementsCovMatInvArray()[idxBin][idxMass][idxWave][jdxWave][0][0] * rhoDiff.real();
							chi2 += rhoDiff.imag() * _fitData->spinDensityMatrixElementsCovMatInvArray()[idxBin][idxMass][idxWave][jdxWave][1][1] * rhoDiff.imag();
						} else if(_useCovariance == useComplexDiagnalElementsOnly) {
							chi2 += rhoDiff.real() * _fitData->spinDensityMatrixElementsCovMatInvArray()[idxBin][idxMass][idxWave][jdxWave][0][0] * rhoDiff.real();
							chi2 += rhoDiff.real() * _fitData->spinDensityMatrixElementsCovMatInvArray()[idxBin][idxMass][idxWave][jdxWave][0][1] * rhoDiff.imag();
							chi2 += rhoDiff.imag() * _fitData->spinDensityMatrixElementsCovMatInvArray()[idxBin][idxMass][idxWave][jdxWave][1][0] * rhoDiff.real();
							chi2 += rhoDiff.imag() * _fitData->spinDensityMatrixElementsCovMatInvArray()[idxBin][idxMass][idxWave][jdxWave][1][1] * rhoDiff.imag();
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
