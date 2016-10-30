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
//      implementation of the master class of the resonance fit
//
//-------------------------------------------------------------------------


#include "massDepFit.h"

#include <algorithm>
#include <complex>
#include <iostream>
#include <vector>
#include <string>

#include <boost/assign/std/vector.hpp>

#include <yaml-cpp/yaml.h>

#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TTree.h>

#include <fitResult.h>
#include <reportingUtils.hpp>
#include <yamlCppUtils.hpp>

#include "components.h"
#include "data.h"
#include "fsmd.h"
#include "function.h"
#include "information.h"
#include "model.h"
#include "resonanceFitHelper.h"


namespace {


	// debug flag for functions in this (anonymous) namespace
	bool debug = false;


	// do some stuff specific to the fit to the production amplitudes
	// * check that the anchor wave is non-zero over the complete fit range
	// * rotate production amplitudes and covariance matrixes such that the
	//   anchor wave is real
	bool
	prepareProductionAmplitudes(const rpwa::resonanceFit::function::useCovarianceMatrix useCovariance,
	                            const size_t idxAnchorWave,
	                            const boost::multi_array<std::pair<size_t, size_t>, 3>& wavePairMassBinLimits,
	                            boost::multi_array<std::complex<double>, 3>& productionAmplitudes,
	                            boost::multi_array<TMatrixT<double>, 2>& productionAmplitudesCovariance)
	{
		const size_t nrBins = *(productionAmplitudes.shape());
		const size_t maxMassBins = *(productionAmplitudes.shape()+1);
		const size_t nrWaves = *(productionAmplitudes.shape()+2);

		// get range of used mass bins (in any wave) for each bin
		std::vector<size_t> idxMassMax(nrBins);
		std::vector<size_t> idxMassMin(nrBins);
		for(size_t idxBin = 0; idxBin < nrBins; ++idxBin) {
			idxMassMin[idxBin] = maxMassBins;
			idxMassMax[idxBin] = 0;
			for(size_t idxWave = 0; idxWave < nrWaves; ++idxWave) {
				idxMassMin[idxBin] = std::min(idxMassMin[idxBin], wavePairMassBinLimits[idxBin][idxWave][idxWave].first);
				idxMassMax[idxBin] = std::max(idxMassMax[idxBin], wavePairMassBinLimits[idxBin][idxWave][idxWave].second);
			}
		}

		// get a list of waves that are zero (those have to be excluded
		// from the inversion of the covariance matrix below) and test
		// that the anchor wave is non-zero over the complete fit range
		bool zeroAnchorWave = false;
		boost::multi_array<std::vector<size_t>, 2> zeroWaves(boost::extents[nrBins][maxMassBins]);
		for(size_t idxBin = 0; idxBin < nrBins; ++idxBin) {
			for(size_t idxMass = idxMassMin[idxBin]; idxMass <= idxMassMax[idxBin]; ++idxMass) {
				for(size_t idxWave = 0; idxWave < nrWaves; ++idxWave) {
					bool zeroThisWave = true;
					zeroThisWave &= (productionAmplitudes[idxBin][idxMass][idxWave].real() == 0.);
					zeroThisWave &= (productionAmplitudes[idxBin][idxMass][idxWave].imag() == 0.);
					zeroThisWave &= (productionAmplitudesCovariance[idxBin][idxMass](2*idxWave,   2*idxWave  ) == 0.);
					zeroThisWave &= (productionAmplitudesCovariance[idxBin][idxMass](2*idxWave,   2*idxWave+1) == 0.);
					zeroThisWave &= (productionAmplitudesCovariance[idxBin][idxMass](2*idxWave+1, 2*idxWave  ) == 0.);
					zeroThisWave &= (productionAmplitudesCovariance[idxBin][idxMass](2*idxWave+1, 2*idxWave+1) == 0.);

					if(zeroThisWave or idxMass < wavePairMassBinLimits[idxBin][idxWave][idxWave].first or idxMass > wavePairMassBinLimits[idxBin][idxWave][idxWave].second) {
						zeroWaves[idxBin][idxMass].push_back(idxWave);
					}

					if(idxWave == idxAnchorWave) {
						zeroAnchorWave |= zeroThisWave;
					}

					// check that a wave is not zero in its fit range
					if(zeroThisWave and idxMass >= wavePairMassBinLimits[idxBin][idxWave][idxWave].first and idxMass <= wavePairMassBinLimits[idxBin][idxWave][idxWave].second) {
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
		for(size_t idxWave = 0; idxWave < nrWaves; ++idxWave) {
			bool realThisWave = true;
			for(size_t idxBin = 0; idxBin < nrBins; ++idxBin) {
				for(size_t idxMass = idxMassMin[idxBin]; idxMass <= idxMassMax[idxBin]; ++idxMass) {
					realThisWave &= (productionAmplitudes[idxBin][idxMass][idxWave].imag() == 0.);
					realThisWave &= (productionAmplitudesCovariance[idxBin][idxMass](2*idxWave,   2*idxWave+1) == 0.);
					realThisWave &= (productionAmplitudesCovariance[idxBin][idxMass](2*idxWave+1, 2*idxWave  ) == 0.);
					realThisWave &= (productionAmplitudesCovariance[idxBin][idxMass](2*idxWave+1, 2*idxWave+1) == 0.);
				}
			}

			if(idxWave == idxAnchorWave) {
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

		for(size_t idxBin = 0; idxBin < nrBins; ++idxBin) {
			for(size_t idxMass = idxMassMin[idxBin]; idxMass <= idxMassMax[idxBin]; ++idxMass) {
				// determine whether the anchor wave should be used in the current bin
				bool skipAnchor = false;
				for(std::vector<size_t>::const_iterator it = zeroWaves[idxBin][idxMass].begin(); it != zeroWaves[idxBin][idxMass].end(); ++it) {
					if(*it == idxAnchorWave) {
						skipAnchor = true;
					}
				}

				// import covariance matrix of production amplitudes
				const size_t matrixSize = 2 * (nrWaves - zeroWaves[idxBin][idxMass].size()) - (skipAnchor ? 0 : 1);
				TMatrixT<double> reducedCovMat(matrixSize, matrixSize);

				if(realAnchorWave) {
					for(size_t idxWave = 0, idxSkip = 0; idxWave < nrWaves; ++idxWave) {
						if(idxSkip < zeroWaves[idxBin][idxMass].size() and zeroWaves[idxBin][idxMass][idxSkip] == idxWave) {
							++idxSkip;
							continue;
						}

						for(size_t jdxWave = 0, jdxSkip = 0; jdxWave < nrWaves; ++jdxWave) {
							if(jdxSkip < zeroWaves[idxBin][idxMass].size() and zeroWaves[idxBin][idxMass][jdxSkip] == jdxWave) {
								++jdxSkip;
								continue;
							}

							const Int_t rowSkip = 2*(idxWave-idxSkip) + ((idxWave>idxAnchorWave and not skipAnchor) ? -1 : 0);
							const Int_t colSkip = 2*(jdxWave-jdxSkip) + ((jdxWave>idxAnchorWave and not skipAnchor) ? -1 : 0);

							reducedCovMat(rowSkip, colSkip) = productionAmplitudesCovariance[idxBin][idxMass](2*idxWave, 2*jdxWave);
							if(jdxWave != idxAnchorWave) {
								reducedCovMat(rowSkip, colSkip+1) = productionAmplitudesCovariance[idxBin][idxMass](2*idxWave, 2*jdxWave+1);
							}
							if(idxWave != idxAnchorWave) {
								reducedCovMat(rowSkip+1, colSkip) = productionAmplitudesCovariance[idxBin][idxMass](2*idxWave+1, 2*jdxWave);
							}
							if(idxWave != idxAnchorWave and jdxWave != idxAnchorWave) {
								reducedCovMat(rowSkip+1, colSkip+1) = productionAmplitudesCovariance[idxBin][idxMass](2*idxWave+1, 2*jdxWave+1);
							}
						}
					}
				} else {
					TMatrixT<double> covariance(matrixSize + (skipAnchor ? 0 : 1), matrixSize + (skipAnchor ? 0 : 1));
					TMatrixT<double> jacobian(matrixSize, matrixSize + (skipAnchor ? 0 : 1));

					for(size_t idxWave = 0, idxSkip = 0; idxWave < nrWaves; ++idxWave) {
						if(idxSkip < zeroWaves[idxBin][idxMass].size() and zeroWaves[idxBin][idxMass][idxSkip] == idxWave) {
							++idxSkip;
							continue;
						}

						for(size_t jdxWave = 0, jdxSkip = 0; jdxWave < nrWaves; ++jdxWave) {
							if(jdxSkip < zeroWaves[idxBin][idxMass].size() and zeroWaves[idxBin][idxMass][jdxSkip] == jdxWave) {
								++jdxSkip;
								continue;
							}

							covariance(2*(idxWave-idxSkip),   2*(jdxWave-jdxSkip)  ) = productionAmplitudesCovariance[idxBin][idxMass](2*idxWave,   2*jdxWave  );
							covariance(2*(idxWave-idxSkip),   2*(jdxWave-jdxSkip)+1) = productionAmplitudesCovariance[idxBin][idxMass](2*idxWave,   2*jdxWave+1);
							covariance(2*(idxWave-idxSkip)+1, 2*(jdxWave-jdxSkip)  ) = productionAmplitudesCovariance[idxBin][idxMass](2*idxWave+1, 2*jdxWave  );
							covariance(2*(idxWave-idxSkip)+1, 2*(jdxWave-jdxSkip)+1) = productionAmplitudesCovariance[idxBin][idxMass](2*idxWave+1, 2*jdxWave+1);

							const double n = abs(productionAmplitudes[idxBin][idxMass][idxAnchorWave]);
							const double n3 = std::pow(n, 3);
							const double xa1 = productionAmplitudes[idxBin][idxMass][idxAnchorWave].real();
							const double xa2 = productionAmplitudes[idxBin][idxMass][idxAnchorWave].imag();
							const double xi1 = productionAmplitudes[idxBin][idxMass][idxWave].real();
							const double xi2 = productionAmplitudes[idxBin][idxMass][idxWave].imag();

							const Int_t rowSkip = 2*(idxWave-idxSkip) + ((idxWave>idxAnchorWave and not skipAnchor) ? -1 : 0);
							if(idxWave == idxAnchorWave and jdxWave == idxAnchorWave) {
								jacobian(rowSkip, 2*(jdxWave-jdxSkip)  ) = xa1 / n;
								jacobian(rowSkip, 2*(jdxWave-jdxSkip)+1) = xa2 / n;
							} else if(jdxWave == idxAnchorWave) {
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
				if(useCovariance != rpwa::resonanceFit::function::useFullCovarianceMatrix) {
					for(size_t idxWave = 0, idxSkip = 0; idxWave < nrWaves; ++idxWave) {
						if(idxSkip < zeroWaves[idxBin][idxMass].size() and zeroWaves[idxBin][idxMass][idxSkip] == idxWave) {
							++idxSkip;
							continue;
						}

						for(size_t jdxWave = 0, jdxSkip = 0; jdxWave < nrWaves; ++jdxWave) {
							if(jdxSkip < zeroWaves[idxBin][idxMass].size() and zeroWaves[idxBin][idxMass][jdxSkip] == jdxWave) {
								++jdxSkip;
								continue;
							}

							const Int_t rowSkip = 2*(idxWave-idxSkip) + ((idxWave>idxAnchorWave and not skipAnchor) ? -1 : 0);
							const Int_t colSkip = 2*(jdxWave-jdxSkip) + ((jdxWave>idxAnchorWave and not skipAnchor) ? -1 : 0);

							if(idxWave == jdxWave) {
								if(useCovariance == rpwa::resonanceFit::function::useDiagnalElementsOnly) {
									if(jdxWave != idxAnchorWave) {
										reducedCovMat(rowSkip, colSkip+1) = 0;
									}
									if(idxWave != idxAnchorWave) {
										reducedCovMat(rowSkip+1, colSkip) = 0;
									}
								}
							} else {
								reducedCovMat(rowSkip, colSkip) = 0;
								if(jdxWave != idxAnchorWave) {
									reducedCovMat(rowSkip, colSkip+1) = 0;
								}
								if(idxWave != idxAnchorWave) {
									reducedCovMat(rowSkip+1, colSkip) = 0;
								}
								if(idxWave != idxAnchorWave and jdxWave != idxAnchorWave) {
									reducedCovMat(rowSkip+1, colSkip+1) = 0;
								}
							}
						}
					}
				}

				reducedCovMat.Invert();

				// import covariance matrix of production amplitudes
				productionAmplitudesCovariance[idxBin][idxMass].Zero();
				for(size_t idxWave = 0, idxSkip = 0; idxWave < nrWaves; ++idxWave) {
					if(idxSkip < zeroWaves[idxBin][idxMass].size() and zeroWaves[idxBin][idxMass][idxSkip] == idxWave) {
						++idxSkip;
						continue;
					}

					for(size_t jdxWave = 0, jdxSkip = 0; jdxWave < nrWaves; ++jdxWave) {
						if(jdxSkip < zeroWaves[idxBin][idxMass].size() and zeroWaves[idxBin][idxMass][jdxSkip] == jdxWave) {
							++jdxSkip;
							continue;
						}

						const Int_t rowSkip = 2*(idxWave-idxSkip) + ((idxWave>idxAnchorWave and not skipAnchor) ? -1 : 0);
						const Int_t colSkip = 2*(jdxWave-jdxSkip) + ((jdxWave>idxAnchorWave and not skipAnchor) ? -1 : 0);

						productionAmplitudesCovariance[idxBin][idxMass](2*idxWave, 2*jdxWave) = reducedCovMat(rowSkip, colSkip);
						if(jdxWave != idxAnchorWave) {
							productionAmplitudesCovariance[idxBin][idxMass](2*idxWave, 2*jdxWave+1) = reducedCovMat(rowSkip, colSkip+1);
						}
						if(idxWave != idxAnchorWave) {
							productionAmplitudesCovariance[idxBin][idxMass](2*idxWave+1, 2*jdxWave) = reducedCovMat(rowSkip+1, colSkip);
						}
						if(idxWave != idxAnchorWave and jdxWave != idxAnchorWave) {
							productionAmplitudesCovariance[idxBin][idxMass](2*idxWave+1, 2*jdxWave+1) = reducedCovMat(rowSkip+1, colSkip+1);
						}
					}
				}
			}
		}

		// modify measured production amplitude such that the anchor wave is always real and positive
		for(size_t idxBin = 0; idxBin < nrBins; ++idxBin) {
			for(size_t idxMass = idxMassMin[idxBin]; idxMass <= idxMassMax[idxBin]; ++idxMass) {
				const std::complex<double> anchorPhase = productionAmplitudes[idxBin][idxMass][idxAnchorWave] / abs(productionAmplitudes[idxBin][idxMass][idxAnchorWave]);
				for(size_t idxWave = 0; idxWave < nrWaves; ++idxWave) {
					productionAmplitudes[idxBin][idxMass][idxWave] /= anchorPhase;
				}
			}
		}

		return true;
	}


	// do some stuff specific to the fit to the spin-density matrix
	bool
	prepareSpinDensityMatrices(const rpwa::resonanceFit::function::useCovarianceMatrix useCovariance,
	                           const boost::multi_array<std::pair<size_t, size_t>, 3>& wavePairMassBinLimits,
	                           const boost::multi_array<std::complex<double>, 4>& spinDensityMatrices,
	                           boost::multi_array<TMatrixT<double>, 2>& spinDensityMatricesCovariance)
	{
		const size_t nrBins = *(spinDensityMatrices.shape());
		const size_t maxMassBins = *(spinDensityMatrices.shape()+1);
		const size_t nrWaves = *(spinDensityMatrices.shape()+2);

		// get range of used mass bins (in any wave) for each bin
		std::vector<size_t> idxMassMax(nrBins);
		std::vector<size_t> idxMassMin(nrBins);
		for(size_t idxBin = 0; idxBin < nrBins; ++idxBin) {
			idxMassMin[idxBin] = maxMassBins;
			idxMassMax[idxBin] = 0;
			for(size_t idxWave = 0; idxWave < nrWaves; ++idxWave) {
				idxMassMin[idxBin] = std::min(idxMassMin[idxBin], wavePairMassBinLimits[idxBin][idxWave][idxWave].first);
				idxMassMax[idxBin] = std::max(idxMassMax[idxBin], wavePairMassBinLimits[idxBin][idxWave][idxWave].second);
			}
		}

		// get a list of waves that are zero (those have to be excluded
		// from the inversion of the covariance matrix below)
		boost::multi_array<std::vector<size_t>, 2> zeroWaves(boost::extents[nrBins][maxMassBins]);
		for(size_t idxBin = 0; idxBin < nrBins; ++idxBin) {
			for(size_t idxMass = idxMassMin[idxBin]; idxMass <= idxMassMax[idxBin]; ++idxMass) {
				for(size_t idxWave = 0; idxWave < nrWaves; ++idxWave) {
					bool zeroThisWave = true;
					for(size_t jdxWave = 0; jdxWave < nrWaves; ++jdxWave) {
						const size_t idx = nrWaves*(nrWaves+1) - ((jdxWave >= idxWave) ? ((nrWaves-idxWave)*(nrWaves-idxWave+1) - 2*(jdxWave-idxWave)) : ((nrWaves-jdxWave)*(nrWaves-jdxWave+1) - 2*(idxWave-jdxWave)));

						zeroThisWave &= (spinDensityMatrices[idxBin][idxMass][idxWave][jdxWave].real() == 0.);
						zeroThisWave &= (spinDensityMatrices[idxBin][idxMass][idxWave][jdxWave].imag() == 0.);
						zeroThisWave &= (spinDensityMatricesCovariance[idxBin][idxMass](idx,   idx  ) == 0.);
						zeroThisWave &= (spinDensityMatricesCovariance[idxBin][idxMass](idx,   idx+1) == 0.);
						zeroThisWave &= (spinDensityMatricesCovariance[idxBin][idxMass](idx+1, idx  ) == 0.);
						zeroThisWave &= (spinDensityMatricesCovariance[idxBin][idxMass](idx+1, idx+1) == 0.);
					}

					if(zeroThisWave or idxMass < wavePairMassBinLimits[idxBin][idxWave][idxWave].first or idxMass > wavePairMassBinLimits[idxBin][idxWave][idxWave].second) {
						zeroWaves[idxBin][idxMass].push_back(idxWave);
					}

					// check that a wave is not zero in its fit range
					if(zeroThisWave and idxMass >= wavePairMassBinLimits[idxBin][idxWave][idxWave].first and idxMass <= wavePairMassBinLimits[idxBin][idxWave][idxWave].second) {
						printErr << "spin-density matrix element of wave " << idxWave << " zero in its fit range (e.g. mass limit in mass-independent fit)." << std::endl;
						return false;
					}
				}
			}
		}

		for(size_t idxBin = 0; idxBin < nrBins; ++idxBin) {
			for(size_t idxMass = idxMassMin[idxBin]; idxMass <= idxMassMax[idxBin]; ++idxMass) {
				// import covariance matrix of spin-density matrix elements
				const size_t reducedMatrixSize((nrWaves - zeroWaves[idxBin][idxMass].size())*(nrWaves - zeroWaves[idxBin][idxMass].size()));
				TMatrixT<double> reducedCovMat(reducedMatrixSize, reducedMatrixSize);

				{
					// i is for loop over rows
					size_t redIdx = 0;
					for(size_t iWave1 = 0, iSkip1 = 0; iWave1 < nrWaves; ++iWave1) {
						if(iSkip1 < zeroWaves[idxBin][idxMass].size() and zeroWaves[idxBin][idxMass][iSkip1] == iWave1) {
							++iSkip1;
							continue;
						}
						for(size_t iWave2 = iWave1, iSkip2 = iSkip1; iWave2 < nrWaves; ++iWave2) {
							if(iSkip2 < zeroWaves[idxBin][idxMass].size() and zeroWaves[idxBin][idxMass][iSkip2] == iWave2) {
								++iSkip2;
								continue;
							}
							const size_t idx = nrWaves*(nrWaves+1) - (nrWaves-iWave1)*(nrWaves-iWave1+1) + 2*(iWave2-iWave1);

							// j is for loop over columns
							size_t redJdx = 0;
							for(size_t jWave1 = 0, jSkip1 = 0; jWave1 < nrWaves; ++jWave1) {
								if(jSkip1 < zeroWaves[idxBin][idxMass].size() and zeroWaves[idxBin][idxMass][jSkip1] == jWave1) {
									++jSkip1;
									continue;
								}
								for(size_t jWave2 = jWave1, jSkip2 = jSkip1; jWave2 < nrWaves; ++jWave2) {
									if(jSkip2 < zeroWaves[idxBin][idxMass].size() and zeroWaves[idxBin][idxMass][jSkip2] == jWave2) {
										++jSkip2;
										continue;
									}
									const size_t jdx = nrWaves*(nrWaves+1) - (nrWaves-jWave1)*(nrWaves-jWave1+1) + 2*(jWave2-jWave1);

									if(iWave1 == iWave2) { // one row
										if(jWave1 == jWave2) { // one column
											if((iWave1 == jWave1 and iWave2 == jWave2) or useCovariance == rpwa::resonanceFit::function::useFullCovarianceMatrix) {
												reducedCovMat(redIdx,   redJdx  ) = spinDensityMatricesCovariance[idxBin][idxMass](idx,   jdx  );
											}
										} else { // two columns
											if(useCovariance == rpwa::resonanceFit::function::useFullCovarianceMatrix) {
												reducedCovMat(redIdx,   redJdx  ) = spinDensityMatricesCovariance[idxBin][idxMass](idx,   jdx  );
												reducedCovMat(redIdx,   redJdx+1) = spinDensityMatricesCovariance[idxBin][idxMass](idx,   jdx+1);
											}
										}
									} else { // two rows
										if(jWave1 == jWave2) { // one column
											if(useCovariance == rpwa::resonanceFit::function::useFullCovarianceMatrix) {
												reducedCovMat(redIdx,   redJdx  ) = spinDensityMatricesCovariance[idxBin][idxMass](idx,   jdx  );
												reducedCovMat(redIdx+1, redJdx  ) = spinDensityMatricesCovariance[idxBin][idxMass](idx+1, jdx  );
											}
										} else { // two columns
											if((iWave1 == jWave1 and iWave2 == jWave2) or useCovariance == rpwa::resonanceFit::function::useFullCovarianceMatrix) {
												reducedCovMat(redIdx,   redJdx  ) = spinDensityMatricesCovariance[idxBin][idxMass](idx,   jdx  );
												reducedCovMat(redIdx+1, redJdx+1) = spinDensityMatricesCovariance[idxBin][idxMass](idx+1, jdx+1);
												if(useCovariance != rpwa::resonanceFit::function::useDiagnalElementsOnly) {
													reducedCovMat(redIdx,   redJdx+1) = spinDensityMatricesCovariance[idxBin][idxMass](idx,   jdx+1);
													reducedCovMat(redIdx+1, redJdx  ) = spinDensityMatricesCovariance[idxBin][idxMass](idx+1, jdx  );
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
				spinDensityMatricesCovariance[idxBin][idxMass].Zero();

				// i is for loop over rows
				size_t redIdx = 0;
				for(size_t iWave1 = 0, iSkip1 = 0; iWave1 < nrWaves; ++iWave1) {
					if(iSkip1 < zeroWaves[idxBin][idxMass].size() and zeroWaves[idxBin][idxMass][iSkip1] == iWave1) {
						++iSkip1;
						continue;
					}
					for(size_t iWave2 = iWave1, iSkip2 = iSkip1; iWave2 < nrWaves; ++iWave2) {
						if(iSkip2 < zeroWaves[idxBin][idxMass].size() and zeroWaves[idxBin][idxMass][iSkip2] == iWave2) {
							++iSkip2;
							continue;
						}
						const size_t idx = nrWaves*(nrWaves+1) - (nrWaves-iWave1)*(nrWaves-iWave1+1) + 2*(iWave2-iWave1);

						// j is for loop over columns
						size_t redJdx = 0;
						for(size_t jWave1 = 0, jSkip1 = 0; jWave1 < nrWaves; ++jWave1) {
							if(jSkip1 < zeroWaves[idxBin][idxMass].size() and zeroWaves[idxBin][idxMass][jSkip1] == jWave1) {
								++jSkip1;
								continue;
							}
							for(size_t jWave2 = jWave1, jSkip2 = jSkip1; jWave2 < nrWaves; ++jWave2) {
								if(jSkip2 < zeroWaves[idxBin][idxMass].size() and zeroWaves[idxBin][idxMass][jSkip2] == jWave2) {
									++jSkip2;
									continue;
								}
								const size_t jdx = nrWaves*(nrWaves+1) - (nrWaves-jWave1)*(nrWaves-jWave1+1) + 2*(jWave2-jWave1);

								if(iWave1 == iWave2) { // one row
									if(jWave1 == jWave2) { // one column
										spinDensityMatricesCovariance[idxBin][idxMass](idx,   jdx  ) = reducedCovMat(redIdx,   redJdx  );
									} else { // two columns
										spinDensityMatricesCovariance[idxBin][idxMass](idx,   jdx  ) = reducedCovMat(redIdx,   redJdx  );
										spinDensityMatricesCovariance[idxBin][idxMass](idx,   jdx+1) = reducedCovMat(redIdx,   redJdx+1);
									}
								} else { // two rows
									if(jWave1 == jWave2) { // one column
										spinDensityMatricesCovariance[idxBin][idxMass](idx,   jdx  ) = reducedCovMat(redIdx,   redJdx  );
										spinDensityMatricesCovariance[idxBin][idxMass](idx+1, jdx  ) = reducedCovMat(redIdx+1, redJdx  );
									} else { // two columns
										spinDensityMatricesCovariance[idxBin][idxMass](idx,   jdx  ) = reducedCovMat(redIdx,   redJdx  );
										spinDensityMatricesCovariance[idxBin][idxMass](idx,   jdx+1) = reducedCovMat(redIdx,   redJdx+1);
										spinDensityMatricesCovariance[idxBin][idxMass](idx+1, jdx  ) = reducedCovMat(redIdx+1, redJdx  );
										spinDensityMatricesCovariance[idxBin][idxMass](idx+1, jdx+1) = reducedCovMat(redIdx+1, redJdx+1);
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

		return true;
	}


	std::map<std::string, double>
	readFitQuality(const YAML::Node& configRoot)
	{
		if(debug) {
			printDebug << "reading 'fitquality'." << std::endl;
		}

		const YAML::Node& configFitQuality = configRoot["fitquality"];

		std::map<std::string, double> fitQuality;
		if(not configFitQuality) {
			// it is perfectly okay for a config file to not
			// contain fit quality information
			return fitQuality;
		}

		if(not configFitQuality.IsMap()) {
			printErr << "'fitquality' is not a YAML map." << std::endl;
			throw;
		}

		for(YAML::const_iterator it = configFitQuality.begin(); it != configFitQuality.end(); ++it) {
			if(not checkVariableType(it->first, rpwa::YamlCppUtils::TypeString) or not checkVariableType(it->second, rpwa::YamlCppUtils::TypeFloat)) {
				printErr << "entries in 'fitquality' must be pairs of 'string' and 'double'." << std::endl;
				throw;
			}

			const std::string key = it->first.as<std::string>();
			const double value = it->second.as<double>();

			if(fitQuality.count(key) != 0) {
				printErr << "variable '" << key << "' of 'fitquality' given multiple times." << std::endl;
				throw;
			}
			fitQuality[key] = value;

			if(debug) {
				printDebug << "read key '" << key << "' with value '" << value << "'." << std::endl;
			}
		}

		return fitQuality;
	}


	void
	writeFitQuality(YAML::Emitter& yamlOutput,
	                const std::map<std::string, double>& fitQuality)
	{
		if(debug) {
			printDebug << "writing 'fitquality'." << std::endl;
		}

		yamlOutput << YAML::Key << "fitquality";
		yamlOutput << YAML::Value;

		yamlOutput << YAML::BeginMap;
		for(std::map<std::string, double>::const_iterator it = fitQuality.begin(); it != fitQuality.end(); ++it) {
			yamlOutput << YAML::Key << it->first;
			yamlOutput << YAML::Value << it->second;
		}
		yamlOutput << YAML::EndMap;
	}


	std::vector<std::string>
	readFreeParameters(const YAML::Node& configRoot)
	{
		if(debug) {
			printDebug << "reading 'freeparameters'." << std::endl;
		}

		const YAML::Node& configFreeParameters = configRoot["freeparameters"];

		std::vector<std::string> freeParameters;
		if(not configFreeParameters) {
			printWarn << "release order of parameters not specified in configuration file, using default one." << std::endl;

			freeParameters.push_back("coupling branching");
			freeParameters.push_back("coupling branching mass m0");
			freeParameters.push_back("*");

			return freeParameters;
		}

		if(not configFreeParameters.IsSequence()) {
			printErr << "'freeparameters' is not a YAML sequence." << std::endl;
			throw;
		}

		const size_t nrItems = configFreeParameters.size();
		if(nrItems == 0) {
			printErr << "'freeparameters' is an empty sequence, when defined it must at least contain one entry." << std::endl;
			throw;
		}

		for(size_t idxItem = 0; idxItem < nrItems; ++idxItem) {
			if(debug) {
				printDebug << "reading of entry " << idxItem << " in 'freeparameters'." << std::endl;
			}

			if(not checkVariableType(configFreeParameters[idxItem], rpwa::YamlCppUtils::TypeString)) {
				printErr << "'freeparameters' entry at index " << idxItem << " is not a string." << std::endl;
				throw;
			}

			freeParameters.push_back(configFreeParameters[idxItem].as<std::string>());

			if(debug) {
				printDebug << "read parameters to release: '" << freeParameters.back() << "'." << std::endl;
			}
		}

		return freeParameters;
	}


	void
	writeFreeParameters(YAML::Emitter& yamlOutput,
	                    const std::vector<std::string>& freeParameters)
	{
		if(debug) {
			printDebug << "writing 'freeparameters'." << std::endl;
		}

		yamlOutput << YAML::Key << "freeparameters";
		yamlOutput << YAML::Value << freeParameters;
	}


	std::vector<rpwa::resonanceFit::information::bin>
	readInformationFitResults(const YAML::Node& configInput)
	{
		if(debug) {
			printDebug << "reading 'fitresults'." << std::endl;
		}

		const YAML::Node& configInputFitResults = configInput["fitresults"];

		if(not configInputFitResults) {
			printErr << "'fitresults' does not exist in 'input'." << std::endl;
			throw;
		}
		if(not configInputFitResults.IsSequence()) {
			printErr << "'fitresults' is not a YAML sequence." << std::endl;
			throw;
		}

		std::vector<rpwa::resonanceFit::information::bin> bins;

		const size_t nrFitResults = configInputFitResults.size();
		for(size_t idxFitResult = 0; idxFitResult < nrFitResults; ++idxFitResult) {
			if(debug) {
				printDebug << "reading of entry " << idxFitResult << " in 'fitresults'." << std::endl;
			}

			const YAML::Node& configInputFitResult = configInputFitResults[idxFitResult];

			std::map<std::string, rpwa::YamlCppUtils::Type> mandatoryArguments;
			boost::assign::insert(mandatoryArguments)
			                     ("name", rpwa::YamlCppUtils::TypeString)
			                     ("tPrimeMean", rpwa::YamlCppUtils::TypeFloat);
			if(not checkIfAllVariablesAreThere(configInputFitResult, mandatoryArguments)) {
				printErr << "'fitresults' entry at index " << idxFitResult << " does not contain all required variables." << std::endl;
				throw;
			}

			const std::string fileName = configInputFitResult["name"].as<std::string>();
			const double tPrimeMean = configInputFitResult["tPrimeMean"].as<double>();

			double rescaleErrors = 1.0;
			if(configInputFitResult["rescaleErrors"]) {
				if(checkVariableType(configInputFitResult["rescaleErrors"], rpwa::YamlCppUtils::TypeFloat)) {
					rescaleErrors = configInputFitResult["rescaleErrors"].as<double>();
				} else {
					printErr << "variable 'rescaleErrors' of 'fitresults' entry at index " << idxFitResult << " is not a floating point number." << std::endl;
					throw;
				}
			}

			if(debug) {
				printDebug << "read file name of fit results of mass-independent fit: '" << fileName << "'." << std::endl;
				printDebug << "read mean t' value: '" << tPrimeMean << "'." << std::endl;
				printDebug << "rescale errors by factor: '" << rescaleErrors << "'." << std::endl;
			}

			std::vector<std::string> sysFileNames;
			// get information for plotting of systematic error
			const YAML::Node& configInputFitResultSystematics = configInputFitResult["systematics"];
			if(configInputFitResultSystematics) {
				if(not configInputFitResultSystematics.IsSequence()) {
					printErr << "'systematics' is not a YAML sequence." << std::endl;
					throw;
				}

				const size_t nrSystematics = configInputFitResultSystematics.size();
				if(debug) {
					printDebug << "going to read information for " << nrSystematics << " files containing information for systematic errors." << std::endl;
				}

				for(size_t idxSystematics = 0; idxSystematics < nrSystematics; ++idxSystematics) {
					if(not checkVariableType(configInputFitResultSystematics[idxSystematics], rpwa::YamlCppUtils::TypeString)) {
						printErr << "'systematics' entry at index " << idxSystematics << " is not a string." << std::endl;
						throw;
					}

					sysFileNames.push_back(configInputFitResultSystematics[idxSystematics].as<std::string>());
				}
			}

			bins.push_back(rpwa::resonanceFit::information::bin(fileName,
			                                                    tPrimeMean,
			                                                    rescaleErrors,
			                                                    sysFileNames));

			if(debug) {
				printDebug << bins.back() << std::endl;
			}
		}

		return bins;
	}


	void
	writeInformationFitResults(YAML::Emitter& yamlOutput,
	                           const std::vector<rpwa::resonanceFit::information::bin>& bins)
	{
		if(debug) {
			printDebug << "writing 'fitresults'." << std::endl;
		}

		yamlOutput << YAML::Key << "fitresults";
		yamlOutput << YAML::Value;

		yamlOutput << YAML::BeginSeq;
		for(std::vector<rpwa::resonanceFit::information::bin>::const_iterator bin = bins.begin(); bin != bins.end(); ++bin) {
			yamlOutput << YAML::BeginMap;

			yamlOutput << YAML::Key << "name";
			yamlOutput << YAML::Value << bin->fileName();

			yamlOutput << YAML::Key << "tPrimeMean";
			yamlOutput << YAML::Value << bin->tPrimeMean();

			if(bin->rescaleErrors() != 1.) {
				yamlOutput << YAML::Key << "rescaleErrors";
				yamlOutput << YAML::Value << bin->rescaleErrors();
			}

			if(bin->sysFileNames().size() > 0) {
				yamlOutput << YAML::Key << "systematics";
				yamlOutput << YAML::Value << bin->sysFileNames();
			}

			yamlOutput << YAML::EndMap;
		}
		yamlOutput << YAML::EndSeq;
	}


	std::vector<rpwa::resonanceFit::information::wave>
	readInformationWaves(const YAML::Node& configInput)
	{
		if(debug) {
			printDebug << "reading 'waves'." << std::endl;
		}

		const YAML::Node& configInputWaves = configInput["waves"];

		if(not configInputWaves) {
			printErr << "'waves' does not exist in 'input'." << std::endl;
			throw;
		}
		if(not configInputWaves.IsSequence()) {
			printErr << "'waves' is not a YAML sequence." << std::endl;
			throw;
		}

		std::vector<rpwa::resonanceFit::information::wave> waves;

		const size_t nrWaves = configInputWaves.size();
		for(size_t idxWave = 0; idxWave < nrWaves; ++idxWave) {
			if(debug) {
				printDebug << "reading of entry " << idxWave << " in 'waves'." << std::endl;
			}

			const YAML::Node& configInputWave = configInputWaves[idxWave];

			std::map<std::string, rpwa::YamlCppUtils::Type> mandatoryArguments;
			boost::assign::insert(mandatoryArguments)
			                     ("name", rpwa::YamlCppUtils::TypeString);
			if(not checkIfAllVariablesAreThere(configInputWave, mandatoryArguments)) {
				printErr << "'waves' entry at index " << idxWave << " does not contain all required variables." << std::endl;
				throw;
			}

			const std::string waveName = configInputWave["name"].as<std::string>();

			// check that wave does not yet exist
			for(size_t idxCheckWave = 0; idxCheckWave < waves.size(); ++idxCheckWave) {
				if(waves[idxCheckWave].waveName() == waveName) {
					printErr << "wave '" << waveName << "' defined twice." << std::endl;
					throw;
				}
				if(std::find(waves[idxCheckWave].waveNameAlternatives().begin(), waves[idxCheckWave].waveNameAlternatives().end(), waveName) != waves[idxCheckWave].waveNameAlternatives().end()) {
					printErr << "wave '" << waveName << "' already defined as alternative name of wave '" << waves[idxCheckWave].waveName() << "'." << std::endl;
					throw;
				}
			}

			double massLower = -1.;
			if(configInputWave["massLower"]) {
				if(checkVariableType(configInputWave["massLower"], rpwa::YamlCppUtils::TypeFloat)) {
					massLower = configInputWave["massLower"].as<double>();
				} else {
					printErr << "variable 'massLower' of 'waves' entry at index " << idxWave << " is not a floating point number." << std::endl;
					throw;
				}
			}
			double massUpper = -1.;
			if(configInputWave["massUpper"]) {
				if(checkVariableType(configInputWave["massUpper"], rpwa::YamlCppUtils::TypeFloat)) {
					massUpper = configInputWave["massUpper"].as<double>();
				} else {
					printErr << "variable 'massUpper' of 'waves' entry at index " << idxWave << " is not a floating point number." << std::endl;
					throw;
				}
			}

			if(debug) {
				printDebug << "read wave name: '" << waveName << "'." << std::endl;
				printDebug << "read mass range: '" << massLower << "' to '" << massUpper << "'." << std::endl;
			}

			std::vector<std::string> waveNameAlternatives;
			if(configInputWave["alternativeNames"]) {
				if(checkVariableType(configInputWave["alternativeNames"], rpwa::YamlCppUtils::TypeSequence)) {
					for(size_t idxAlt = 0; idxAlt < configInputWave["alternativeNames"].size(); ++idxAlt) {
						if(not checkVariableType(configInputWave["alternativeNames"][idxAlt], rpwa::YamlCppUtils::TypeString)) {
							printErr << "element " << idxAlt << " of variable 'alternativeNames' of 'waves' entry at index " << idxWave << " is not a string." << std::endl;
							throw;
						}
						const std::string waveNameAlternative = configInputWave["alternativeNames"][idxAlt].as<std::string>();

						// check that the alternative name does not yet exist
						if(waveNameAlternative == waveName) {
							printErr << "alternative name '" << waveNameAlternative << "' is equal to name of wave '" << waveName << "'." << std::endl;
							throw;
						}
						if(std::find(waveNameAlternatives.begin(), waveNameAlternatives.end(), waveNameAlternative) != waveNameAlternatives.end()) {
							printErr << "alternative name '" << waveNameAlternative << "' of wave '" << waveName << "' defined twice." << std::endl;
							throw;
						}
						for(size_t idxCheckWave = 0; idxCheckWave < waves.size(); ++idxCheckWave) {
							if(waves[idxCheckWave].waveName() == waveNameAlternative) {
								printErr << "alternative name '" << waveNameAlternative << "' of wave '" << waveName << "' already defined as separate wave." << std::endl;
								throw;
							}
							if(std::find(waves[idxCheckWave].waveNameAlternatives().begin(), waves[idxCheckWave].waveNameAlternatives().end(), waveNameAlternative) != waves[idxCheckWave].waveNameAlternatives().end()) {
								printErr << "alternative name '" << waveNameAlternative << "' of wave '" << waveName << "' already defined as alternative name of wave '" << waves[idxCheckWave].waveName() << "'." << std::endl;
								throw;
							}
						}

						waveNameAlternatives.push_back(waveNameAlternative);
					}
				} else {
					printErr << "variable 'alternativeNames' of 'waves' entry at index " << idxWave << " is not a sequence." << std::endl;
					throw;
				}
			}

			waves.push_back(rpwa::resonanceFit::information::wave(waveName,
			                                                      std::make_pair(massLower, massUpper),
			                                                      waveNameAlternatives));

			if(debug) {
				printDebug << waves.back() << std::endl;
			}
		}

		std::ostringstream output;
		for(size_t idxWave = 0; idxWave < nrWaves; ++idxWave) {
			output << "    " << waves[idxWave].waveName() << std::endl;
		}
		printInfo << nrWaves << " waves to be used in fit:" << std::endl
		          << output.str();

		return waves;
	}


	void
	writeInformationWaves(YAML::Emitter& yamlOutput,
	                      const std::vector<rpwa::resonanceFit::information::wave>& waves)
	{
		if(debug) {
			printDebug << "writing 'waves'." << std::endl;
		}

		yamlOutput << YAML::Key << "waves";
		yamlOutput << YAML::Value;

		yamlOutput << YAML::BeginSeq;
		for(std::vector<rpwa::resonanceFit::information::wave>::const_iterator wave = waves.begin(); wave != waves.end(); ++wave) {
			yamlOutput << YAML::BeginMap;

			yamlOutput << YAML::Key << "name";
			yamlOutput << YAML::Value << wave->waveName();

			if(wave->massLimits().first >= 0) {
				yamlOutput << YAML::Key << "massLower";
				yamlOutput << YAML::Value << wave->massLimits().first;
			}
			if(wave->massLimits().second >= 0) {
				yamlOutput << YAML::Key << "massUpper";
				yamlOutput << YAML::Value << wave->massLimits().second;
			}

			if(wave->waveNameAlternatives().size() > 0) {
				yamlOutput << YAML::Key << "alternativeNames";
				yamlOutput << YAML::Value;

				yamlOutput << YAML::Flow;
				yamlOutput << wave->waveNameAlternatives();
				yamlOutput << YAML::Block;
			}

			yamlOutput << YAML::EndMap;
		}
		yamlOutput << YAML::EndSeq;
	}


	rpwa::resonanceFit::informationConstPtr
	readInformation(const YAML::Node& configRoot)
	{
		if(debug) {
			printDebug << "reading 'input'." << std::endl;
		}

		const YAML::Node& configInput = configRoot["input"];

		if(not configInput) {
			printErr << "'input' does not exist in configuration file." << std::endl;
			throw;
		}

		const std::vector<rpwa::resonanceFit::information::bin> bins = readInformationFitResults(configInput);
		const std::vector<rpwa::resonanceFit::information::wave> waves = readInformationWaves(configInput);

		return std::make_shared<rpwa::resonanceFit::information>(bins, waves);
	}


	void
	writeInformation(YAML::Emitter& yamlOutput,
	                 const rpwa::resonanceFit::informationConstPtr& fitInformation)
	{
		if(debug) {
			printDebug << "writing 'input'." << std::endl;
		}

		yamlOutput << YAML::Key << "input";
		yamlOutput << YAML::Value;

		yamlOutput << YAML::BeginMap;

		writeInformationFitResults(yamlOutput, fitInformation->bins());
		writeInformationWaves(yamlOutput, fitInformation->waves());

		yamlOutput << YAML::EndMap;
	}


	void
	readFitResultMassBins(TTree* tree,
	                      rpwa::fitResult* fit,
	                      size_t& nrMassBins,
	                      boost::multi_array<double, 1>& massBinCenters)
	{
		if(not tree or not fit) {
			printErr << "'tree' or 'fit' is not a pointer to a valid object." << std::endl;
			throw;
		}

		// extract data from tree
		const Long64_t nrEntries = tree->GetEntries();

		if(debug) {
			printDebug << "getting center of mass bins from " << nrEntries << " entries in tree." << std::endl;
		}

		nrMassBins = 0;
		for(Long64_t idx = 0; idx < nrEntries; ++idx) {
			if(tree->GetEntry(idx) == 0) {
				printErr << "error while reading entry " << idx << " from tree." << std::endl;
				throw;
			}
			const double newMass = fit->massBinCenter();

			if(debug) {
				printDebug << "entry " << idx << ": center of mass bin at " << newMass << " GeV/c^2" << std::endl;
			}

			bool found = false;
			for(size_t idxMass = 0; idxMass < nrMassBins; ++idxMass) {
				if(std::abs(massBinCenters[idxMass]-newMass) < 1000.*std::numeric_limits<double>::epsilon()) {
					found = true;
					if(debug) {
						printDebug << "this center of mass bin already was encountered before." << std::endl;
					}
					break;
				}
			}

			if(not found) {
				rpwa::resonanceFit::adjustSizeAndSet(massBinCenters, nrMassBins++, newMass);
			}
		} // end loop over entries in tree

		// sort mass bins
		std::sort(massBinCenters.data(), massBinCenters.data() + nrMassBins);

		printInfo << "found " << nrMassBins << " mass bins, center of first and last mass bins: "
		          << massBinCenters[0] << " and " << massBinCenters[nrMassBins - 1] << " GeV/c^2." << std::endl;

		const double massStep = (massBinCenters[nrMassBins-1] - massBinCenters[0]) / (nrMassBins - 1);
		for(size_t idxMass = 1; idxMass < nrMassBins; ++idxMass) {
			if(std::abs(massBinCenters[idxMass]-massBinCenters[idxMass-1] - massStep) > 1000.*std::numeric_limits<double>::epsilon()) {
				printErr << "mass distance between bins " << idxMass-1 << " (" << massBinCenters[idxMass-1] << " GeV/c^2) and "
				         << idxMass << " (" << massBinCenters[idxMass] << " GeV/c^2) does not agree with nominal distance "
				         << massStep << " GeV/c^2" << std::endl;
				throw;
			}
		}
	}


	void
	checkFitResultMassBins(TTree* tree,
	                       rpwa::fitResult* fit,
	                       const size_t nrMassBins,
	                       const boost::multi_array<double, 1>& massBinCenters,
	                       std::vector<Long64_t>& mapping)
	{
		if(not tree or not fit) {
			printErr << "'tree' or 'fit' is not a pointer to a valid object." << std::endl;
			throw;
		}

		// reset mapping
		mapping.assign(nrMassBins, std::numeric_limits<Long64_t>::max());

		// extract data from tree
		const Long64_t nrEntries = tree->GetEntries();

		if(debug) {
			printDebug << "check that the centers of mass bins of " << nrEntries << " entries in tree are at a known place, "
			           << "and map the " << nrMassBins << " mass bins to those entries." << std::endl;
		}

		for(Long64_t idx = 0; idx < nrEntries; ++idx) {
			if(tree->GetEntry(idx) == 0) {
				printErr << "error while reading entry " << idx << " from tree." << std::endl;
				throw;
			}
			//FIXME: this would also be the place to select the best fit in case one file contains more than one fit result per mass bin
			const double mass = fit->massBinCenter();

			if(debug) {
				printDebug << "entry " << idx << ": center of mass bin at " << mass << " GeV/c^2" << std::endl;
			}

			bool found = false;
			size_t idxMass = 0;
			while(idxMass < nrMassBins) {
				if(std::abs(massBinCenters[idxMass]-mass) < 1000.*std::numeric_limits<double>::epsilon()) {
					found = true;
					break;
				}
				++idxMass;
			}

			if(not found) {
				printErr << "could not map mass bin centered at " << mass << " GeV/c^2 to a known mass bin." << std::endl;
				throw;
			}

			if(mapping[idxMass] != std::numeric_limits<Long64_t>::max()) {
				printErr << "cannot map tree entry " << idx << " to mass bin " << idxMass << " (" << massBinCenters[idxMass] << " GeV/c^2)  "
				         << "which is already mapped to tree entry " << mapping[idxMass] << "." << std::endl;
				throw;
			}

			if(debug) {
				printDebug << "mapping mass bin " << idxMass << " (" << massBinCenters[idxMass] << " GeV/c^2) to tree entry " << idx << "." << std::endl;
			}
			mapping[idxMass] = idx;
		} // end loop over entries in tree

		// check that all mass bins are mapped
		for(size_t idx = 0; idx < mapping.size(); ++idx) {
			if(mapping[idx] == std::numeric_limits<Long64_t>::max()) {
				printErr << "mass bin " << idx << " (" << massBinCenters[idx] << " GeV/c^2) not mapped." << std::endl;
				throw;
			}
		}

		if(debug) {
			std::ostringstream output;
			for(size_t idx = 0; idx < mapping.size(); ++idx) {
				output << " " << idx << "->" << mapping[idx];
			}
			printDebug << "established mapping:" << output.str() << std::endl;
		}
	}


	void
	readFitResultWaveNames(const rpwa::resonanceFit::informationConstPtr& fitInformation,
	                       TTree* tree,
	                       rpwa::fitResult* fit,
	                       boost::multi_array<std::string, 1>& waveNames)
	{
		if(not tree or not fit) {
			printErr << "'tree' or 'fit' is not a pointer to a valid object." << std::endl;
			throw;
		}

		if(debug) {
			printDebug << "getting wave names to use for current fit result." << std::endl;
		}

		// read wave names from first fit result in tree
		waveNames.resize(boost::extents[fitInformation->nrWaves()]);
		if(tree->GetEntry(0) == 0) {
			printErr << "error while reading entry " << 0 << " from tree." << std::endl;
			throw;
		}
		for(size_t idxWave = 0; idxWave < fitInformation->nrWaves(); ++idxWave) {
			int idx = fit->waveIndex(fitInformation->getWave(idxWave).waveName());
			// try alternative wave names
			for(size_t idxAlt = 0; idxAlt < fitInformation->getWave(idxWave).waveNameAlternatives().size(); ++idxAlt) {
				const int altIdx = fit->waveIndex(fitInformation->getWave(idxWave).waveNameAlternatives()[idxAlt]);
				if(altIdx != -1) {
					if(idx != -1) {
						printErr << "more than one wave name or alternative wave name is matching wave in fit result for wave '" << fitInformation->getWave(idxWave).waveName() << "'." << std::endl;
						throw;
					}
					idx = altIdx;
				}
			}
			if(idx == -1) {
				printErr << "wave '" << fitInformation->getWave(idxWave).waveName() << "' not in fit result." << std::endl;
				throw;
			}
			waveNames[idxWave] = fit->waveName(idx);
		}
	}


	void
	readFitResultMatrices(const rpwa::resonanceFit::informationConstPtr& fitInformation,
	                      TTree* tree,
	                      rpwa::fitResult* fit,
	                      const std::vector<Long64_t>& mapping,
	                      const double rescaleErrors,
	                      const boost::multi_array<std::string, 1>& waveNames,
	                      boost::multi_array<std::complex<double>, 2>& productionAmplitudes,
	                      boost::multi_array<TMatrixT<double>, 1>& productionAmplitudesCovariance,
	                      boost::multi_array<std::complex<double>, 3>& spinDensityMatrices,
	                      boost::multi_array<TMatrixT<double>, 1>& spinDensityCovarianceMatrices,
	                      boost::multi_array<std::pair<double, double>, 2>& plottingIntensities,
	                      boost::multi_array<std::pair<double, double>, 3>& plottingSpinDensityMatrixElementsReal,
	                      boost::multi_array<std::pair<double, double>, 3>& plottingSpinDensityMatrixElementsImag,
	                      boost::multi_array<std::pair<double, double>, 3>& plottingPhases)
	{
		if(not tree or not fit) {
			printErr << "'tree' or 'fit' is not a pointer to a valid object." << std::endl;
			throw;
		}

		if(debug) {
			printDebug << "reading spin-density matrices for " << fitInformation->nrWaves() << " waves from fit result." << std::endl;
		}

		productionAmplitudes.resize(boost::extents[mapping.size()][waveNames.size()]);
		productionAmplitudesCovariance.resize(boost::extents[mapping.size()]);

		spinDensityMatrices.resize(boost::extents[mapping.size()][waveNames.size()][waveNames.size()]);
		spinDensityCovarianceMatrices.resize(boost::extents[mapping.size()]);

		plottingIntensities.resize(boost::extents[mapping.size()][waveNames.size()]);
		plottingSpinDensityMatrixElementsReal.resize(boost::extents[mapping.size()][waveNames.size()][waveNames.size()]);
		plottingSpinDensityMatrixElementsImag.resize(boost::extents[mapping.size()][waveNames.size()][waveNames.size()]);
		plottingPhases.resize(boost::extents[mapping.size()][waveNames.size()][waveNames.size()]);

		for(size_t idxMass = 0; idxMass < mapping.size(); ++idxMass) {
			if(debug) {
				printDebug << "reading entry " << mapping[idxMass] << " for mass bin " << idxMass << " from tree." << std::endl;
			}
			// FIXME: in case of reading the fit result for a systematic tree this might happen, so this should be allowed in certain cases
			if(tree->GetEntry(mapping[idxMass]) == 0) {
				printErr << "error while reading entry " << mapping[idxMass] << " from tree." << std::endl;
				throw;
			}

			spinDensityCovarianceMatrices[idxMass].ResizeTo(waveNames.size() * (waveNames.size() + 1), waveNames.size() * (waveNames.size() + 1));
			for(size_t idxWave = 0; idxWave < waveNames.size(); ++idxWave) {
				const int idx = fit->waveIndex(waveNames[idxWave]);
				if(idx == -1) {
					printErr << "wave '" << fitInformation->getWave(idxWave).waveName() << "' not in fit result." << std::endl;
					throw;
				}

				plottingIntensities[idxMass][idxWave] = std::make_pair(fit->intensity(idx),
				                                                       fit->intensityErr(idx) * sqrt(rescaleErrors));

				for(size_t jdxWave = 0; jdxWave < waveNames.size(); ++jdxWave) {
					const int jdx = fit->waveIndex(waveNames[jdxWave]);
					if(jdx == -1) {
						printErr << "wave '" << fitInformation->getWave(jdxWave).waveName() << "' not in fit result." << std::endl;
						throw;
					}

					plottingSpinDensityMatrixElementsReal[idxMass][idxWave][jdxWave] = std::make_pair(fit->spinDensityMatrixElem(idx, jdx).real(),
					                                                                                  sqrt(fit->spinDensityMatrixElemCov(idx, jdx)(0, 0) * rescaleErrors));
					plottingSpinDensityMatrixElementsImag[idxMass][idxWave][jdxWave] = std::make_pair(fit->spinDensityMatrixElem(idx, jdx).imag(),
					                                                                                  sqrt(fit->spinDensityMatrixElemCov(idx, jdx)(1, 1) * rescaleErrors));
					plottingPhases[idxMass][idxWave][jdxWave] = std::make_pair(fit->phase(idx, jdx),
					                                                           fit->phaseErr(idx, jdx) * sqrt(rescaleErrors));

					spinDensityMatrices[idxMass][idxWave][jdxWave] = fit->spinDensityMatrixElem(idx, jdx);

					if(jdxWave >= idxWave) {
						const TMatrixT<double> spinDensityMatrixElemCov = fit->spinDensityMatrixElemCov(idx, jdx) * rescaleErrors;

						const size_t idxCov = waveNames.size() * (waveNames.size() + 1) - (waveNames.size() - idxWave) * ( waveNames.size() - idxWave + 1) + 2 * (jdxWave - idxWave);
						spinDensityCovarianceMatrices[idxMass].SetSub(idxCov, idxCov, spinDensityMatrixElemCov);
					}
				}
			}

			// for the production amplitudes loop over the production
			// amplitudes of the fit result
			std::vector<unsigned int> prodAmpIndicesForCov(waveNames.size());
			for(unsigned int idxProdAmp = 0; idxProdAmp < fit->nmbProdAmps(); ++idxProdAmp) {
				const std::string waveName = fit->waveNameForProdAmp(idxProdAmp);

				const boost::multi_array<std::string, 1>::const_iterator it = std::find(waveNames.begin(), waveNames.end(), waveName);
				// most of the waves are ignored
				if(it == waveNames.end()) {
					continue;
				}
				size_t idxWave = it - waveNames.begin();

				int rank = fit->rankOfProdAmp(idxProdAmp);
				// TODO: multiple ranks, in that case also check that rank is not -1
				if(rank != 0) {
					printErr << "can only handle rank-1 fit (production amplitude '" << fit->prodAmpName(idxProdAmp)
					         << "' of wave '" << waveName << "' has rank " << rank << ")." << std::endl;
					throw;
				}

				productionAmplitudes[idxMass][idxWave] = fit->prodAmp(idxProdAmp);

				prodAmpIndicesForCov[idxWave] = idxProdAmp;
			}
			const TMatrixT<double> prodAmpCov = fit->prodAmpCov(prodAmpIndicesForCov) * rescaleErrors;
			productionAmplitudesCovariance[idxMass].ResizeTo(prodAmpCov);
			productionAmplitudesCovariance[idxMass] = prodAmpCov;

			if(debug) {
				std::ostringstream outputProdAmp;
				std::ostringstream outputProdAmpCovariance;
				std::ostringstream output;
				std::ostringstream outputCovariance;

				outputProdAmp << " (";
				for(size_t idxWave = 0; idxWave < waveNames.size(); ++idxWave) {
					outputProdAmp << " " << productionAmplitudes[idxMass][idxWave];

					outputProdAmpCovariance << " (";
					output << " (";
					for(size_t jdxWave = 0; jdxWave < waveNames.size(); ++jdxWave) {
						output << " " << spinDensityMatrices[idxMass][idxWave][jdxWave];

						outputProdAmpCovariance << " (";
						for(size_t idx = 0; idx < 2; ++idx) {
							outputProdAmpCovariance << " (";
							for(size_t jdx = 0; jdx < 2; ++jdx) {
								outputProdAmpCovariance << " " << productionAmplitudesCovariance[idxMass](idxWave+idx, jdxWave+jdx);
							}
							outputProdAmpCovariance << " )";
						}
						outputProdAmpCovariance << " )";

						if(jdxWave >= idxWave) {
							const size_t idxCov = waveNames.size()*(waveNames.size()+1) - (waveNames.size()-idxWave)*(waveNames.size()-idxWave+1) + 2*(jdxWave-idxWave);
							outputCovariance << " (";
							for(size_t idx = 0; idx < 2; ++idx) {
								outputCovariance << " (";
								for(size_t jdx = 0; jdx < 2; ++jdx) {
									outputCovariance << " " << spinDensityCovarianceMatrices[idxMass](idxCov+idx, idxCov+jdx);
								}
								outputCovariance << " )";
							}
							outputCovariance << " )";
						}
					}
					outputProdAmpCovariance << " )";
					output << " )";
				}
				outputProdAmp << " )";

				printDebug << "production amplitudes: " << outputProdAmp.str() << std::endl;
				printDebug << "production amplitudes covariances: " << outputProdAmpCovariance.str() << std::endl;
				printDebug << "spin-density matrix: " << output.str() << std::endl;
				printDebug << "spin-density covariance matrix (diagonal elements): " << outputCovariance.str() << std::endl;
			}
		} // end loop over mass bins
	}


	void
	readFitResultIntegrals(const rpwa::resonanceFit::informationConstPtr& fitInformation,
	                       TTree* tree,
	                       rpwa::fitResult* fit,
	                       const std::vector<Long64_t>& mapping,
	                       const boost::multi_array<std::string, 1>& waveNames,
	                       boost::multi_array<double, 2>& phaseSpaceIntegrals)
	{
		if(not tree or not fit) {
			printErr << "'tree' or 'fit' is not a pointer to a valid object." << std::endl;
			throw;
		}

		phaseSpaceIntegrals.resize(boost::extents[mapping.size()][waveNames.size()]);

		if(debug) {
			printDebug << "reading phase-space integrals for " << waveNames.size() << " waves from fit result." << std::endl;
		}

		for(size_t idxMass = 0; idxMass < mapping.size(); ++idxMass) {
			if(debug) {
				printDebug << "reading entry " << mapping[idxMass] << " for mass bin " << idxMass << " from tree." << std::endl;
			}
			if(tree->GetEntry(mapping[idxMass]) == 0) {
				printErr << "error while reading entry " << mapping[idxMass] << " from tree." << std::endl;
				throw;
			}

			for(size_t idxWave = 0; idxWave < waveNames.size(); ++idxWave) {
				const double ps = fit->phaseSpaceIntegral(waveNames[idxWave]);
				phaseSpaceIntegrals[idxMass][idxWave] = ps;
			}
		}

		if(debug) {
			for(size_t idxWave = 0; idxWave < waveNames.size(); ++idxWave) {
				std::ostringstream output;
				for(size_t idxMass = 0; idxMass < mapping.size(); ++idxMass) {
					output << " " << phaseSpaceIntegrals[idxMass][idxWave];
				}
				printDebug << "phase-space integrals for wave '" << fitInformation->getWave(idxWave).waveName() << "' (" << idxWave << "):" << output.str() << std::endl;
			}
		}
	}


	void
	readInFile(const rpwa::resonanceFit::informationConstPtr& fitInformation,
	           const rpwa::resonanceFit::information::bin& bin,
	           boost::multi_array<std::string, 1>& waveNames,
	           size_t& nrMassBins,
	           boost::multi_array<double, 1>& massBinCenters,
	           boost::multi_array<double, 2>& phaseSpaceIntegrals,
	           boost::multi_array<std::complex<double>, 2>& productionAmplitudes,
	           boost::multi_array<TMatrixT<double>, 1>& productionAmplitudesCovariance,
	           boost::multi_array<std::complex<double>, 3>& spinDensityMatrices,
	           boost::multi_array<TMatrixT<double>, 1>& spinDensityMatricesCovariance,
	           boost::multi_array<std::pair<double, double>, 2>& plottingIntensities,
	           boost::multi_array<std::pair<double, double>, 3>& plottingSpinDensityMatrixElementsReal,
	           boost::multi_array<std::pair<double, double>, 3>& plottingSpinDensityMatrixElementsImag,
	           boost::multi_array<std::pair<double, double>, 3>& plottingPhases,
	           const std::string& valTreeName = "pwa",
	           const std::string& valBranchName = "fitResult_v2")
	{
		if(debug) {
			printDebug << "reading fit result from file '" << bin.fileName() << "'." << std::endl;
		}

		std::unique_ptr<TFile> inFile(TFile::Open(bin.fileName().c_str()));
		if(not inFile) {
			printErr << "input file '" << bin.fileName() << "' not found."<< std::endl;
			throw;
		}
		if(inFile->IsZombie()) {
			printErr << "error while reading input file '" << bin.fileName() << "'."<< std::endl;
			throw;
		}

		if(debug) {
			printDebug << "searching for tree '" << valTreeName << "' in file '" << bin.fileName() << "'." << std::endl;
		}

		TTree* inTree;
		inFile->GetObject(valTreeName.c_str(), inTree);
		if(not inTree) {
			printErr << "input tree '" << valTreeName << "' not found in input file '" << bin.fileName() << "'."<< std::endl;
			throw;
		}

		if(debug) {
			printDebug << "searching for branch '" << valBranchName << "' in tree '" << valTreeName << "'." << std::endl;
		}

		rpwa::fitResult* inFit = NULL;
		if(inTree->SetBranchAddress(valBranchName.c_str(), &inFit)) {
			printErr << "branch '" << valBranchName << "' not found in input tree '" << valTreeName << "'." << std::endl;
			throw;
		}

		readFitResultMassBins(inTree,
		                      inFit,
		                      nrMassBins,
		                      massBinCenters);

		std::vector<Long64_t> inMapping;
		checkFitResultMassBins(inTree,
		                       inFit,
		                       nrMassBins,
		                       massBinCenters,
		                       inMapping);

		readFitResultWaveNames(fitInformation,
		                       inTree,
		                       inFit,
		                       waveNames);

		readFitResultMatrices(fitInformation,
		                      inTree,
		                      inFit,
		                      inMapping,
		                      bin.rescaleErrors(),
		                      waveNames,
		                      productionAmplitudes,
		                      productionAmplitudesCovariance,
		                      spinDensityMatrices,
		                      spinDensityMatricesCovariance,
		                      plottingIntensities,
		                      plottingSpinDensityMatrixElementsReal,
		                      plottingSpinDensityMatrixElementsImag,
		                      plottingPhases);

		readFitResultIntegrals(fitInformation,
		                       inTree,
		                       inFit,
		                       inMapping,
		                       waveNames,
		                       phaseSpaceIntegrals);
	}


	void
	readSystematicsFile(const rpwa::resonanceFit::informationConstPtr& fitInformation,
	                    const size_t idxBin,
	                    const rpwa::resonanceFit::information::bin& bin,
	                    const size_t idxSystematics,
	                    const boost::multi_array<std::string, 1>& waveNames,
	                    const size_t nrMassBins,
	                    const boost::multi_array<double, 1>& massBinCenters,
	                    const boost::multi_array<std::pair<double, double>, 3>& plottingPhases,
	                    boost::multi_array<std::pair<double, double>, 2>& sysPlottingIntensities,
	                    boost::multi_array<std::pair<double, double>, 3>& sysPlottingSpinDensityMatrixElementsReal,
	                    boost::multi_array<std::pair<double, double>, 3>& sysPlottingSpinDensityMatrixElementsImag,
	                    boost::multi_array<std::pair<double, double>, 3>& sysPlottingPhases,
	                    const std::string& valTreeName = "pwa",
	                    const std::string& valBranchName = "fitResult_v2")
	{
		if(debug) {
			printDebug << "reading fit result for systematics for bin " << idxBin << " from file at index " << idxSystematics << ": '" << bin.sysFileNames()[idxSystematics] << "'." << std::endl;
		}

		std::unique_ptr<TFile> sysFile(TFile::Open(bin.sysFileNames()[idxSystematics].c_str()));
		if(not sysFile) {
			printErr << "input file '" << bin.sysFileNames()[idxSystematics] << "' not found."<< std::endl;
			throw;
		}
		if(sysFile->IsZombie()) {
			printErr << "error while reading input file '" << bin.sysFileNames()[idxSystematics] << "'."<< std::endl;
			throw;
		}

		if(debug) {
			printDebug << "searching for tree '" << valTreeName << "' in file '" << bin.sysFileNames()[idxSystematics] << "'." << std::endl;
		}

		TTree* sysTree;
		sysFile->GetObject(valTreeName.c_str(), sysTree);
		if(not sysTree) {
			printErr << "input tree '" << valTreeName << "' not found in input file '" << bin.sysFileNames()[idxSystematics] << "'."<< std::endl;
			throw;
		}

		if(debug) {
			printDebug << "searching for branch '" << valBranchName << "' in tree '" << valTreeName << "'." << std::endl;
		}

		rpwa::fitResult* sysFit = NULL;
		if(sysTree->SetBranchAddress(valBranchName.c_str(), &sysFit)) {
			printErr << "branch '" << valBranchName << "' not found in input tree '" << valTreeName << "'." << std::endl;
			throw;
		}

		std::vector<Long64_t> sysMapping;
		checkFitResultMassBins(sysTree,
		                       sysFit,
		                       nrMassBins,
		                       massBinCenters,
		                       sysMapping);

		boost::multi_array<std::complex<double>, 2> tempProductionAmplitudes;
		boost::multi_array<TMatrixT<double>, 1> tempProductionAmplitudesCovariance;
		boost::multi_array<std::complex<double>, 3> tempSpinDensityMatrices;
		boost::multi_array<TMatrixT<double>, 1> tempSpinDensityCovarianceMatrices;
		boost::multi_array<std::pair<double, double>, 2> tempSysPlottingIntensities;
		boost::multi_array<std::pair<double, double>, 3> tempSysPlottingSpinDensityMatrixElementsReal;
		boost::multi_array<std::pair<double, double>, 3> tempSysPlottingSpinDensityMatrixElementsImag;
		boost::multi_array<std::pair<double, double>, 3> tempSysPlottingPhases;
		readFitResultMatrices(fitInformation,
		                      sysTree,
		                      sysFit,
		                      sysMapping,
		                      bin.rescaleErrors(),
		                      waveNames,
		                      tempProductionAmplitudes,
		                      tempProductionAmplitudesCovariance,
		                      tempSpinDensityMatrices,
		                      tempSpinDensityCovarianceMatrices,
		                      tempSysPlottingIntensities,
		                      tempSysPlottingSpinDensityMatrixElementsReal,
		                      tempSysPlottingSpinDensityMatrixElementsImag,
		                      tempSysPlottingPhases);

		for(size_t idxMass = 0; idxMass < nrMassBins; ++idxMass) {
			for(size_t idxWave = 0; idxWave < fitInformation->nrWaves(); ++idxWave) {
				sysPlottingIntensities[idxMass][idxWave].first = std::min(sysPlottingIntensities[idxMass][idxWave].first,
				                                                          tempSysPlottingIntensities[idxMass][idxWave].first);
				sysPlottingIntensities[idxMass][idxWave].second = std::max(sysPlottingIntensities[idxMass][idxWave].second,
				                                                           tempSysPlottingIntensities[idxMass][idxWave].first);

				for(size_t jdxWave = 0; jdxWave < fitInformation->nrWaves(); ++jdxWave) {
					sysPlottingSpinDensityMatrixElementsReal[idxMass][idxWave][jdxWave].first = std::min(sysPlottingSpinDensityMatrixElementsReal[idxMass][idxWave][jdxWave].first,
					                                                                                     tempSysPlottingSpinDensityMatrixElementsReal[idxMass][idxWave][jdxWave].first);
					sysPlottingSpinDensityMatrixElementsReal[idxMass][idxWave][jdxWave].second = std::max(sysPlottingSpinDensityMatrixElementsReal[idxMass][idxWave][jdxWave].second,
					                                                                                      tempSysPlottingSpinDensityMatrixElementsReal[idxMass][idxWave][jdxWave].first);

					sysPlottingSpinDensityMatrixElementsImag[idxMass][idxWave][jdxWave].first = std::min(sysPlottingSpinDensityMatrixElementsImag[idxMass][idxWave][jdxWave].first,
					                                                                                     tempSysPlottingSpinDensityMatrixElementsImag[idxMass][idxWave][jdxWave].first);
					sysPlottingSpinDensityMatrixElementsImag[idxMass][idxWave][jdxWave].second = std::max(sysPlottingSpinDensityMatrixElementsImag[idxMass][idxWave][jdxWave].second,
					                                                                                      tempSysPlottingSpinDensityMatrixElementsImag[idxMass][idxWave][jdxWave].first);

					// rotate phase by +- 360 degrees to be as
					// close as possible to the data used in the
					// fit
					if(std::abs(tempSysPlottingPhases[idxMass][idxWave][jdxWave].first+360. - plottingPhases[idxMass][idxWave][jdxWave].first) < std::abs(tempSysPlottingPhases[idxMass][idxWave][jdxWave].first - plottingPhases[idxMass][idxWave][jdxWave].first)) {
						tempSysPlottingPhases[idxMass][idxWave][jdxWave].first += 360.;
					} else if(std::abs(tempSysPlottingPhases[idxMass][idxWave][jdxWave].first-360. - plottingPhases[idxMass][idxWave][jdxWave].first) < std::abs(tempSysPlottingPhases[idxMass][idxWave][jdxWave].first - plottingPhases[idxMass][idxWave][jdxWave].first)) {
						tempSysPlottingPhases[idxMass][idxWave][jdxWave].first -= 360.;
					}

					sysPlottingPhases[idxMass][idxWave][jdxWave].first = std::min(sysPlottingPhases[idxMass][idxWave][jdxWave].first,
					                                                              tempSysPlottingPhases[idxMass][idxWave][jdxWave].first);
					sysPlottingPhases[idxMass][idxWave][jdxWave].second = std::max(sysPlottingPhases[idxMass][idxWave][jdxWave].second,
					                                                               tempSysPlottingPhases[idxMass][idxWave][jdxWave].first);
				}
			}
		}
	}


	void
	readSystematicsFiles(const rpwa::resonanceFit::informationConstPtr& fitInformation,
	                     const size_t idxBin,
	                     const rpwa::resonanceFit::information::bin& bin,
	                     const boost::multi_array<std::string, 1>& waveNames,
	                     const size_t nrMassBins,
	                     const boost::multi_array<double, 1>& massBinCenters,
	                     const boost::multi_array<std::pair<double, double>, 3>& plottingPhases,
	                     boost::multi_array<std::pair<double, double>, 2>& sysPlottingIntensities,
	                     boost::multi_array<std::pair<double, double>, 3>& sysPlottingSpinDensityMatrixElementsReal,
	                     boost::multi_array<std::pair<double, double>, 3>& sysPlottingSpinDensityMatrixElementsImag,
	                     boost::multi_array<std::pair<double, double>, 3>& sysPlottingPhases,
	                     const std::string& valTreeName = "pwa",
	                     const std::string& valBranchName = "fitResult_v2")
	{
		if(debug) {
			printDebug << "reading fit results for systematic errors for bin " << idxBin << " from " << bin.sysFileNames().size() << " files." << std::endl;
		}

		for(size_t idxSystematics = 0; idxSystematics < bin.sysFileNames().size(); ++idxSystematics) {
			readSystematicsFile(fitInformation,
			                    idxBin,
			                    bin,
			                    idxSystematics,
			                    waveNames,
			                    nrMassBins,
			                    massBinCenters,
			                    plottingPhases,
			                    sysPlottingIntensities,
			                    sysPlottingSpinDensityMatrixElementsReal,
			                    sysPlottingSpinDensityMatrixElementsImag,
			                    sysPlottingPhases,
			                    valTreeName,
			                    valBranchName);
		}
	}


	void
	readInFiles(const rpwa::resonanceFit::informationConstPtr& fitInformation,
	            boost::multi_array<std::string, 2>& waveNames,
	            std::vector<size_t>& nrMassBins,
	            boost::multi_array<double, 2>& massBinCenters,
	            boost::multi_array<double, 3>& phaseSpaceIntegrals,
	            boost::multi_array<std::complex<double>, 3>& productionAmplitudes,
	            boost::multi_array<TMatrixT<double>, 2>& productionAmplitudesCovariance,
	            boost::multi_array<std::complex<double>, 4>& spinDensityMatrices,
	            boost::multi_array<TMatrixT<double>, 2>& spinDensityMatricesCovariance,
	            boost::multi_array<std::pair<double, double>, 3>& plottingIntensities,
	            boost::multi_array<std::pair<double, double>, 4>& plottingSpinDensityMatrixElementsReal,
	            boost::multi_array<std::pair<double, double>, 4>& plottingSpinDensityMatrixElementsImag,
	            boost::multi_array<std::pair<double, double>, 4>& plottingPhases,
	            boost::multi_array<std::pair<double, double>, 3>& sysPlottingIntensities,
	            boost::multi_array<std::pair<double, double>, 4>& sysPlottingSpinDensityMatrixElementsReal,
	            boost::multi_array<std::pair<double, double>, 4>& sysPlottingSpinDensityMatrixElementsImag,
	            boost::multi_array<std::pair<double, double>, 4>& sysPlottingPhases,
	            const std::string& valTreeName = "pwa",
	            const std::string& valBranchName = "fitResult_v2")
	{
		for(size_t idxBin = 0; idxBin < fitInformation->nrBins(); ++idxBin) {
			const rpwa::resonanceFit::information::bin& bin = fitInformation->getBin(idxBin);

			boost::multi_array<std::string, 1> tempWaveNames;
			size_t tempNrMassBins;
			boost::multi_array<double, 1> tempMassBinCenters;
			boost::multi_array<double, 2> tempPhaseSpaceIntegrals;
			boost::multi_array<std::complex<double>, 2> tempProductionAmplitudes;
			boost::multi_array<TMatrixT<double>, 1> tempProductionAmplitudesCovariance;
			boost::multi_array<std::complex<double>, 3> tempSpinDensityMatrices;
			boost::multi_array<TMatrixT<double>, 1> tempSpinDensityMatricesCovariance;
			boost::multi_array<std::pair<double, double>, 2> tempPlottingIntensities;
			boost::multi_array<std::pair<double, double>, 3> tempPlottingSpinDensityMatrixElementsReal;
			boost::multi_array<std::pair<double, double>, 3> tempPlottingSpinDensityMatrixElementsImag;
			boost::multi_array<std::pair<double, double>, 3> tempPlottingPhases;

			readInFile(fitInformation,
			           bin,
			           tempWaveNames,
			           tempNrMassBins,
			           tempMassBinCenters,
			           tempPhaseSpaceIntegrals,
			           tempProductionAmplitudes,
			           tempProductionAmplitudesCovariance,
			           tempSpinDensityMatrices,
			           tempSpinDensityMatricesCovariance,
			           tempPlottingIntensities,
			           tempPlottingSpinDensityMatrixElementsReal,
			           tempPlottingSpinDensityMatrixElementsImag,
			           tempPlottingPhases,
			           valTreeName,
			           valBranchName);

			rpwa::resonanceFit::adjustSizeAndSet(waveNames, idxBin, tempWaveNames);
			rpwa::resonanceFit::adjustSizeAndSet(nrMassBins, idxBin, tempNrMassBins);
			rpwa::resonanceFit::adjustSizeAndSet(massBinCenters, idxBin, tempMassBinCenters);
			rpwa::resonanceFit::adjustSizeAndSet(phaseSpaceIntegrals, idxBin, tempPhaseSpaceIntegrals);
			rpwa::resonanceFit::adjustSizeAndSet(productionAmplitudes, idxBin, tempProductionAmplitudes);
			rpwa::resonanceFit::adjustSizeAndSet(productionAmplitudesCovariance, idxBin, tempProductionAmplitudesCovariance);
			rpwa::resonanceFit::adjustSizeAndSet(spinDensityMatrices, idxBin, tempSpinDensityMatrices);
			rpwa::resonanceFit::adjustSizeAndSet(spinDensityMatricesCovariance, idxBin, tempSpinDensityMatricesCovariance);
			rpwa::resonanceFit::adjustSizeAndSet(plottingIntensities, idxBin, tempPlottingIntensities);
			rpwa::resonanceFit::adjustSizeAndSet(plottingSpinDensityMatrixElementsReal, idxBin, tempPlottingSpinDensityMatrixElementsReal);
			rpwa::resonanceFit::adjustSizeAndSet(plottingSpinDensityMatrixElementsImag, idxBin, tempPlottingSpinDensityMatrixElementsImag);
			rpwa::resonanceFit::adjustSizeAndSet(plottingPhases, idxBin, tempPlottingPhases);

			// extract information for systematic errors
			// initialize with real fit result
			boost::multi_array<std::pair<double, double>, 2> tempSysPlottingIntensities(std::vector<size_t>(tempPlottingIntensities.shape(), tempPlottingIntensities.shape()+tempPlottingIntensities.num_dimensions()));
			boost::multi_array<std::pair<double, double>, 3> tempSysPlottingSpinDensityMatrixElementsReal(std::vector<size_t>(tempPlottingSpinDensityMatrixElementsReal.shape(), tempPlottingSpinDensityMatrixElementsReal.shape()+tempPlottingSpinDensityMatrixElementsReal.num_dimensions()));
			boost::multi_array<std::pair<double, double>, 3> tempSysPlottingSpinDensityMatrixElementsImag(std::vector<size_t>(tempPlottingSpinDensityMatrixElementsImag.shape(), tempPlottingSpinDensityMatrixElementsImag.shape()+tempPlottingSpinDensityMatrixElementsImag.num_dimensions()));
			boost::multi_array<std::pair<double, double>, 3> tempSysPlottingPhases(std::vector<size_t>(tempPlottingPhases.shape(), tempPlottingPhases.shape()+tempPlottingPhases.num_dimensions()));

			for(size_t idxMass = 0; idxMass < nrMassBins[idxBin]; ++idxMass) {
				for(size_t idxWave = 0; idxWave < fitInformation->nrWaves(); ++idxWave) {
					tempSysPlottingIntensities[idxMass][idxWave] = std::make_pair(tempPlottingIntensities[idxMass][idxWave].first,
					                                                              tempPlottingIntensities[idxMass][idxWave].first);

					for(size_t jdxWave = 0; jdxWave < fitInformation->nrWaves(); ++jdxWave) {
						tempSysPlottingSpinDensityMatrixElementsReal[idxMass][idxWave][jdxWave] = std::make_pair(tempPlottingSpinDensityMatrixElementsReal[idxMass][idxWave][jdxWave].first,
						                                                                                         tempPlottingSpinDensityMatrixElementsReal[idxMass][idxWave][jdxWave].first);
						tempSysPlottingSpinDensityMatrixElementsImag[idxMass][idxWave][jdxWave] = std::make_pair(tempPlottingSpinDensityMatrixElementsImag[idxMass][idxWave][jdxWave].first,
						                                                                                         tempPlottingSpinDensityMatrixElementsImag[idxMass][idxWave][jdxWave].first);
						tempSysPlottingPhases[idxMass][idxWave][jdxWave] = std::make_pair(tempPlottingPhases[idxMass][idxWave][jdxWave].first,
						                                                                  tempPlottingPhases[idxMass][idxWave][jdxWave].first);
					}
				}
			}

			if(bin.sysFileNames().size() > 0) {
				readSystematicsFiles(fitInformation,
				                     idxBin,
				                     bin,
				                     waveNames[idxBin],
				                     nrMassBins[idxBin],
				                     massBinCenters[idxBin],
				                     tempPlottingPhases,
				                     tempSysPlottingIntensities,
				                     tempSysPlottingSpinDensityMatrixElementsReal,
				                     tempSysPlottingSpinDensityMatrixElementsImag,
				                     tempSysPlottingPhases,
				                     valTreeName,
				                     valBranchName);
			}

			rpwa::resonanceFit::adjustSizeAndSet(sysPlottingIntensities, idxBin, tempSysPlottingIntensities);
			rpwa::resonanceFit::adjustSizeAndSet(sysPlottingSpinDensityMatrixElementsReal, idxBin, tempSysPlottingSpinDensityMatrixElementsReal);
			rpwa::resonanceFit::adjustSizeAndSet(sysPlottingSpinDensityMatrixElementsImag, idxBin, tempSysPlottingSpinDensityMatrixElementsImag);
			rpwa::resonanceFit::adjustSizeAndSet(sysPlottingPhases, idxBin, tempSysPlottingPhases);
		}
	}


	void
	readModelAnchors(const YAML::Node& configModel,
	                 std::string& anchorWaveName,
	                 std::string& anchorComponentName)
	{
		if(debug) {
			printDebug << "reading 'anchorwave'." << std::endl;
		}

		const YAML::Node& configAnchors = configModel["anchorwave"];
		if(not configAnchors) {
			printErr << "'anchorwave' is not a valid YAML node." << std::endl;
			throw;
		}
		if(not configAnchors.IsMap()) {
			printErr << "'anchorwave' is not a YAML map." << std::endl;
			throw;
		}

		std::map<std::string, rpwa::YamlCppUtils::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("name", rpwa::YamlCppUtils::TypeString)
		                     ("resonance", rpwa::YamlCppUtils::TypeString);
		if(not checkIfAllVariablesAreThere(configAnchors, mandatoryArguments)) {
			printErr << "'anchorwave' does not contain all required variables." << std::endl;
			throw;
		}

		anchorWaveName = configAnchors["name"].as<std::string>();
		anchorComponentName = configAnchors["resonance"].as<std::string>();
	}


	void
	writeModelAnchors(YAML::Emitter& yamlOutput,
	                  const std::string& anchorWaveName,
	                  const std::string& anchorComponentName)
	{
		if(debug) {
			printDebug << "writing 'anchorwave'." << std::endl;
		}

		yamlOutput << YAML::Key << "anchorwave";
		yamlOutput << YAML::Value;

		yamlOutput << YAML::BeginMap;

		yamlOutput << YAML::Key << "name";
		yamlOutput << YAML::Value << anchorWaveName;

		yamlOutput << YAML::Key << "resonance";
		yamlOutput << YAML::Value << anchorComponentName;

		yamlOutput << YAML::EndMap;
	}


	void
	readModelParameter(const YAML::Node& configComponent,
	                   const std::string& parameterName,
	                   rpwa::resonanceFit::parameter& parameter)
	{
		if(debug) {
			printDebug << "reading parameter '" << parameterName << "'." << std::endl;
		}

		const YAML::Node& configParameter = configComponent[parameterName];
		if(not configParameter) {
			printErr << "final-state mass-dependence does not define parameter '" << parameterName << "'." << std::endl;
			throw;
		}

		parameter.setName(parameterName);

		std::map<std::string, rpwa::YamlCppUtils::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("val", rpwa::YamlCppUtils::TypeFloat)
		                     ("fix", rpwa::YamlCppUtils::TypeBoolean);
		if(not checkIfAllVariablesAreThere(configParameter, mandatoryArguments)) {
			printErr << "'" << parameterName << "' of final-state mass-dependence does not contain all required variables." << std::endl;
			throw;
		}

		const double startValue = configParameter["val"].as<double>();
		parameter.setStartValue(startValue);

		if(configParameter["error"]) {
			if(checkVariableType(configParameter["error"], rpwa::YamlCppUtils::TypeFloat)) {
				parameter.setStartError(configParameter["error"].as<double>());
			} else if(checkVariableType(configParameter["error"], rpwa::YamlCppUtils::TypeString) and configParameter["error"].as<std::string>() == "nan") {
				// some systems cannot convert the string 'nan'
				// to a floating point number, so this is to be
				// done manually. in any case this does not
				// really matter as the error is usually not
				// used.
				parameter.setStartError(std::numeric_limits<double>::has_quiet_NaN ? std::numeric_limits<double>::quiet_NaN() : 0.0);
			} else {
				printErr << "variable 'error' for parameter '" << parameterName << "' of final-state mass-dependence defined, but not a floating point number." << std::endl;
				throw;
			}
		}

		parameter.setFixed(configParameter["fix"].as<bool>());

		if(configParameter["lower"]) {
			if(checkVariableType(configParameter["lower"], rpwa::YamlCppUtils::TypeFloat)) {
				parameter.setLimitedLower(true);
				parameter.setLimitLower(configParameter["lower"].as<double>());
			} else {
				printErr << "variable 'lower' for parameter '" << parameterName << "' of final-state mass-dependence defined, but not a floating point number." << std::endl;
				throw;
			}
		} else {
			parameter.setLimitedLower(false);
		}
		if(configParameter["upper"]) {
			if(checkVariableType(configParameter["upper"], rpwa::YamlCppUtils::TypeFloat)) {
				parameter.setLimitedUpper(true);
				parameter.setLimitUpper(configParameter["upper"].as<double>());
			} else {
				printErr << "variable 'upper' for parameter '" << parameterName << "' of final-state mass-dependence defined, but not a floating point number." << std::endl;
				throw;
			}
		} else {
			parameter.setLimitedUpper(false);
		}

		if(configParameter["step"]) {
			if(checkVariableType(configParameter["step"], rpwa::YamlCppUtils::TypeFloat)) {
				parameter.setStep(configParameter["step"].as<double>());
			} else {
				printErr << "variable 'step' for parameter '" << parameterName << "' of final-state mass-dependence defined, but not a floating point number." << std::endl;
				throw;
			}
		}
	}


	void
	writeModelParameter(YAML::Emitter& yamlOutput,
	                    const rpwa::resonanceFit::parameter& parameter,
	                    const double value,
	                    const double error)
	{
		if(debug) {
			printDebug << "writing parameter '" << parameter.name() << "'." << std::endl;
		}

		yamlOutput << YAML::Key << parameter.name();
		yamlOutput << YAML::Value;

		yamlOutput << YAML::BeginMap;

		yamlOutput << YAML::Key << "val";
		yamlOutput << YAML::Value << value;

		yamlOutput << YAML::Key << "error";
		yamlOutput << YAML::Value << error;

		if(parameter.limitedLower()) {
			yamlOutput << YAML::Key << "lower";
			yamlOutput << YAML::Value << parameter.limitLower();
		}

		if(parameter.limitedUpper()) {
			yamlOutput << YAML::Key << "upper";
			yamlOutput << YAML::Value << parameter.limitUpper();
		}

		yamlOutput << YAML::Key << "step";
		yamlOutput << YAML::Value << parameter.step();

		yamlOutput << YAML::Key << "fix";
		yamlOutput << YAML::Value << parameter.fixed();

		yamlOutput << YAML::EndMap;
	}


	double
	readModelComponentDecayChannelBranchingRatio(const YAML::Node& configDecayChannel)
	{
		std::map<std::string, rpwa::YamlCppUtils::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("branchingRatio", rpwa::YamlCppUtils::TypeFloat);
		if(not checkIfAllVariablesAreThere(configDecayChannel, mandatoryArguments)) {
			printErr << "decay channel does not contain all required variables." << std::endl;
			throw;
		}

		return configDecayChannel["branchingRatio"].as<double>();
	}


	template<typename T>
	void
	writeModelComponentDecayChannelBranchingRatio(YAML::Emitter& yamlOutput,
	                                              const std::shared_ptr<const T>& component,
	                                              const size_t idxDecayChannel)
	{
		yamlOutput << YAML::Key << "branchingRatio";
		yamlOutput << YAML::Value << component->branchingRatio()[idxDecayChannel];
	}


	double
	readModelComponentDecayChannelExponent(const YAML::Node& configDecayChannel)
	{
		std::map<std::string, rpwa::YamlCppUtils::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("exponent", rpwa::YamlCppUtils::TypeFloat);
		if(not checkIfAllVariablesAreThere(configDecayChannel, mandatoryArguments)) {
			printErr << "decay channel does not contain all required variables." << std::endl;
			throw;
		}

		return configDecayChannel["exponent"].as<double>();
	}


	double
	readModelComponentDecayChannelExponent(const YAML::Node& configDecayChannel,
	                                       const double defaultExponent)
	{
		if(configDecayChannel["exponent"]) {
			if(checkVariableType(configDecayChannel["exponent"], rpwa::YamlCppUtils::TypeFloat)) {
				return configDecayChannel["exponent"].as<double>();
			} else {
				printErr << "variable 'exponent' defined, but not a floating-point number." << std::endl;
				throw;
			}
		}

		printInfo << "variable 'exponent' not defined, using default value " << defaultExponent << "." << std::endl;
		return defaultExponent;
	}


	template<typename T>
	void
	writeModelComponentDecayChannelExponent(YAML::Emitter& yamlOutput,
	                                        const std::shared_ptr<const T>& component)
	{
		yamlOutput << YAML::Key << "exponent";
		yamlOutput << YAML::Value << component->exponent();
	}


	void
	readModelComponentDecayChannelIntegral(const YAML::Node& configDecayChannel,
	                                       std::vector<double>& masses,
	                                       std::vector<double>& values)
	{
		std::map<std::string, rpwa::YamlCppUtils::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("integral", rpwa::YamlCppUtils::TypeSequence);
		if(not checkIfAllVariablesAreThere(configDecayChannel, mandatoryArguments)) {
			printErr << "decay channel does not contain all required variables." << std::endl;
			throw;
		}

		const YAML::Node& integrals = configDecayChannel["integral"];

		const size_t nrValues = integrals.size();
		if(nrValues < 2) {
			printErr << "phase-space integral has to contain at least two points." << std::endl;
			throw;
		}

		masses.clear();
		values.clear();

		for(size_t idx = 0; idx < nrValues; ++idx) {
			const YAML::Node& integral = integrals[idx];
			if(not integral.IsSequence()) {
				printErr << "phase-space integral has to consist of arrays of two floating point numbers." << std::endl;
				throw;
			}
			if(integral.size() != 2) {
				printErr << "phase-space integral has to consist of arrays of two floating point numbers." << std::endl;
				throw;
			}
			if(not checkVariableType(integral[0], rpwa::YamlCppUtils::TypeFloat)) {
				printErr << "phase-space integral has to consist of arrays of two floating point numbers." << std::endl;
				throw;
			}
			if(not checkVariableType(integral[1], rpwa::YamlCppUtils::TypeFloat)) {
				printErr << "phase-space integral has to consist of arrays of two floating point numbers." << std::endl;
				throw;
			}

			const double mass = integral[0].as<double>();
			const double val = integral[1].as<double>();

			if(masses.size() > 0 and masses.back() > mass) {
				printErr << "masses of phase-space integral have to be strictly ordered." << std::endl;
				throw;
			}

			masses.push_back(mass);
			values.push_back(val);
		}
	}


	template<typename T>
	void
	writeModelComponentDecayChannelIntegral(YAML::Emitter& yamlOutput,
	                                        const std::shared_ptr<const T>& component)
	{
		yamlOutput << YAML::Key << "integral";
		yamlOutput << YAML::Value;

		yamlOutput << YAML::Flow;
		yamlOutput << YAML::BeginSeq;

		const std::vector<double> masses =  component->masses();
		const std::vector<double> values =  component->values();
		for(size_t idx = 0; idx < masses.size(); ++idx) {
			yamlOutput << YAML::BeginSeq;
			yamlOutput << masses[idx];
			yamlOutput << values[idx];
			yamlOutput << YAML::EndSeq;
		}

		yamlOutput << YAML::EndSeq;
		yamlOutput << YAML::Block;
	}


	template<typename T>
	void
	writeModelComponentDecayChannelIntegral(YAML::Emitter& yamlOutput,
	                                        const std::shared_ptr<const T>& component,
	                                        const size_t idxDecayChannel)
	{
		yamlOutput << YAML::Key << "integral";
		yamlOutput << YAML::Value;

		yamlOutput << YAML::Flow;
		yamlOutput << YAML::BeginSeq;

		const std::vector<double> masses =  component->masses()[idxDecayChannel];
		const std::vector<double> values =  component->values()[idxDecayChannel];
		for(size_t idx = 0; idx < masses.size(); ++idx) {
			yamlOutput << YAML::BeginSeq;
			yamlOutput << masses[idx];
			yamlOutput << values[idx];
			yamlOutput << YAML::EndSeq;
		}

		yamlOutput << YAML::EndSeq;
		yamlOutput << YAML::Block;
	}


	double
	readModelComponentDecayChannelMIsobar1(const YAML::Node& configDecayChannel)
	{
		std::map<std::string, rpwa::YamlCppUtils::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("mIsobar1", rpwa::YamlCppUtils::TypeFloat);
		if(not checkIfAllVariablesAreThere(configDecayChannel, mandatoryArguments)) {
			printErr << "decay channel does not contain all required variables." << std::endl;
			throw;
		}

		return configDecayChannel["mIsobar1"].as<double>();
	}


	template<typename T>
	void
	writeModelComponentDecayChannelMIsobar1(YAML::Emitter& yamlOutput,
	                                        const std::shared_ptr<const T>& component)
	{
		yamlOutput << YAML::Key << "mIsobar1";
		yamlOutput << YAML::Value << component->mIsobar1();
	}


	template<typename T>
	void
	writeModelComponentDecayChannelMIsobar1(YAML::Emitter& yamlOutput,
	                                        const std::shared_ptr<const T>& component,
	                                        const size_t idxDecayChannel)
	{
		yamlOutput << YAML::Key << "mIsobar1";
		yamlOutput << YAML::Value << component->mIsobar1()[idxDecayChannel];
	}


	double
	readModelComponentDecayChannelMIsobar2(const YAML::Node& configDecayChannel)
	{
		std::map<std::string, rpwa::YamlCppUtils::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("mIsobar2", rpwa::YamlCppUtils::TypeFloat);
		if(not checkIfAllVariablesAreThere(configDecayChannel, mandatoryArguments)) {
			printErr << "decay channel does not contain all required variables." << std::endl;
			throw;
		}

		return configDecayChannel["mIsobar2"].as<double>();
	}


	template<typename T>
	void
	writeModelComponentDecayChannelMIsobar2(YAML::Emitter& yamlOutput,
	                                        const std::shared_ptr<const T>& component)
	{
		yamlOutput << YAML::Key << "mIsobar2";
		yamlOutput << YAML::Value << component->mIsobar2();
	}


	template<typename T>
	void
	writeModelComponentDecayChannelMIsobar2(YAML::Emitter& yamlOutput,
	                                        const std::shared_ptr<const T>& component,
	                                        const size_t idxDecayChannel)
	{
		yamlOutput << YAML::Key << "mIsobar2";
		yamlOutput << YAML::Value << component->mIsobar2()[idxDecayChannel];
	}


	int
	readModelComponentDecayChannelRelAngularMom(const YAML::Node& configDecayChannel)
	{
		std::map<std::string, rpwa::YamlCppUtils::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("relAngularMom", rpwa::YamlCppUtils::TypeInt);
		if(not checkIfAllVariablesAreThere(configDecayChannel, mandatoryArguments)) {
			printErr << "decay channel does not contain all required variables." << std::endl;
			throw;
		}

		return configDecayChannel["relAngularMom"].as<int>();
	}


	int
	readModelComponentDecayChannelRelAngularMom(const YAML::Node& configDecayChannel,
	                                            const int defaultRelAngularMom)
	{
		if(configDecayChannel["relAngularMom"]) {
			if(checkVariableType(configDecayChannel["relAngularMom"], rpwa::YamlCppUtils::TypeInt)) {
				return configDecayChannel["relAngularMom"].as<int>();
			} else {
				printErr << "variable 'relAngularMom' defined, but not an integer." << std::endl;
				throw;
			}
		}

		printInfo << "variable 'relAngularMom' not defined, using default value " << defaultRelAngularMom << "." << std::endl;
		return defaultRelAngularMom;
	}


	template<typename T>
	void
	writeModelComponentDecayChannelRelAngularMom(YAML::Emitter& yamlOutput,
	                                             const std::shared_ptr<const T>& component)
	{
		yamlOutput << YAML::Key << "relAngularMom";
		yamlOutput << YAML::Value << component->relAngularMom();
	}


	template<typename T>
	void
	writeModelComponentDecayChannelRelAngularMom(YAML::Emitter& yamlOutput,
	                                             const std::shared_ptr<const T>& component,
	                                             const size_t idxDecayChannel)
	{
		yamlOutput << YAML::Key << "relAngularMom";
		yamlOutput << YAML::Value << component->relAngularMom()[idxDecayChannel];
	}


	std::vector<rpwa::resonanceFit::component::channel>
	readModelComponentDecayChannels(const YAML::Node& configComponent,
	                                const std::string& componentName,
	                                const rpwa::resonanceFit::informationConstPtr& fitInformation,
	                                const boost::multi_array<std::string, 2>& waveNames,
	                                const std::vector<size_t>& nrMassBins,
	                                const boost::multi_array<double, 2>& massBinCenters,
	                                const boost::multi_array<double, 3>& phaseSpaceIntegrals)
	{
		if(debug) {
			printDebug << "reading 'decaychannels'." << std::endl;
		}

		std::map<std::string, size_t> waveIndices;
		for(size_t idxWave = 0; idxWave < fitInformation->nrWaves(); ++idxWave) {
			const rpwa::resonanceFit::information::wave& wave = fitInformation->getWave(idxWave);

			waveIndices[wave.waveName()] = idxWave;
			for(size_t idxAlt = 0; idxAlt < wave.waveNameAlternatives().size(); ++idxAlt) {
				waveIndices[wave.waveNameAlternatives()[idxAlt]] = idxWave;
			}
		}

		std::map<std::string, std::vector<size_t> > waveBins;
		for(size_t idxBin = 0; idxBin < fitInformation->nrBins(); ++idxBin) {
			for(size_t idxWave = 0; idxWave < fitInformation->nrWaves(); ++idxWave) {
				waveBins[waveNames[idxBin][idxWave]].push_back(idxBin);
			}
		}

		std::map<std::string, rpwa::YamlCppUtils::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("decaychannels", rpwa::YamlCppUtils::TypeSequence);
		if(not checkIfAllVariablesAreThere(configComponent, mandatoryArguments)) {
			printErr << "'components' entry of component '" << componentName << "' does not contain all required variables." << std::endl;
			throw;
		}

		const YAML::Node& configDecayChannels = configComponent["decaychannels"];

		std::vector<rpwa::resonanceFit::component::channel> decayChannels;
		const size_t nrDecayChannels = configDecayChannels.size();
		for(size_t idxDecayChannel = 0; idxDecayChannel < nrDecayChannels; ++idxDecayChannel) {
			const YAML::Node& configDecayChannel = configDecayChannels[idxDecayChannel];

			std::map<std::string, rpwa::YamlCppUtils::Type> mandatoryArguments;
			boost::assign::insert(mandatoryArguments)
			                     ("amp", rpwa::YamlCppUtils::TypeString);
			if(not checkIfAllVariablesAreThere(configDecayChannel, mandatoryArguments)) {
				printErr << "one of the decay channels of component '" << componentName << "' does not contain all required variables." << std::endl;
				throw;
			}

			const std::string waveName = configDecayChannel["amp"].as<std::string>();

			// check that a wave with this wave is not yet in the decay channels
			for(size_t idxChannel = 0; idxChannel < decayChannels.size(); ++idxChannel) {
				if(decayChannels[idxChannel].getWaveName() == waveName) {
					printErr << "wave '" << waveName << "' defined twice in the decay channels of component '" << componentName << "'." << std::endl;
					throw;
				}
			}

			// get index of wave in array for wave names, phase-space integrals, ...
			const std::map<std::string, size_t>::const_iterator it = waveIndices.find(waveName);
			if(it == waveIndices.end()) {
				printErr << "wave '" << waveName << "' not in fit, but used as decay channel in component '" << componentName << "'." << std::endl;
				throw;
			}
			const size_t waveIdx = it->second;

			// get list of bins this wave is defined in
			const std::map<std::string, std::vector<size_t> >::const_iterator it2 = waveBins.find(waveName);
			if(it2 == waveBins.end()) {
				printErr << "wave '" << waveName << "' not in fit, but used as decay channel in component '" << componentName << "'." << std::endl;
				throw;
			}
			const std::vector<size_t>& binsForWave = it2->second;

			// get a view for the current wave for all bins and all mass bins
			boost::multi_array<double, 3>::const_array_view<2>::type view = phaseSpaceIntegrals[boost::indices[boost::multi_array<double, 3>::index_range()][boost::multi_array<double, 3>::index_range()][waveIdx]];
			decayChannels.push_back(rpwa::resonanceFit::component::channel(waveIdx,
			                                                               waveName,
			                                                               binsForWave,
			                                                               nrMassBins,
			                                                               massBinCenters,
			                                                               view));
		}

		return decayChannels;
	}


	template<typename T>
	rpwa::resonanceFit::componentPtr
	readModelComponent(const YAML::Node& /*configComponent*/,
	                   const rpwa::resonanceFit::informationConstPtr& /*fitInformation*/,
	                   const size_t id,
	                   const std::string& name,
	                   const std::vector<rpwa::resonanceFit::parameter>& parameters,
	                   const std::vector<rpwa::resonanceFit::component::channel>& decayChannels,
	                   const std::vector<size_t>& nrMassBins,
	                   const boost::multi_array<double, 2>& massBinCenters,
	                   const bool useBranchings)
	{
		if(debug) {
			printDebug << "reading component of non-specialized type." << std::endl;
		}

		return std::make_shared<T>(id,
		                           name,
		                           parameters,
		                           decayChannels,
		                           nrMassBins,
		                           massBinCenters,
		                           useBranchings);
	}


	template<>
	rpwa::resonanceFit::componentPtr
	readModelComponent<rpwa::resonanceFit::dynamicWidthBreitWigner>(const YAML::Node& configComponent,
	                                                                const rpwa::resonanceFit::informationConstPtr& /*fitInformation*/,
	                                                                const size_t id,
	                                                                const std::string& name,
	                                                                const std::vector<rpwa::resonanceFit::parameter>& parameters,
	                                                                const std::vector<rpwa::resonanceFit::component::channel>& decayChannels,
	                                                                const std::vector<size_t>& nrMassBins,
	                                                                const boost::multi_array<double, 2>& massBinCenters,
	                                                                const bool useBranchings)
	{
		if(debug) {
			printDebug << "reading component of type 'dynamicWidthBreitWigner'." << std::endl;
		}

		std::vector<double> branchingRatio;
		std::vector<int> relAngularMom;
		std::vector<double> mIsobar1;
		std::vector<double> mIsobar2;

		const YAML::Node& configDecayChannels = configComponent["decaychannels"];
		if(not configDecayChannels) {
			printErr << "a component of type 'dynamicWidthBreitWigner' has no decay channels." << std::endl;
			throw;
		}

		const size_t nrDecayChannels = configDecayChannels.size();
		for(size_t idxDecayChannel = 0; idxDecayChannel < nrDecayChannels; ++idxDecayChannel) {
			const YAML::Node& configDecayChannel = configDecayChannels[idxDecayChannel];

			branchingRatio.push_back(readModelComponentDecayChannelBranchingRatio(configDecayChannel));
			relAngularMom.push_back(readModelComponentDecayChannelRelAngularMom(configDecayChannel));
			mIsobar1.push_back(readModelComponentDecayChannelMIsobar1(configDecayChannel));
			mIsobar2.push_back(readModelComponentDecayChannelMIsobar2(configDecayChannel));
		}

		const YAML::Node& configExtraDecayChannels = configComponent["extradecaychannels"];
		const size_t nrExtraDecayChannels = configExtraDecayChannels ? configExtraDecayChannels.size() : 0;
		for(size_t idxExtraDecayChannel = 0; idxExtraDecayChannel < nrExtraDecayChannels; ++idxExtraDecayChannel) {
			const YAML::Node& configExtraDecayChannel = configExtraDecayChannels[idxExtraDecayChannel];

			branchingRatio.push_back(readModelComponentDecayChannelBranchingRatio(configExtraDecayChannel));
			relAngularMom.push_back(readModelComponentDecayChannelRelAngularMom(configExtraDecayChannel));
			mIsobar1.push_back(readModelComponentDecayChannelMIsobar1(configExtraDecayChannel));
			mIsobar2.push_back(readModelComponentDecayChannelMIsobar2(configExtraDecayChannel));
		}

		return std::make_shared<rpwa::resonanceFit::dynamicWidthBreitWigner>(id,
		                                                                     name,
		                                                                     parameters,
		                                                                     decayChannels,
		                                                                     nrMassBins,
		                                                                     massBinCenters,
		                                                                     useBranchings,
		                                                                     branchingRatio,
		                                                                     relAngularMom,
		                                                                     mIsobar1,
		                                                                     mIsobar2);
	}


	template<>
	rpwa::resonanceFit::componentPtr
	readModelComponent<rpwa::resonanceFit::integralWidthBreitWigner>(const YAML::Node& configComponent,
	                                                                 const rpwa::resonanceFit::informationConstPtr& /*fitInformation*/,
	                                                                 const size_t id,
	                                                                 const std::string& name,
	                                                                 const std::vector<rpwa::resonanceFit::parameter>& parameters,
	                                                                 const std::vector<rpwa::resonanceFit::component::channel>& decayChannels,
	                                                                 const std::vector<size_t>& nrMassBins,
	                                                                 const boost::multi_array<double, 2>& massBinCenters,
	                                                                 const bool useBranchings)
	{
		if(debug) {
			printDebug << "reading component of type 'integralWidthBreitWigner'." << std::endl;
		}

		std::vector<double> branchingRatio;
		std::vector<std::vector<double> > masses;
		std::vector<std::vector<double> > values;

		const YAML::Node& configDecayChannels = configComponent["decaychannels"];
		if(not configDecayChannels) {
			printErr << "a component of type 'integralWidthBreitWigner' has no decay channels." << std::endl;
			throw;
		}

		const size_t nrDecayChannels = configDecayChannels.size();
		for(size_t idxDecayChannel = 0; idxDecayChannel < nrDecayChannels; ++idxDecayChannel) {
			const YAML::Node& configDecayChannel = configDecayChannels[idxDecayChannel];

			branchingRatio.push_back(readModelComponentDecayChannelBranchingRatio(configDecayChannel));

			std::vector<double> tempMasses;
			std::vector<double> tempValues;
			readModelComponentDecayChannelIntegral(configDecayChannel, tempMasses, tempValues);
			masses.push_back(tempMasses);
			values.push_back(tempValues);
		}

		const YAML::Node& configExtraDecayChannels = configComponent["extradecaychannels"];
		const size_t nrExtraDecayChannels = configExtraDecayChannels ? configExtraDecayChannels.size() : 0;
		for(size_t idxExtraDecayChannel = 0; idxExtraDecayChannel < nrExtraDecayChannels; ++idxExtraDecayChannel) {
			const YAML::Node& configExtraDecayChannel = configExtraDecayChannels[idxExtraDecayChannel];

			branchingRatio.push_back(readModelComponentDecayChannelBranchingRatio(configExtraDecayChannel));

			std::vector<double> tempMasses;
			std::vector<double> tempValues;
			readModelComponentDecayChannelIntegral(configExtraDecayChannel, tempMasses, tempValues);
			masses.push_back(tempMasses);
			values.push_back(tempValues);
		}

		return std::make_shared<rpwa::resonanceFit::integralWidthBreitWigner>(id,
		                                                                      name,
		                                                                      parameters,
		                                                                      decayChannels,
		                                                                      nrMassBins,
		                                                                      massBinCenters,
		                                                                      useBranchings,
		                                                                      branchingRatio,
		                                                                      masses,
		                                                                      values);
	}


	template<>
	rpwa::resonanceFit::componentPtr
	readModelComponent<rpwa::resonanceFit::exponentialBackground>(const YAML::Node& configComponent,
	                                                              const rpwa::resonanceFit::informationConstPtr& /*fitInformation*/,
	                                                              const size_t id,
	                                                              const std::string& name,
	                                                              const std::vector<rpwa::resonanceFit::parameter>& parameters,
	                                                              const std::vector<rpwa::resonanceFit::component::channel>& decayChannels,
	                                                              const std::vector<size_t>& nrMassBins,
	                                                              const boost::multi_array<double, 2>& massBinCenters,
	                                                              const bool useBranchings)
	{
		if(debug) {
			printDebug << "reading component of type 'exponentialBackground'." << std::endl;
		}

		const int relAngularMom = readModelComponentDecayChannelRelAngularMom(configComponent, 0);
		const double mIsobar1 = readModelComponentDecayChannelMIsobar1(configComponent);
		const double mIsobar2 = readModelComponentDecayChannelMIsobar2(configComponent);
		const double exponent = readModelComponentDecayChannelExponent(configComponent, 2.0);

		return std::make_shared<rpwa::resonanceFit::exponentialBackground>(id,
		                                                                   name,
		                                                                   parameters,
		                                                                   decayChannels,
		                                                                   nrMassBins,
		                                                                   massBinCenters,
		                                                                   useBranchings,
		                                                                   relAngularMom,
		                                                                   mIsobar1,
		                                                                   mIsobar2,
		                                                                   exponent);
	}


	template<>
	rpwa::resonanceFit::componentPtr
	readModelComponent<rpwa::resonanceFit::tPrimeDependentBackground>(const YAML::Node& configComponent,
	                                                                  const rpwa::resonanceFit::informationConstPtr& fitInformation,
	                                                                  const size_t id,
	                                                                  const std::string& name,
	                                                                  const std::vector<rpwa::resonanceFit::parameter>& parameters,
	                                                                  const std::vector<rpwa::resonanceFit::component::channel>& decayChannels,
	                                                                  const std::vector<size_t>& nrMassBins,
	                                                                  const boost::multi_array<double, 2>& massBinCenters,
	                                                                  const bool useBranchings)
	{
		if(debug) {
			printDebug << "reading component of type 'tPrimeDependentBackground'." << std::endl;
		}

		std::vector<double> tPrimeMeans;
		for(size_t idxBin = 0; idxBin < fitInformation->nrBins(); ++idxBin) {
			tPrimeMeans.push_back(fitInformation->getBin(idxBin).tPrimeMean());
		}

		const int relAngularMom = readModelComponentDecayChannelRelAngularMom(configComponent, 0);
		const double mIsobar1 = readModelComponentDecayChannelMIsobar1(configComponent);
		const double mIsobar2 = readModelComponentDecayChannelMIsobar2(configComponent);
		const double exponent = readModelComponentDecayChannelExponent(configComponent, 2.0);

		return std::make_shared<rpwa::resonanceFit::tPrimeDependentBackground>(id,
		                                                                       name,
		                                                                       parameters,
		                                                                       decayChannels,
		                                                                       nrMassBins,
		                                                                       massBinCenters,
		                                                                       useBranchings,
		                                                                       tPrimeMeans,
		                                                                       relAngularMom,
		                                                                       mIsobar1,
		                                                                       mIsobar2,
		                                                                       exponent);
	}


	template<>
	rpwa::resonanceFit::componentPtr
	readModelComponent<rpwa::resonanceFit::exponentialBackgroundIntegral>(const YAML::Node& configComponent,
	                                                                      const rpwa::resonanceFit::informationConstPtr& /*fitInformation*/,
	                                                                      const size_t id,
	                                                                      const std::string& name,
	                                                                      const std::vector<rpwa::resonanceFit::parameter>& parameters,
	                                                                      const std::vector<rpwa::resonanceFit::component::channel>& decayChannels,
	                                                                      const std::vector<size_t>& nrMassBins,
	                                                                      const boost::multi_array<double, 2>& massBinCenters,
	                                                                      const bool useBranchings)
	{
		if(debug) {
			printDebug << "reading component of type 'exponentialBackgroundIntegral'." << std::endl;
		}

		std::vector<double> masses;
		std::vector<double> values;
		readModelComponentDecayChannelIntegral(configComponent, masses, values);
		const double exponent = readModelComponentDecayChannelExponent(configComponent);

		return std::make_shared<rpwa::resonanceFit::exponentialBackgroundIntegral>(id,
		                                                                           name,
		                                                                           parameters,
		                                                                           decayChannels,
		                                                                           nrMassBins,
		                                                                           massBinCenters,
		                                                                           useBranchings,
		                                                                           masses,
		                                                                           values,
		                                                                           exponent);
	}


	template<>
	rpwa::resonanceFit::componentPtr
	readModelComponent<rpwa::resonanceFit::tPrimeDependentBackgroundIntegral>(const YAML::Node& configComponent,
	                                                                          const rpwa::resonanceFit::informationConstPtr& fitInformation,
	                                                                          const size_t id,
	                                                                          const std::string& name,
	                                                                          const std::vector<rpwa::resonanceFit::parameter>& parameters,
	                                                                          const std::vector<rpwa::resonanceFit::component::channel>& decayChannels,
	                                                                          const std::vector<size_t>& nrMassBins,
	                                                                          const boost::multi_array<double, 2>& massBinCenters,
	                                                                          const bool useBranchings)
	{
		if(debug) {
			printDebug << "reading component of type 'tPrimeDependentBackgroundIntegral'." << std::endl;
		}

		std::vector<double> tPrimeMeans;
		for(size_t idxBin = 0; idxBin < fitInformation->nrBins(); ++idxBin) {
			tPrimeMeans.push_back(fitInformation->getBin(idxBin).tPrimeMean());
		}

		std::vector<double> masses;
		std::vector<double> values;
		readModelComponentDecayChannelIntegral(configComponent, masses, values);
		const double exponent = readModelComponentDecayChannelExponent(configComponent);

		return std::make_shared<rpwa::resonanceFit::tPrimeDependentBackgroundIntegral>(id,
		                                                                               name,
		                                                                               parameters,
		                                                                               decayChannels,
		                                                                               nrMassBins,
		                                                                               massBinCenters,
		                                                                               useBranchings,
		                                                                               tPrimeMeans,
		                                                                               masses,
		                                                                               values,
		                                                                               exponent);
	}


	template<typename T>
	rpwa::resonanceFit::componentPtr
	readModelComponent(const YAML::Node& configComponent,
	                   const rpwa::resonanceFit::informationConstPtr& fitInformation,
	                   const size_t id,
	                   const boost::multi_array<std::string, 2>& waveNames,
	                   const std::vector<size_t>& nrMassBins,
	                   const boost::multi_array<double, 2>& massBinCenters,
	                   const boost::multi_array<double, 3>& phaseSpaceIntegrals,
	                   const bool useBranchings)
	{
		if(debug) {
			printDebug << "reading component of type 'component'." << std::endl;
		}

		std::map<std::string, rpwa::YamlCppUtils::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("name", rpwa::YamlCppUtils::TypeString);
		if(not checkIfAllVariablesAreThere(configComponent, mandatoryArguments)) {
			printErr << "'components' entry does not contain all required variables." << std::endl;
			throw;
		}

		const std::string name = configComponent["name"].as<std::string>();

		std::vector<rpwa::resonanceFit::parameter> parameters = T::getDefaultParameters();
		for(size_t idxParameter = 0; idxParameter < parameters.size(); ++idxParameter) {
			const std::string parameterName = parameters[idxParameter].name();
			readModelParameter(configComponent,
			                   parameterName,
			                   parameters[idxParameter]);
		}

		const std::vector<rpwa::resonanceFit::component::channel> decayChannels = readModelComponentDecayChannels(configComponent,
		                                                                                                          name,
		                                                                                                          fitInformation,
		                                                                                                          waveNames,
		                                                                                                          nrMassBins,
		                                                                                                          massBinCenters,
		                                                                                                          phaseSpaceIntegrals);

		return readModelComponent<T>(configComponent,
		                             fitInformation,
		                             id,
		                             name,
		                             parameters,
		                             decayChannels,
		                             nrMassBins,
		                             massBinCenters,
		                             useBranchings);
	}


	rpwa::resonanceFit::componentPtr
	readModelComponent(const YAML::Node& configComponent,
	                   const rpwa::resonanceFit::informationConstPtr& fitInformation,
	                   const size_t id,
	                   const boost::multi_array<std::string, 2>& waveNames,
	                   const std::vector<size_t>& nrMassBins,
	                   const boost::multi_array<double, 2>& massBinCenters,
	                   const boost::multi_array<double, 3>& phaseSpaceIntegrals,
	                   const bool useBranchings)
	{
		if(debug) {
			printDebug << "reading component." << std::endl;
		}

		std::string type = "fixedWidthBreitWigner";
		if(configComponent["type"]) {
			if(not checkVariableType(configComponent["type"], rpwa::YamlCppUtils::TypeString)) {
				printErr << "a component has a type that is not a string." << std::endl;
				throw;
			}
			type = configComponent["type"].as<std::string>();
		}

		if(debug) {
			printDebug << "found component of type '" << type << "'." << std::endl;
		}

		rpwa::resonanceFit::componentPtr component;
		if(type == "fixedWidthBreitWigner") {
			component = readModelComponent<rpwa::resonanceFit::fixedWidthBreitWigner>(configComponent,
			                                                                          fitInformation,
			                                                                          id,
			                                                                          waveNames,
			                                                                          nrMassBins,
			                                                                          massBinCenters,
			                                                                          phaseSpaceIntegrals,
			                                                                          useBranchings);
		} else if(type == "dynamicWidthBreitWigner") {
			component = readModelComponent<rpwa::resonanceFit::dynamicWidthBreitWigner>(configComponent,
			                                                                            fitInformation,
			                                                                            id,
			                                                                            waveNames,
			                                                                            nrMassBins,
			                                                                            massBinCenters,
			                                                                            phaseSpaceIntegrals,
			                                                                            useBranchings);
		} else if(type == "integralWidthBreitWigner") {
			component = readModelComponent<rpwa::resonanceFit::integralWidthBreitWigner>(configComponent,
			                                                                             fitInformation,
			                                                                             id,
			                                                                             waveNames,
			                                                                             nrMassBins,
			                                                                             massBinCenters,
			                                                                             phaseSpaceIntegrals,
			                                                                             useBranchings);
		} else if(type == "constantBackground") {
			component = readModelComponent<rpwa::resonanceFit::constantBackground>(configComponent,
			                                                                       fitInformation,
			                                                                       id,
			                                                                       waveNames,
			                                                                       nrMassBins,
			                                                                       massBinCenters,
			                                                                       phaseSpaceIntegrals,
			                                                                       useBranchings);
		} else if(type == "exponentialBackground") {
			component = readModelComponent<rpwa::resonanceFit::exponentialBackground>(configComponent,
			                                                                          fitInformation,
			                                                                          id,
			                                                                          waveNames,
			                                                                          nrMassBins,
			                                                                          massBinCenters,
			                                                                          phaseSpaceIntegrals,
			                                                                          useBranchings);
		} else if(type == "tPrimeDependentBackground") {
			component = readModelComponent<rpwa::resonanceFit::tPrimeDependentBackground>(configComponent,
			                                                                              fitInformation,
			                                                                              id,
			                                                                              waveNames,
			                                                                              nrMassBins,
			                                                                              massBinCenters,
			                                                                              phaseSpaceIntegrals,
			                                                                              useBranchings);
		} else if(type == "exponentialBackgroundIntegral") {
			component = readModelComponent<rpwa::resonanceFit::exponentialBackgroundIntegral>(configComponent,
			                                                                                  fitInformation,
			                                                                                  id,
			                                                                                  waveNames,
			                                                                                  nrMassBins,
			                                                                                  massBinCenters,
			                                                                                  phaseSpaceIntegrals,
			                                                                                  useBranchings);
		} else if(type == "tPrimeDependentBackgroundIntegral") {
			component = readModelComponent<rpwa::resonanceFit::tPrimeDependentBackgroundIntegral>(configComponent,
			                                                                                      fitInformation,
			                                                                                      id,
			                                                                                      waveNames,
			                                                                                      nrMassBins,
			                                                                                      massBinCenters,
			                                                                                      phaseSpaceIntegrals,
			                                                                                      useBranchings);
		} else {
			printErr << "unknown type '" << type << "'." << std::endl;
			throw;
		}

		if(debug) {
			component->print(printDebug);
		}

		return component;
	}


	rpwa::resonanceFit::componentPtr
	readModelComponent(const YAML::Node& configComponent,
	                   const rpwa::resonanceFit::informationConstPtr& fitInformation,
	                   rpwa::resonanceFit::parameters& fitParameters,
	                   rpwa::resonanceFit::parameters& fitParametersError,
	                   const size_t id,
	                   const boost::multi_array<std::string, 2>& waveNames,
	                   const std::vector<size_t>& nrMassBins,
	                   const boost::multi_array<double, 2>& massBinCenters,
	                   const boost::multi_array<double, 3>& phaseSpaceIntegrals,
	                   const bool useBranchings)
	{
		if(debug) {
			printDebug << "reading component and its parameters." << std::endl;
		}

		const rpwa::resonanceFit::componentPtr& component = readModelComponent(configComponent,
		                                                                       fitInformation,
		                                                                       id,
		                                                                       waveNames,
		                                                                       nrMassBins,
		                                                                       massBinCenters,
		                                                                       phaseSpaceIntegrals,
		                                                                       useBranchings);

		fitParameters.resize(id+1, component->getNrChannels(), component->getNrParameters(), fitInformation->nrBins());
		fitParametersError.resize(id+1, component->getNrChannels(), component->getNrParameters(), fitInformation->nrBins());

		for(size_t idxParameter = 0; idxParameter < component->getNrParameters(); ++idxParameter) {
			fitParameters.setParameter(id, idxParameter, component->getParameter(idxParameter).startValue());
			fitParametersError.setParameter(id, idxParameter, component->getParameter(idxParameter).startError());
		}

		for(size_t idxDecayChannel = 0; idxDecayChannel < component->getNrChannels(); ++idxDecayChannel) {
			const rpwa::resonanceFit::component::channel& channel = component->getChannel(idxDecayChannel);
			const YAML::Node& configDecayChannel = configComponent["decaychannels"][idxDecayChannel];

			if(component->mapCouplingToMasterChannel(component->mapChannelToCoupling(idxDecayChannel)) == idxDecayChannel) {
				const YAML::Node& configCouplings = configDecayChannel["couplings"];
				if(not configCouplings) {
					printErr << "decay channel '" << channel.getWaveName() << "' of component '" << component->getName() << "' has no couplings." << std::endl;
					throw;
				}
				if(not configCouplings.IsSequence()) {
					printErr << "decay channel '" << channel.getWaveName() << "' of component '" << component->getName() << "' has no couplings." << std::endl;
					throw;
				}

				const size_t nrCouplings = configCouplings.size();
				if(nrCouplings != channel.getBins().size()) {
					printErr << "decay channel '" << channel.getWaveName() << "' of component '" << component->getName() << "' has " << nrCouplings << " couplings, not " << channel.getBins().size() << "." << std::endl;
					throw;
				}

				for(size_t idxCoupling = 0; idxCoupling < nrCouplings; ++idxCoupling) {
					const YAML::Node& configCoupling = configCouplings[idxCoupling];
					if(not checkVariableType(configCoupling, rpwa::YamlCppUtils::TypeSequence)) {
						printErr << "one of the couplings of the decay channel '" << channel.getWaveName() << "' of the component '" << component->getName() << "' is not a YAML sequence." << std::endl;
						throw;
					}
					if(configCoupling.size() != 2) {
						printErr << "one of the couplings of the decay channel '" << channel.getWaveName() << "' of the component '" << component->getName() << "' does not contain exactly two entries." << std::endl;
						throw;
					}

					if(not checkVariableType(configCoupling[0], rpwa::YamlCppUtils::TypeFloat)) {
						printErr << "real part of one of the couplings of the decay channel '" << channel.getWaveName() << "' of the component '" << component->getName() << "' is not a floating point number." << std::endl;
						throw;
					}
					if(not checkVariableType(configCoupling[1], rpwa::YamlCppUtils::TypeFloat)) {
						printErr << "imaginary part of one of the couplings of the decay channel '" << channel.getWaveName() << "' of the component '" << component->getName() << "' is not a floating point number." << std::endl;
						throw;
					}

					const double couplingReal = configCoupling[0].as<double>();
					const double couplingImag = configCoupling[1].as<double>();

					const std::complex<double> coupling(couplingReal, couplingImag);
					fitParameters.setCoupling(id, component->mapChannelToCoupling(idxDecayChannel), channel.getBins()[idxCoupling], coupling);
				}
			}

			if(component->getNrBranchings() > 1) {
				if(component->mapBranchingToMasterChannel(component->mapChannelToBranching(idxDecayChannel)) == idxDecayChannel) {
					const YAML::Node& configBranching = configDecayChannel["branching"];
					if(not configBranching) {
						printErr << "decay channel '" << channel.getWaveName() << "' of component '" << component->getName() << "' has no branching." << std::endl;
						throw;
					}
					if(not configBranching.IsSequence()) {
						printErr << "decay channel '" << channel.getWaveName() << "' of component '" << component->getName() << "' has no branching." << std::endl;
						throw;
					}

					if(configBranching.size() != 2) {
						printErr << "branching of the decay channel '" << channel.getWaveName() << "' of the component '" << component->getName() << "' does not contain exactly two entries." << std::endl;
						throw;
					}

					if(not checkVariableType(configBranching[0], rpwa::YamlCppUtils::TypeFloat)) {
						printErr << "real part of branching of the decay channel '" << channel.getWaveName() << "' of the component '" << component->getName() << "' is not a floating point number." << std::endl;
						throw;
					}
					if(not checkVariableType(configBranching[1], rpwa::YamlCppUtils::TypeFloat)) {
						printErr << "imaginary part of branching of the decay channel '" << channel.getWaveName() << "' of the component '" << component->getName() << "' is not a floating point number." << std::endl;
						throw;
					}

					double branchingReal = configBranching[0].as<double>();
					double branchingImag = configBranching[1].as<double>();

					if(component->isBranchingFixed(component->mapChannelToBranching(idxDecayChannel))) {
						// the first branching should always be 1.
						if(branchingReal != 1.0 or branchingImag != 0.0) {
							printWarn << "branching of the decay channel '" << channel.getWaveName() << "' of the component '" << component->getName() << "' forced to 1." << std::endl;
							branchingReal = 1.0;
							branchingImag = 0.0;
						}
					}

					const std::complex<double> branching(branchingReal, branchingImag);
					fitParameters.setBranching(id, component->mapChannelToBranching(idxDecayChannel), branching);
				}
			} else {
				const std::complex<double> branching(1.0, 0.0);
				fitParameters.setBranching(id, 0, branching);
			}
		}

		return component;
	}


	void
	writeModelComponent(YAML::Emitter& yamlOutput,
	                    const rpwa::resonanceFit::componentConstPtr& component,
	                    const rpwa::resonanceFit::parameters& fitParameters,
	                    const rpwa::resonanceFit::parameters& fitParametersError)
	{
		if(debug) {
			printDebug << "writing component and its parameters." << std::endl;
		}

		yamlOutput << YAML::BeginMap;

		yamlOutput << YAML::Key << "name";
		yamlOutput << YAML::Value << component->getName();

		yamlOutput << YAML::Key << "type";
		yamlOutput << YAML::Value << component->getType();

		for(size_t idxParameter = 0; idxParameter < component->getNrParameters(); ++idxParameter) {
			writeModelParameter(yamlOutput,
			                    component->getParameter(idxParameter),
			                    fitParameters.getParameter(component->getId(), idxParameter),
			                    fitParametersError.getParameter(component->getId(), idxParameter));
		}

		yamlOutput << YAML::Key << "decaychannels";
		yamlOutput << YAML::Value;
		yamlOutput << YAML::BeginSeq;

		for(size_t idxDecayChannel = 0; idxDecayChannel < component->getNrChannels(); ++idxDecayChannel) {
			yamlOutput << YAML::BeginMap;

			yamlOutput << YAML::Key << "amp";
			yamlOutput << YAML::Value << component->getChannel(idxDecayChannel).getWaveName();

			if(component->mapCouplingToMasterChannel(component->mapChannelToCoupling(idxDecayChannel)) == idxDecayChannel) {
				yamlOutput << YAML::Key << "couplings";
				yamlOutput << YAML::Value;
				yamlOutput << YAML::BeginSeq;

				const std::vector<size_t>& bins = component->getChannel(idxDecayChannel).getBins();
				for(size_t i = 0; i < bins.size(); ++i) {
					const size_t idxBin = bins[i];
					yamlOutput << YAML::Flow;
					yamlOutput << YAML::BeginSeq;

					yamlOutput << fitParameters.getCoupling(component->getId(), component->mapChannelToCoupling(idxDecayChannel), idxBin).real();
					yamlOutput << fitParameters.getCoupling(component->getId(), component->mapChannelToCoupling(idxDecayChannel), idxBin).imag();

					yamlOutput << YAML::EndSeq;
					yamlOutput << YAML::Block;
				}

				yamlOutput << YAML::EndSeq;
			}

			if(component->getNrBranchings() > 1) {
				if(component->mapBranchingToMasterChannel(component->mapChannelToBranching(idxDecayChannel)) == idxDecayChannel) {
					yamlOutput << YAML::Key << "branching";
					yamlOutput << YAML::Value;

					yamlOutput << YAML::Flow;
					yamlOutput << YAML::BeginSeq;

					yamlOutput << fitParameters.getBranching(component->getId(), component->mapChannelToBranching(idxDecayChannel)).real();
					yamlOutput << fitParameters.getBranching(component->getId(), component->mapChannelToBranching(idxDecayChannel)).imag();

					yamlOutput << YAML::EndSeq;
					yamlOutput << YAML::Block;
				}
			}

			if(std::dynamic_pointer_cast<const rpwa::resonanceFit::dynamicWidthBreitWigner>(component)) {
				writeModelComponentDecayChannelMIsobar1(yamlOutput,
				                                        std::dynamic_pointer_cast<const rpwa::resonanceFit::dynamicWidthBreitWigner>(component),
				                                        idxDecayChannel);
				writeModelComponentDecayChannelMIsobar2(yamlOutput,
				                                        std::dynamic_pointer_cast<const rpwa::resonanceFit::dynamicWidthBreitWigner>(component),
				                                        idxDecayChannel);
				writeModelComponentDecayChannelRelAngularMom(yamlOutput,
				                                             std::dynamic_pointer_cast<const rpwa::resonanceFit::dynamicWidthBreitWigner>(component),
				                                             idxDecayChannel);
				writeModelComponentDecayChannelBranchingRatio(yamlOutput,
				                                              std::dynamic_pointer_cast<const rpwa::resonanceFit::dynamicWidthBreitWigner>(component),
				                                              idxDecayChannel);
			}
			if(std::dynamic_pointer_cast<const rpwa::resonanceFit::integralWidthBreitWigner>(component)) {
				writeModelComponentDecayChannelIntegral(yamlOutput,
				                                        std::dynamic_pointer_cast<const rpwa::resonanceFit::integralWidthBreitWigner>(component),
				                                        idxDecayChannel);
				writeModelComponentDecayChannelBranchingRatio(yamlOutput,
				                                              std::dynamic_pointer_cast<const rpwa::resonanceFit::integralWidthBreitWigner>(component),
				                                              idxDecayChannel);
			}

			yamlOutput << YAML::EndMap;
		}

		yamlOutput << YAML::EndSeq;

		if(component->getTotalNrChannels() > component->getNrChannels()) {
			yamlOutput << YAML::Key << "extradecaychannels";
			yamlOutput << YAML::Value;
			yamlOutput << YAML::BeginSeq;

			for(size_t idxDecayChannel = component->getNrChannels(); idxDecayChannel < component->getTotalNrChannels(); ++idxDecayChannel) {
				yamlOutput << YAML::BeginMap;

				if(std::dynamic_pointer_cast<const rpwa::resonanceFit::dynamicWidthBreitWigner>(component)) {
					writeModelComponentDecayChannelMIsobar1(yamlOutput,
					                                        std::dynamic_pointer_cast<const rpwa::resonanceFit::dynamicWidthBreitWigner>(component),
					                                        idxDecayChannel);
					writeModelComponentDecayChannelMIsobar2(yamlOutput,
					                                        std::dynamic_pointer_cast<const rpwa::resonanceFit::dynamicWidthBreitWigner>(component),
					                                        idxDecayChannel);
					writeModelComponentDecayChannelRelAngularMom(yamlOutput,
					                                             std::dynamic_pointer_cast<const rpwa::resonanceFit::dynamicWidthBreitWigner>(component),
					                                             idxDecayChannel);
					writeModelComponentDecayChannelBranchingRatio(yamlOutput,
					                                              std::dynamic_pointer_cast<const rpwa::resonanceFit::dynamicWidthBreitWigner>(component),
					                                              idxDecayChannel);
				}
				if(std::dynamic_pointer_cast<const rpwa::resonanceFit::integralWidthBreitWigner>(component)) {
					writeModelComponentDecayChannelIntegral(yamlOutput,
					                                        std::dynamic_pointer_cast<const rpwa::resonanceFit::integralWidthBreitWigner>(component),
					                                        idxDecayChannel);
					writeModelComponentDecayChannelBranchingRatio(yamlOutput,
					                                              std::dynamic_pointer_cast<const rpwa::resonanceFit::integralWidthBreitWigner>(component),
					                                              idxDecayChannel);
				}

				yamlOutput << YAML::EndMap;
			}

			yamlOutput << YAML::EndSeq;
		}

		if(std::dynamic_pointer_cast<const rpwa::resonanceFit::exponentialBackground>(component)) {
			writeModelComponentDecayChannelMIsobar1(yamlOutput,
			                                        std::dynamic_pointer_cast<const rpwa::resonanceFit::exponentialBackground>(component));
			writeModelComponentDecayChannelMIsobar2(yamlOutput,
			                                        std::dynamic_pointer_cast<const rpwa::resonanceFit::exponentialBackground>(component));
			writeModelComponentDecayChannelRelAngularMom(yamlOutput,
			                                             std::dynamic_pointer_cast<const rpwa::resonanceFit::exponentialBackground>(component));
			writeModelComponentDecayChannelExponent(yamlOutput,
			                                        std::dynamic_pointer_cast<const rpwa::resonanceFit::exponentialBackground>(component));
		}
		if(std::dynamic_pointer_cast<const rpwa::resonanceFit::tPrimeDependentBackground>(component)) {
			writeModelComponentDecayChannelMIsobar1(yamlOutput,
			                                        std::dynamic_pointer_cast<const rpwa::resonanceFit::tPrimeDependentBackground>(component));
			writeModelComponentDecayChannelMIsobar2(yamlOutput,
			                                        std::dynamic_pointer_cast<const rpwa::resonanceFit::tPrimeDependentBackground>(component));
			writeModelComponentDecayChannelRelAngularMom(yamlOutput,
			                                             std::dynamic_pointer_cast<const rpwa::resonanceFit::tPrimeDependentBackground>(component));
			writeModelComponentDecayChannelExponent(yamlOutput,
			                                        std::dynamic_pointer_cast<const rpwa::resonanceFit::tPrimeDependentBackground>(component));
		}
		if(std::dynamic_pointer_cast<const rpwa::resonanceFit::exponentialBackgroundIntegral>(component)) {
			writeModelComponentDecayChannelIntegral(yamlOutput,
			                                        std::dynamic_pointer_cast<const rpwa::resonanceFit::exponentialBackgroundIntegral>(component));
			writeModelComponentDecayChannelExponent(yamlOutput,
			                                        std::dynamic_pointer_cast<const rpwa::resonanceFit::exponentialBackgroundIntegral>(component));
		}
		if(std::dynamic_pointer_cast<const rpwa::resonanceFit::tPrimeDependentBackgroundIntegral>(component)) {
			writeModelComponentDecayChannelIntegral(yamlOutput,
			                                        std::dynamic_pointer_cast<const rpwa::resonanceFit::tPrimeDependentBackgroundIntegral>(component));
			writeModelComponentDecayChannelExponent(yamlOutput,
			                                        std::dynamic_pointer_cast<const rpwa::resonanceFit::tPrimeDependentBackgroundIntegral>(component));
		}

		yamlOutput << YAML::EndMap;
	}


	std::vector<rpwa::resonanceFit::componentPtr>
	readModelComponents(const YAML::Node& configModel,
	                    const rpwa::resonanceFit::informationConstPtr& fitInformation,
	                    rpwa::resonanceFit::parameters& fitParameters,
	                    rpwa::resonanceFit::parameters& fitParametersError,
	                    const boost::multi_array<std::string, 2>& waveNames,
	                    const std::vector<size_t>& nrMassBins,
	                    const boost::multi_array<double, 2>& massBinCenters,
	                    const boost::multi_array<double, 3>& phaseSpaceIntegrals,
	                    const bool useBranchings)
	{
		if(debug) {
			printDebug << "reading 'components'." << std::endl;
		}

		const YAML::Node& configComponents = configModel["components"];
		if(not configComponents) {
			printErr << "'components' is not a valid YAML node." << std::endl;
			throw;
		}
		if(not configComponents.IsSequence()) {
			printErr << "'components' is not a YAML sequence." << std::endl;
			throw;
		}

		const size_t nrComponents = configComponents.size();
		if(debug) {
			printDebug << "reading " << nrComponents << " components from configuration file." << std::endl;
		}

		std::vector<rpwa::resonanceFit::componentPtr> components;
		for(size_t idxComponent = 0; idxComponent < nrComponents; ++idxComponent) {
			const YAML::Node& configComponent = configComponents[idxComponent];

			const rpwa::resonanceFit::componentPtr& component = readModelComponent(configComponent,
			                                                                       fitInformation,
			                                                                       fitParameters,
			                                                                       fitParametersError,
			                                                                       components.size(),
			                                                                       waveNames,
			                                                                       nrMassBins,
			                                                                       massBinCenters,
			                                                                       phaseSpaceIntegrals,
			                                                                       useBranchings);

			for(size_t idx = 0; idx < components.size(); ++idx) {
				if(components[idx]->getName() == component->getName()) {
					printErr << "component '" << component->getName() << "' defined twice." << std::endl;
					throw;
				}
			}

			components.push_back(component);
		}

		return components;
	}


	void
	writeModelComponents(YAML::Emitter& yamlOutput,
	                     const std::vector<rpwa::resonanceFit::componentConstPtr>& components,
	                     const rpwa::resonanceFit::parameters& fitParameters,
	                     const rpwa::resonanceFit::parameters& fitParametersError)
	{
		if(debug) {
			printDebug << "writing 'components'." << std::endl;
		}

		yamlOutput << YAML::Key << "components";
		yamlOutput << YAML::Value;

		yamlOutput << YAML::BeginSeq;

		for(size_t idxComponent = 0; idxComponent < components.size(); ++idxComponent) {
			writeModelComponent(yamlOutput,
			                    components[idxComponent],
			                    fitParameters,
			                    fitParametersError);
		}

		yamlOutput << YAML::EndSeq;
	}


	void
	readModelFsmdBin(const YAML::Node& configFsmd,
	                 std::shared_ptr<TFormula>& function,
	                 boost::multi_array<rpwa::resonanceFit::parameter, 1>& parameters)
	{
		if(debug) {
			printDebug << "reading 'finalStateMassDependence' for an individual bin." << std::endl;
		}

		if(not configFsmd.IsMap()) {
			printErr << "'finalStateMassDependence' for an individual bin is not a YAML map." << std::endl;
			throw;
		}

		std::map<std::string, rpwa::YamlCppUtils::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("formula", rpwa::YamlCppUtils::TypeString);
		if(not checkIfAllVariablesAreThere(configFsmd, mandatoryArguments)) {
			printErr << "'finalStateMassDependence' for an individual bin does not contain all required variables." << std::endl;
			return throw;
		}

		const std::string formula = configFsmd["formula"].as<std::string>();
		function.reset(new TFormula("finalStateMassDependence", formula.c_str()));

		const size_t nrParameters = function->GetNpar();
		parameters.resize(boost::extents[nrParameters]);

		for(size_t idxParameter = 0; idxParameter < nrParameters; ++idxParameter) {
			const std::string parameterName = function->GetParName(idxParameter);

			// set default value for step size
			parameters[idxParameter].setStep(0.0001);

			readModelParameter(configFsmd,
			                   parameterName,
			                   parameters[idxParameter]);
		}
	}


	void
	writeModelFsmdBin(YAML::Emitter& yamlOutput,
	                  const rpwa::resonanceFit::fsmdConstPtr& fsmd,
	                  const size_t idxBin,
	                  const rpwa::resonanceFit::parameters& fitParameters,
	                  const rpwa::resonanceFit::parameters& fitParametersError)
	{
		if(debug) {
			printDebug << "writing 'finalStateMassDependence' for an individual bin." << std::endl;
		}

		yamlOutput << YAML::BeginMap;

		yamlOutput << YAML::Key << "formula";
		yamlOutput << YAML::Value << fsmd->getFunction(idxBin)->GetTitle();

		for(size_t idxParameter = 0; idxParameter < fsmd->getNrParameters(idxBin); ++idxParameter) {
			writeModelParameter(yamlOutput,
			                    fsmd->getParameter(idxBin, idxParameter),
			                    fitParameters.getParameter(fsmd->getId(), fsmd->getParameterIndex(idxBin)+idxParameter),
			                    fitParametersError.getParameter(fsmd->getId(), fsmd->getParameterIndex(idxBin)+idxParameter));
		}

		yamlOutput << YAML::EndMap;
	}


	rpwa::resonanceFit::fsmdPtr
	readModelFsmd(const YAML::Node& configModel,
	              rpwa::resonanceFit::parameters& fitParameters,
	              rpwa::resonanceFit::parameters& fitParametersError,
	              const size_t id,
	              const std::vector<size_t>& nrMassBins,
	              const boost::multi_array<double, 2>& massBinCenters)
	{
		if(debug) {
			printDebug << "reading 'finalStateMassDependence'." << std::endl;
		}

		const YAML::Node& configFsmd = configModel["finalStateMassDependence"];
		if(not configFsmd) {
			// final-state mass-dependence might not be specified
			return rpwa::resonanceFit::fsmdPtr();
		}

		if(not configFsmd.IsMap() and not configFsmd.IsSequence()) {
			printErr << "'finalStateMassDependence' is not a YAML map or sequence." << std::endl;
			throw;
		}

		rpwa::resonanceFit::fsmdPtr fsmd;
		if(configFsmd.IsMap()) {
			// a single final-state mass-dependence is given
			std::shared_ptr<TFormula> function;
			boost::multi_array<rpwa::resonanceFit::parameter, 1> parameters;

			readModelFsmdBin(configFsmd,
			                 function,
			                 parameters);

			const size_t nrParameters = function->GetNpar();
			fitParameters.resize(id+1, 0, nrParameters, 0);
			fitParametersError.resize(id+1, 0, nrParameters, 0);

			for(size_t idxParameter = 0; idxParameter < nrParameters; ++idxParameter) {
				fitParameters.setParameter(id, idxParameter, parameters[idxParameter].startValue());
				fitParametersError.setParameter(id, idxParameter, parameters[idxParameter].startError());
			}

			fsmd.reset(new rpwa::resonanceFit::fsmd(id,
			                                        nrMassBins,
			                                        massBinCenters,
			                                        function,
			                                        parameters));
		} else {
			// a final-state mass-dependence for each bin is given
			std::vector<std::shared_ptr<TFormula> > functions;
			boost::multi_array<rpwa::resonanceFit::parameter, 2> parameters;

			size_t nrParametersSum = 0;
			const size_t nrBins = configFsmd.size();
			for(size_t idxBin = 0; idxBin < nrBins; ++idxBin) {
				std::shared_ptr<TFormula> tempFunction;
				boost::multi_array<rpwa::resonanceFit::parameter, 1> tempParameters;

				readModelFsmdBin(configFsmd[idxBin],
				                 tempFunction,
				                 tempParameters);

				rpwa::resonanceFit::adjustSizeAndSet(functions, idxBin, tempFunction);
				rpwa::resonanceFit::adjustSizeAndSet(parameters, idxBin, tempParameters);

				const size_t idxParameterFirst = nrParametersSum;
				const size_t nrParameters = tempFunction->GetNpar();
				nrParametersSum += nrParameters;

				fitParameters.resize(id+1, 0, nrParametersSum, nrBins);
				fitParametersError.resize(id+1, 0, nrParametersSum, nrBins);

				for(size_t idxParameter = 0; idxParameter < nrParameters; ++idxParameter) {
					fitParameters.setParameter(id, idxParameterFirst+idxParameter, tempParameters[idxParameter].startValue());
					fitParametersError.setParameter(id, idxParameterFirst+idxParameter, tempParameters[idxParameter].startError());
				}
			}

			fsmd.reset(new rpwa::resonanceFit::fsmd(id,
			                                        nrMassBins,
			                                        massBinCenters,
			                                        functions,
			                                        parameters));
		}

		if(debug) {
			printDebug << *fsmd;
		}

		return fsmd;
	}


	void
	writeModelFsmd(YAML::Emitter& yamlOutput,
	               const rpwa::resonanceFit::fsmdConstPtr& fsmd,
	               const rpwa::resonanceFit::parameters& fitParameters,
	               const rpwa::resonanceFit::parameters& fitParametersError)
	{
		if(debug) {
			printDebug << "writing 'finalStateMassDependence'." << std::endl;
		}

		yamlOutput << YAML::Key << "finalStateMassDependence";
		yamlOutput << YAML::Value;

		if(not fsmd->isSameFunctionForAllBins()) {
			yamlOutput << YAML::BeginSeq;
		}

		const size_t maxNrBins = fsmd->isSameFunctionForAllBins() ? 1 : fsmd->getNrBins();
		for(size_t idxBin = 0; idxBin < maxNrBins; ++idxBin) {
			writeModelFsmdBin(yamlOutput,
			                  fsmd,
			                  idxBin,
			                  fitParameters,
			                  fitParametersError);
		}

		if(not fsmd->isSameFunctionForAllBins()) {
			yamlOutput << YAML::EndSeq;
		}
	}


	rpwa::resonanceFit::modelConstPtr
	readModel(const YAML::Node& configRoot,
	          const rpwa::resonanceFit::informationConstPtr& fitInformation,
	          rpwa::resonanceFit::parameters& fitParameters,
	          rpwa::resonanceFit::parameters& fitParametersError,
	          const boost::multi_array<std::string, 2>& waveNames,
	          const std::vector<size_t>& nrMassBins,
	          const boost::multi_array<double, 2>& massBinCenters,
	          const boost::multi_array<double, 3>& phaseSpaceIntegrals,
	          const bool useBranchings)
	{
		if(debug) {
			printDebug << "reading 'model'." << std::endl;
		}

		const YAML::Node& configModel = configRoot["model"];

		if(not configModel) {
			printErr << "'model' does not exist in configuration file." << std::endl;
			throw;
		}

		std::string anchorWaveName;
		std::string anchorComponentName;
		readModelAnchors(configModel,
		                 anchorWaveName,
		                 anchorComponentName);

		const std::vector<rpwa::resonanceFit::componentPtr> components = readModelComponents(configModel,
		                                                                                     fitInformation,
		                                                                                     fitParameters,
		                                                                                     fitParametersError,
		                                                                                     waveNames,
		                                                                                     nrMassBins,
		                                                                                     massBinCenters,
		                                                                                     phaseSpaceIntegrals,
		                                                                                     useBranchings);

		// get information for creating the final-state mass-dependence
		const rpwa::resonanceFit::fsmdPtr& fsmd = readModelFsmd(configModel,
		                                                        fitParameters,
		                                                        fitParametersError,
		                                                        components.size(),
		                                                        nrMassBins,
		                                                        massBinCenters);

		return std::make_shared<rpwa::resonanceFit::model>(fitInformation,
		                                                   components,
		                                                   fsmd,
		                                                   anchorWaveName,
		                                                   anchorComponentName);
	}


	void
	writeModel(YAML::Emitter& yamlOutput,
	           const rpwa::resonanceFit::modelConstPtr& fitModel,
	           const rpwa::resonanceFit::parameters& fitParameters,
	           const rpwa::resonanceFit::parameters& fitParametersError)
	{
		yamlOutput << YAML::Key << "model";
		yamlOutput << YAML::Value;

		yamlOutput << YAML::BeginMap;

		writeModelAnchors(yamlOutput,
		                  fitModel->getAnchorWaveName(),
		                  fitModel->getAnchorComponentName());

		std::vector<rpwa::resonanceFit::componentConstPtr> components;
		for(size_t idxComponent = 0; idxComponent < fitModel->getNrComponents(); ++idxComponent) {
			components.push_back(fitModel->getComponent(idxComponent));
		}
		writeModelComponents(yamlOutput,
		                     components,
		                     fitParameters,
		                     fitParametersError);

		writeModelFsmd(yamlOutput,
		               fitModel->getFsmd(),
		               fitParameters,
		               fitParametersError);

		yamlOutput << YAML::EndMap;
	}


	void
	prepareMassLimit(const rpwa::resonanceFit::informationConstPtr& fitInformation,
	                 const size_t nrMassBins,
	                 const boost::multi_array<double, 1>& massBinCenters,
	                 boost::multi_array<std::pair<size_t, size_t>, 2>& wavePairMassBinLimits)
	{
		if(debug) {
			printDebug << "determine which mass bins to use in the fit for " << nrMassBins << " mass bins, center of first and last mass bins: "
			           << massBinCenters[0] << " and " << massBinCenters[nrMassBins - 1] << " GeV/c^2." << std::endl;
		}

		// determine which mass bins to use for a single wave
		std::vector<std::pair<size_t, size_t> > waveMassBinLimits(fitInformation->nrWaves());
		for(size_t idxWave = 0; idxWave < fitInformation->nrWaves(); ++idxWave) {
			const rpwa::resonanceFit::information::wave& wave = fitInformation->getWave(idxWave);

			size_t binFirst = 0;
			size_t binLast = nrMassBins-1;
			for(size_t idxMass = 0; idxMass < nrMassBins; ++idxMass) {
				if(massBinCenters[idxMass] < wave.massLimits().first) {
					binFirst = idxMass+1;
				}
				if(massBinCenters[idxMass] == wave.massLimits().first) {
					binFirst = idxMass;
				}
				if(massBinCenters[idxMass] <= wave.massLimits().second) {
					binLast = idxMass;
				}
			}
			if(wave.massLimits().first < 0) {
				binFirst = 0;
			}
			if(wave.massLimits().second < 0) {
				binLast = nrMassBins-1;
			}

			const double massStep = (massBinCenters[nrMassBins-1] - massBinCenters[0]) / (nrMassBins - 1);
			const double massFirst = massBinCenters[binFirst] - massStep/2.;
			const double massLast = massBinCenters[binLast] + massStep/2.;
			if(debug) {
				printDebug << idxWave << ": " << wave.waveName() << ": "
				           << "mass range: " << massFirst << "-" << massLast << " GeV/c^2, "
				           << "bin range: " << binFirst << "-" << binLast << std::endl;
			}
			waveMassBinLimits[idxWave] = std::make_pair(binFirst, binLast);
		}

		// determine which mass bins to use for each pair of waves
		wavePairMassBinLimits.resize(boost::extents[fitInformation->nrWaves()][fitInformation->nrWaves()]);
		for(size_t idxWave = 0; idxWave < fitInformation->nrWaves(); ++idxWave) {
			for(size_t jdxWave = 0; jdxWave < fitInformation->nrWaves(); ++jdxWave) {
				wavePairMassBinLimits[idxWave][jdxWave] = std::make_pair(std::max(waveMassBinLimits[idxWave].first,  waveMassBinLimits[jdxWave].first),
				                                                         std::min(waveMassBinLimits[idxWave].second, waveMassBinLimits[jdxWave].second));
			}
		}

		if(debug) {
			printDebug << "waves and mass limits:" << std::endl;
			for(size_t idxWave = 0; idxWave < fitInformation->nrWaves(); ++idxWave) {
				const rpwa::resonanceFit::information::wave& wave = fitInformation->getWave(idxWave);

				std::ostringstream output;
				for(size_t jdxWave = 0; jdxWave < fitInformation->nrWaves(); ++jdxWave) {
					output << wavePairMassBinLimits[idxWave][jdxWave].first << "-" << wavePairMassBinLimits[idxWave][jdxWave].second << " ";
				}
				printDebug << wave.waveName() << " " << waveMassBinLimits[idxWave].first << "-" << waveMassBinLimits[idxWave].second
				           << ": " << output.str() << std::endl;
			}
		}
	}


	void
	prepareMassLimits(const rpwa::resonanceFit::informationConstPtr& fitInformation,
	                  const std::vector<size_t>& nrMassBins,
	                  const boost::multi_array<double, 2>& massBinCenters,
	                  boost::multi_array<std::pair<size_t, size_t>, 3>& wavePairMassBinLimits)
	{
		if(debug) {
			printDebug << "determine which mass bins to use in the fit." << std::endl;
		}

		for(size_t idxBin = 0; idxBin < nrMassBins.size(); ++idxBin) {
			boost::multi_array<std::pair<size_t, size_t>, 2> tempWavePairMassBinLimits;

			prepareMassLimit(fitInformation,
			                 nrMassBins[idxBin],
			                 massBinCenters[idxBin],
			                 tempWavePairMassBinLimits);

			rpwa::resonanceFit::adjustSizeAndSet(wavePairMassBinLimits, idxBin, tempWavePairMassBinLimits);
		}
	}


	void
	createPlotsWave(const rpwa::resonanceFit::informationConstPtr& fitInformation,
	                const rpwa::resonanceFit::dataConstPtr& fitData,
	                const rpwa::resonanceFit::modelConstPtr& fitModel,
	                const rpwa::resonanceFit::parameters& fitParameters,
	                rpwa::resonanceFit::cache& cache,
	                TDirectory* outDirectory,
	                const bool rangePlotting,
	                const size_t extraBinning,
	                const size_t idxWave,
	                const size_t idxBin)
	{
		const rpwa::resonanceFit::information::wave& wave = fitInformation->getWave(idxWave);
		if(debug) {
			printDebug << "start creating plots for wave '" << wave.waveName() << "' in bin " << idxBin << "." << std::endl;
		}

		TMultiGraph graphs;
		graphs.SetName(wave.waveName().c_str());
		graphs.SetTitle(wave.waveName().c_str());

		TGraphErrors* systematics = NULL;
		if(fitInformation->getBin(idxBin).sysFileNames().size() > 0) {
			systematics = new TGraphErrors;
			systematics->SetName((wave.waveName() + "__sys").c_str());
			systematics->SetTitle((wave.waveName() + "__sys").c_str());
			systematics->SetLineColor(kAzure-9);
			systematics->SetFillColor(kAzure-9);
			graphs.Add(systematics, "2");
		}

		TGraphErrors* data = new TGraphErrors;
		data->SetName((wave.waveName() + "__data").c_str());
		data->SetTitle((wave.waveName() + "__data").c_str());
		graphs.Add(data, "P");

		TGraph* fit = new TGraph;
		fit->SetName((wave.waveName() + "__fit").c_str());
		fit->SetTitle((wave.waveName() + "__fit").c_str());
		fit->SetLineColor(kRed);
		fit->SetLineWidth(2);
		fit->SetMarkerColor(kRed);
		graphs.Add(fit, "L");

		TGraph* phaseSpace = new TGraph;
		phaseSpace->SetName((wave.waveName() + "__ps").c_str());
		phaseSpace->SetTitle((wave.waveName() + "__ps").c_str());
		graphs.Add(phaseSpace, "L");

		const std::vector<std::pair<size_t, size_t> >& compChannel = fitModel->getComponentChannel(idxBin, idxWave);
		std::vector<TGraph*> components;
		for(size_t idxComponents = 0; idxComponents < compChannel.size(); ++idxComponents) {
			const size_t idxComponent = compChannel[idxComponents].first;
			TGraph* component = new TGraph;
			component->SetName((wave.waveName() + "__" + fitModel->getComponent(idxComponent)->getName()).c_str());
			component->SetTitle((wave.waveName() + "__" + fitModel->getComponent(idxComponent)->getName()).c_str());

			Color_t color = kBlue;
			if(fitModel->getComponent(idxComponent)->getName().find("bkg") != std::string::npos) {
				color = kMagenta;
			}
			component->SetLineColor(color);
			component->SetMarkerColor(color);

			graphs.Add(component, "L");
			components.push_back(component);
		}

		// plot data
		double maxIE = -std::numeric_limits<double>::max();
		for(size_t point = 0; point <= (fitData->nrMassBins()[idxBin]-1); ++point) {
			const size_t idxMass = point;
			const double mass = fitData->massBinCenters()[idxBin][idxMass];
			const double halfBin = 0.5 * (fitData->massBinCenters()[idxBin][fitData->nrMassBins()[idxBin]-1] - fitData->massBinCenters()[idxBin][0]) / (fitData->nrMassBins()[idxBin] - 1);

			data->SetPoint(point, mass, fitData->plottingIntensities()[idxBin][idxMass][idxWave].first);
			data->SetPointError(point, halfBin, fitData->plottingIntensities()[idxBin][idxMass][idxWave].second);
			maxIE = std::max(maxIE, fitData->plottingIntensities()[idxBin][idxMass][idxWave].first+fitData->plottingIntensities()[idxBin][idxMass][idxWave].second);

			if(fitInformation->getBin(idxBin).sysFileNames().size() > 0) {
				const double minSI = fitData->sysPlottingIntensities()[idxBin][idxMass][idxWave].first;
				const double maxSI = fitData->sysPlottingIntensities()[idxBin][idxMass][idxWave].second;
				systematics->SetPoint(point, mass, (maxSI+minSI)/2.);
				systematics->SetPointError(point, halfBin, (maxSI-minSI)/2.);
				maxIE = std::max(maxIE, maxSI);
			}
		}

		// plot fit, either over full or limited mass range
		const size_t firstPoint = rangePlotting ? (extraBinning*fitData->wavePairMassBinLimits()[idxBin][idxWave][idxWave].first) : 0;
		const size_t lastPoint = rangePlotting ? (extraBinning*fitData->wavePairMassBinLimits()[idxBin][idxWave][idxWave].second) : (extraBinning*(fitData->nrMassBins()[idxBin]-1));
		for(size_t point = firstPoint; point <= lastPoint; ++point) {
			const size_t idxMass = (point%extraBinning == 0) ? (point/extraBinning) : std::numeric_limits<size_t>::max();
			const double massStep = (fitData->massBinCenters()[idxBin][fitData->nrMassBins()[idxBin]-1] - fitData->massBinCenters()[idxBin][0]) / (fitData->nrMassBins()[idxBin] - 1) / extraBinning;
			const double mass = (idxMass != std::numeric_limits<size_t>::max()) ? fitData->massBinCenters()[idxBin][idxMass] : (fitData->massBinCenters()[idxBin][point/extraBinning] + (point%extraBinning) * massStep);

			const double intensity = fitModel->intensity(fitParameters, cache, idxWave, idxBin, mass, idxMass);
			fit->SetPoint(point-firstPoint, mass, intensity);
			maxIE = std::max(maxIE, intensity);

			for(size_t idxComponents = 0; idxComponents < compChannel.size(); ++idxComponents) {
				const size_t idxComponent = compChannel[idxComponents].first;
				const size_t idxChannel = compChannel[idxComponents].second;

				std::complex<double> prodAmp = fitModel->getComponent(idxComponent)->val(fitParameters, cache, idxBin, mass, idxMass);
				prodAmp *= fitModel->getComponent(idxComponent)->getCouplingPhaseSpace(fitParameters, cache, idxChannel, idxBin, mass, idxMass);
				if(fitModel->getFsmd()) {
					prodAmp *= fitModel->getFsmd()->val(fitParameters, cache, idxBin, mass, idxMass);
				}

				components[idxComponents]->SetPoint(point-firstPoint, mass, norm(prodAmp));
				maxIE = std::max(maxIE, norm(prodAmp));
			}
		}

		boost::multi_array<double, 2>::const_array_view<1>::type viewM = fitData->massBinCenters()[boost::indices[idxBin][boost::multi_array<double, 2>::index_range(0, fitData->nrMassBins()[idxBin])]];
		boost::multi_array<double, 3>::const_array_view<1>::type viewInt = fitData->phaseSpaceIntegrals()[boost::indices[idxBin][boost::multi_array<double, 3>::index_range(0, fitData->nrMassBins()[idxBin])][idxWave]];
		ROOT::Math::Interpolator phaseSpaceInterpolator(std::vector<double>(viewM.begin(), viewM.end()), std::vector<double>(viewInt.begin(), viewInt.end()), ROOT::Math::Interpolation::kLINEAR);

		// plot phase-space
		double maxP = -std::numeric_limits<double>::max();
		for(size_t point = 0; point <= (extraBinning*(fitData->nrMassBins()[idxBin]-1)); ++point) {
			const size_t idxMass = (point%extraBinning == 0) ? (point/extraBinning) : std::numeric_limits<size_t>::max();
			const double massStep = (fitData->massBinCenters()[idxBin][fitData->nrMassBins()[idxBin]-1] - fitData->massBinCenters()[idxBin][0]) / (fitData->nrMassBins()[idxBin] - 1) / extraBinning;
			const double mass = (idxMass != std::numeric_limits<size_t>::max()) ? fitData->massBinCenters()[idxBin][idxMass] : (fitData->massBinCenters()[idxBin][point/extraBinning] + (point%extraBinning) * massStep);

			double ps = pow((idxMass != std::numeric_limits<size_t>::max()) ? fitData->phaseSpaceIntegrals()[idxBin][idxMass][idxWave] : phaseSpaceInterpolator.Eval(mass), 2);
			if(fitModel->getFsmd()) {
				ps *= std::norm(fitModel->getFsmd()->val(fitParameters, cache, idxBin, mass, idxMass));
			}
			phaseSpace->SetPoint(point, mass, ps);
			maxP = std::max(maxP, ps);
		}

		// scale phase-space graph to half-height of intensity graphs
		for(Int_t idx = 0; idx < phaseSpace->GetN(); ++idx) {
			double x, y;
			phaseSpace->GetPoint(idx, x, y);
			phaseSpace->SetPoint(idx, x, y * 0.5 * maxIE/maxP);
		}

		outDirectory->cd();
		graphs.Write();
	}


	void
	createPlotsWaveSum(const rpwa::resonanceFit::informationConstPtr& fitInformation,
	                   const rpwa::resonanceFit::dataConstPtr& fitData,
	                   const rpwa::resonanceFit::modelConstPtr& fitModel,
	                   const rpwa::resonanceFit::parameters& fitParameters,
	                   rpwa::resonanceFit::cache& cache,
	                   TDirectory* outDirectory,
	                   const bool rangePlotting,
	                   const size_t extraBinning,
	                   const size_t idxWave)
	{
		const rpwa::resonanceFit::information::wave& wave = fitInformation->getWave(idxWave);
		if(debug) {
			printDebug << "start creating plots for wave '" << wave.waveName() << "' for sum over all bins." << std::endl;
		}

		// all mass binnings must be the same to be able to create the sum plots
		if(not fitData->hasSameMassBinning() or not fitModel->isMappingEqualInAllBins()) {
			printErr << "cannot create plots for wave '" << wave.waveName() << "' for sum over all bins if the bins used different mass binnings." << std::endl;
			throw;
		}
		const size_t idxBin = 0;

		TMultiGraph graphs;
		graphs.SetName(wave.waveName().c_str());
		graphs.SetTitle(wave.waveName().c_str());

		TGraphErrors* data = new TGraphErrors;
		data->SetName((wave.waveName() + "__data").c_str());
		data->SetTitle((wave.waveName() + "__data").c_str());
		graphs.Add(data, "P");

		TGraph* fit = new TGraph;
		fit->SetName((wave.waveName() + "__fit").c_str());
		fit->SetTitle((wave.waveName() + "__fit").c_str());
		fit->SetLineColor(kRed);
		fit->SetLineWidth(2);
		fit->SetMarkerColor(kRed);
		graphs.Add(fit, "L");

		const std::vector<std::pair<size_t, size_t> >& compChannel = fitModel->getComponentChannel(idxBin, idxWave);
		std::vector<TGraph*> components;
		for(size_t idxComponents = 0; idxComponents < compChannel.size(); ++idxComponents) {
			const size_t idxComponent = compChannel[idxComponents].first;
			TGraph* component = new TGraph;
			component->SetName((wave.waveName() + "__" + fitModel->getComponent(idxComponent)->getName()).c_str());
			component->SetTitle((wave.waveName() + "__" + fitModel->getComponent(idxComponent)->getName()).c_str());

			Color_t color = kBlue;
			if(fitModel->getComponent(idxComponent)->getName().find("bkg") != std::string::npos) {
				color = kMagenta;
			}
			component->SetLineColor(color);
			component->SetMarkerColor(color);

			graphs.Add(component, "L");
			components.push_back(component);
		}

		// plot data
		for(size_t point = 0; point <= (fitData->nrMassBins()[idxBin]-1); ++point) {
			const size_t idxMass = point;
			const double mass = fitData->massBinCenters()[idxBin][idxMass];
			const double halfBin = 0.5 * (fitData->massBinCenters()[idxBin][fitData->nrMassBins()[idxBin]-1] - fitData->massBinCenters()[idxBin][0]) / (fitData->nrMassBins()[idxBin] - 1);

			double sum = 0.;
			double error2 = 0.;
			for(size_t idxBin = 0; idxBin < fitInformation->nrBins(); ++idxBin) {
				sum += fitData->plottingIntensities()[idxBin][idxMass][idxWave].first;
				error2 += std::pow(fitData->plottingIntensities()[idxBin][idxMass][idxWave].second, 2);
			}
			data->SetPoint(point, mass, sum);
			data->SetPointError(point, halfBin, sqrt(error2));
		}

		// plot fit, either over full or limited mass range
		const size_t firstPoint = rangePlotting ? (extraBinning*fitData->wavePairMassBinLimits()[idxBin][idxWave][idxWave].first) : 0;
		const size_t lastPoint = rangePlotting ? (extraBinning*fitData->wavePairMassBinLimits()[idxBin][idxWave][idxWave].second) : (extraBinning*(fitData->nrMassBins()[idxBin]-1));
		for(size_t point = firstPoint; point <= lastPoint; ++point) {
			const size_t idxMass = (point%extraBinning == 0) ? (point/extraBinning) : std::numeric_limits<size_t>::max();
			const double massStep = (fitData->massBinCenters()[idxBin][fitData->nrMassBins()[idxBin]-1] - fitData->massBinCenters()[idxBin][0]) / (fitData->nrMassBins()[idxBin] - 1) / extraBinning;
			const double mass = (idxMass != std::numeric_limits<size_t>::max()) ? fitData->massBinCenters()[idxBin][idxMass] : (fitData->massBinCenters()[idxBin][point/extraBinning] + (point%extraBinning) * massStep);

			double sum = 0.;
			for(size_t idxBin = 0; idxBin < fitInformation->nrBins(); ++idxBin) {
				sum += fitModel->intensity(fitParameters, cache, idxWave, idxBin, mass, idxMass);
			}
			fit->SetPoint(point-firstPoint, mass, sum);

			for(size_t idxComponents = 0; idxComponents < compChannel.size(); ++idxComponents) {
				const size_t idxComponent = compChannel[idxComponents].first;
				const size_t idxChannel = compChannel[idxComponents].second;

				double sum = 0.;
				for(size_t idxBin = 0; idxBin < fitInformation->nrBins(); ++idxBin) {
					std::complex<double> prodAmp = fitModel->getComponent(idxComponent)->val(fitParameters, cache, idxBin, mass, idxMass);
					prodAmp *= fitModel->getComponent(idxComponent)->getCouplingPhaseSpace(fitParameters, cache, idxChannel, idxBin, mass, idxMass);
					if(fitModel->getFsmd()) {
						prodAmp *= fitModel->getFsmd()->val(fitParameters, cache, idxBin, mass, idxMass);
					}
					sum += norm(prodAmp);
				}
				components[idxComponents]->SetPoint(point-firstPoint, mass, sum);
			}
		}

		outDirectory->cd();
		graphs.Write();
	}


	void
	createPlotsWavePair(const rpwa::resonanceFit::informationConstPtr& fitInformation,
	                    const rpwa::resonanceFit::dataConstPtr& fitData,
	                    const rpwa::resonanceFit::modelConstPtr& fitModel,
	                    const rpwa::resonanceFit::parameters& fitParameters,
	                    rpwa::resonanceFit::cache& cache,
	                    TDirectory* outDirectory,
	                    const bool rangePlotting,
	                    const size_t extraBinning,
	                    const size_t idxWave,
	                    const size_t jdxWave,
	                    const size_t idxBin)
	{
		const rpwa::resonanceFit::information::wave& waveI = fitInformation->getWave(idxWave);
		const rpwa::resonanceFit::information::wave& waveJ = fitInformation->getWave(jdxWave);
		if(debug) {
			printDebug << "start creating plots for wave pair '" << waveI.waveName() << "' and '" << waveJ.waveName() << "' in bin " << idxBin << "." << std::endl;
		}

		const std::string realName = waveI.waveName() + "__" + waveJ.waveName() + "__real";
		const std::string imagName = waveI.waveName() + "__" + waveJ.waveName() + "__imag";
		const std::string phaseName = waveI.waveName() + "__" + waveJ.waveName() + "__phase";

		TMultiGraph real;
		real.SetName(realName.c_str());
		real.SetTitle(realName.c_str());

		TMultiGraph imag;
		imag.SetName(imagName.c_str());
		imag.SetTitle(imagName.c_str());

		TMultiGraph phase;
		phase.SetName(phaseName.c_str());
		phase.SetTitle(phaseName.c_str());

		TGraphErrors* realSystematics = NULL;
		TGraphErrors* imagSystematics = NULL;
		TGraphErrors* phaseSystematics = NULL;
		if(fitInformation->getBin(idxBin).sysFileNames().size() > 0) {
			realSystematics = new TGraphErrors;
			realSystematics->SetName((realName + "__sys").c_str());
			realSystematics->SetTitle((realName + "__sys").c_str());
			realSystematics->SetLineColor(kAzure-9);
			realSystematics->SetFillColor(kAzure-9);
			real.Add(realSystematics, "2");

			imagSystematics = new TGraphErrors;
			imagSystematics->SetName((imagName + "__sys").c_str());
			imagSystematics->SetTitle((imagName + "__sys").c_str());
			imagSystematics->SetLineColor(kAzure-9);
			imagSystematics->SetFillColor(kAzure-9);
			imag.Add(imagSystematics, "2");

			phaseSystematics = new TGraphErrors;
			phaseSystematics->SetName((phaseName + "__sys").c_str());
			phaseSystematics->SetTitle((phaseName + "__sys").c_str());
			phaseSystematics->SetLineColor(kAzure-9);
			phaseSystematics->SetFillColor(kAzure-9);
			phase.Add(phaseSystematics, "2");
		}

		TGraphErrors* realData = new TGraphErrors;
		realData->SetName((realName + "__data").c_str());
		realData->SetTitle((realName + "__data").c_str());
		real.Add(realData, "P");

		TGraphErrors* imagData = new TGraphErrors;
		imagData->SetName((imagName + "__data").c_str());
		imagData->SetTitle((imagName + "__data").c_str());
		imag.Add(imagData, "P");

		TGraphErrors* phaseData = new TGraphErrors;
		phaseData->SetName((phaseName + "__data").c_str());
		phaseData->SetTitle((phaseName + "__data").c_str());
		phase.Add(phaseData, "P");

		TGraph* realFit = new TGraph;
		realFit->SetName((realName + "__fit").c_str());
		realFit->SetTitle((realName + "__fit").c_str());
		realFit->SetLineColor(kRed);
		realFit->SetLineWidth(2);
		realFit->SetMarkerColor(kRed);
		real.Add(realFit, "L");

		TGraph* imagFit = new TGraph;
		imagFit->SetName((imagName + "__fit").c_str());
		imagFit->SetTitle((imagName + "__fit").c_str());
		imagFit->SetLineColor(kRed);
		imagFit->SetLineWidth(2);
		imagFit->SetMarkerColor(kRed);
		imag.Add(imagFit, "L");

		TGraph* phaseFit = new TGraph;
		phaseFit->SetName((phaseName + "__fit").c_str());
		phaseFit->SetTitle((phaseName + "__fit").c_str());
		phaseFit->SetLineColor(kRed);
		phaseFit->SetLineWidth(2);
		phaseFit->SetMarkerColor(kRed);
		phase.Add(phaseFit, "L");

		// plot data
		for(size_t point = 0; point <= (fitData->nrMassBins()[idxBin]-1); ++point) {
			const size_t idxMass = point;
			const double mass = fitData->massBinCenters()[idxBin][idxMass];
			const double halfBin = 0.5 * (fitData->massBinCenters()[idxBin][fitData->nrMassBins()[idxBin]-1] - fitData->massBinCenters()[idxBin][0]) / (fitData->nrMassBins()[idxBin] - 1);

			realData->SetPoint(point, mass, fitData->plottingSpinDensityMatrixElementsReal()[idxBin][idxMass][idxWave][jdxWave].first);
			realData->SetPointError(point, halfBin, fitData->plottingSpinDensityMatrixElementsReal()[idxBin][idxMass][idxWave][jdxWave].second);

			imagData->SetPoint(point, mass, fitData->plottingSpinDensityMatrixElementsImag()[idxBin][idxMass][idxWave][jdxWave].first);
			imagData->SetPointError(point, halfBin, fitData->plottingSpinDensityMatrixElementsImag()[idxBin][idxMass][idxWave][jdxWave].second);

			phaseData->SetPoint(point, mass, fitData->plottingPhases()[idxBin][idxMass][idxWave][jdxWave].first);
			phaseData->SetPointError(point, halfBin, fitData->plottingPhases()[idxBin][idxMass][idxWave][jdxWave].second);

			if(fitInformation->getBin(idxBin).sysFileNames().size() > 0) {
				const double minSR = fitData->sysPlottingSpinDensityMatrixElementsReal()[idxBin][idxMass][idxWave][jdxWave].first;
				const double maxSR = fitData->sysPlottingSpinDensityMatrixElementsReal()[idxBin][idxMass][idxWave][jdxWave].second;
				realSystematics->SetPoint(point, mass, (maxSR+minSR)/2.);
				realSystematics->SetPointError(point, halfBin, (maxSR-minSR)/2.);

				const double minSI = fitData->sysPlottingSpinDensityMatrixElementsImag()[idxBin][idxMass][idxWave][jdxWave].first;
				const double maxSI = fitData->sysPlottingSpinDensityMatrixElementsImag()[idxBin][idxMass][idxWave][jdxWave].second;
				imagSystematics->SetPoint(point, mass, (maxSI+minSI)/2.);
				imagSystematics->SetPointError(point, halfBin, (maxSI-minSI)/2.);

				const double minSP = fitData->sysPlottingPhases()[idxBin][idxMass][idxWave][jdxWave].first;
				const double maxSP = fitData->sysPlottingPhases()[idxBin][idxMass][idxWave][jdxWave].second;
				phaseSystematics->SetPoint(point, mass, (maxSP+minSP)/2.);
				phaseSystematics->SetPointError(point, halfBin, (maxSP-minSP)/2.);
			}
		}

		// plot fit, either over full or limited mass range
		const size_t firstPoint = rangePlotting ? (extraBinning*fitData->wavePairMassBinLimits()[idxBin][idxWave][jdxWave].first) : 0;
		const size_t lastPoint = rangePlotting ? (extraBinning*fitData->wavePairMassBinLimits()[idxBin][idxWave][jdxWave].second) : (extraBinning*(fitData->nrMassBins()[idxBin]-1));
		for(size_t point = firstPoint; point <= lastPoint; ++point) {
			const size_t idxMass = (point%extraBinning == 0) ? (point/extraBinning) : std::numeric_limits<size_t>::max();
			const double massStep = (fitData->massBinCenters()[idxBin][fitData->nrMassBins()[idxBin]-1] - fitData->massBinCenters()[idxBin][0]) / (fitData->nrMassBins()[idxBin] - 1) / extraBinning;
			const double mass = (idxMass != std::numeric_limits<size_t>::max()) ? fitData->massBinCenters()[idxBin][idxMass] : (fitData->massBinCenters()[idxBin][point/extraBinning] + (point%extraBinning) * massStep);

			const std::complex<double> element = fitModel->spinDensityMatrix(fitParameters, cache, idxWave, jdxWave, idxBin, mass, idxMass);
			realFit->SetPoint(point-firstPoint, mass, element.real());
			imagFit->SetPoint(point-firstPoint, mass, element.imag());
		}

		// keep track of phase over full mass range
		TGraph phaseFitAll;
		for(size_t point = 0; point <= (extraBinning*(fitData->nrMassBins()[idxBin]-1)); ++point) {
			const size_t idxMass = (point%extraBinning == 0) ? (point/extraBinning) : std::numeric_limits<size_t>::max();
			const double massStep = (fitData->massBinCenters()[idxBin][fitData->nrMassBins()[idxBin]-1] - fitData->massBinCenters()[idxBin][0]) / (fitData->nrMassBins()[idxBin] - 1) / extraBinning;
			const double mass = (idxMass != std::numeric_limits<size_t>::max()) ? fitData->massBinCenters()[idxBin][idxMass] : (fitData->massBinCenters()[idxBin][point/extraBinning] + (point%extraBinning) * massStep);

			const double phase = fitModel->phase(fitParameters, cache, idxWave, jdxWave, idxBin, mass, idxMass) * TMath::RadToDeg();

			if(point != 0) {
				int bestOffs = 0;
				double bestDiff = std::numeric_limits<double>::max();

				double x;
				double prev;
				phaseFitAll.GetPoint(point-1, x, prev);
				for(int offs = -5; offs < 6; ++offs) {
					if(std::abs(phase + offs*360. - prev) < bestDiff) {
						bestDiff = std::abs(phase + offs*360. - prev);
						bestOffs = offs;
					}
				}

				phaseFitAll.SetPoint(point, mass, phase + bestOffs*360.);
			} else {
				phaseFitAll.SetPoint(point, mass, phase);
			}
		}

		// rectify phase graphs
		for(size_t point = 0; point <= (extraBinning*(fitData->nrMassBins()[idxBin]-1)); ++point) {
			const size_t idxMass = (point%extraBinning == 0) ? (point/extraBinning) : std::numeric_limits<size_t>::max();
			const double massStep = (fitData->massBinCenters()[idxBin][fitData->nrMassBins()[idxBin]-1] - fitData->massBinCenters()[idxBin][0]) / (fitData->nrMassBins()[idxBin] - 1) / extraBinning;
			const double mass = (idxMass != std::numeric_limits<size_t>::max()) ? fitData->massBinCenters()[idxBin][idxMass] : (fitData->massBinCenters()[idxBin][point/extraBinning] + (point%extraBinning) * massStep);

			double x;
			double valueFit;
			phaseFitAll.GetPoint(point, x, valueFit);

			if(idxMass != std::numeric_limits<size_t>::max()) {
				int bestOffs = 0;
				double bestDiff = std::numeric_limits<double>::max();

				double data;
				phaseData->GetPoint(idxMass, x, data);
				for(int offs = -5; offs < 6; ++offs) {
					if(std::abs(data + offs*360. - valueFit) < bestDiff) {
						bestDiff = std::abs(data + offs*360. - valueFit);
						bestOffs = offs;
					}
				}

				phaseData->SetPoint(idxMass, x, data + bestOffs*360.);
				if(fitInformation->getBin(idxBin).sysFileNames().size() > 0) {
					phaseSystematics->GetPoint(idxMass, x, data);
					phaseSystematics->SetPoint(idxMass, x, data + bestOffs*360.);
				}
			}

			// check that this mass bin should be taken into account for this
			// combination of waves
			if(point < firstPoint or point > lastPoint) {
				continue;
			}

			phaseFit->SetPoint(point-firstPoint, mass, valueFit);
		}

		outDirectory->cd();
		real.Write();
		imag.Write();
		phase.Write();
	}


	void
	createPlotsFsmd(const rpwa::resonanceFit::dataConstPtr& fitData,
	                const rpwa::resonanceFit::modelConstPtr& fitModel,
	                const rpwa::resonanceFit::parameters& fitParameters,
	                rpwa::resonanceFit::cache& cache,
	                TDirectory* outDirectory,
	                const size_t extraBinning,
	                const size_t idxBin)
	{
		if(debug) {
			printDebug << "start creating plots for final-state mass-dependence." << std::endl;
		}

		TGraph graph;
		graph.SetName("finalStateMassDependence");
		graph.SetTitle("finalStateMassDependence");

		for(size_t point = 0; point <= (extraBinning*(fitData->nrMassBins()[idxBin]-1)); ++point) {
			const size_t idxMass = (point%extraBinning == 0) ? (point/extraBinning) : std::numeric_limits<size_t>::max();
			const double massStep = (fitData->massBinCenters()[idxBin][fitData->nrMassBins()[idxBin]-1] - fitData->massBinCenters()[idxBin][0]) / (fitData->nrMassBins()[idxBin] - 1) / extraBinning;
			const double mass = (idxMass != std::numeric_limits<size_t>::max()) ? fitData->massBinCenters()[idxBin][idxMass] : (fitData->massBinCenters()[idxBin][point/extraBinning] + (point%extraBinning) * massStep);

			graph.SetPoint(point, mass, std::norm(fitModel->getFsmd()->val(fitParameters, cache, idxBin, mass, idxMass)));
		}

		outDirectory->cd();
		graph.Write();
	}


}


void
rpwa::resonanceFit::createPlots(const rpwa::resonanceFit::informationConstPtr& fitInformation,
                                const rpwa::resonanceFit::dataConstPtr& fitData,
                                const rpwa::resonanceFit::modelConstPtr& fitModel,
                                const rpwa::resonanceFit::parameters& fitParameters,
                                rpwa::resonanceFit::cache& cache,
                                TFile* outFile,
                                const bool rangePlotting,
                                const size_t extraBinning)
{
	if(debug) {
		printDebug << "start creating plots." << std::endl;
	}

	const bool sameMassBinning = fitData->hasSameMassBinning();
	for(size_t idxBin = 0; idxBin < fitInformation->nrBins(); ++idxBin) {
		TDirectory* outDirectory = NULL;
		if(fitInformation->nrBins() == 1) {
			outDirectory = outFile;
		} else {
			std::ostringstream name;
			name << "bin" << idxBin;
			outDirectory = outFile->mkdir(name.str().c_str());
		}

		for(size_t idxWave = 0; idxWave < fitInformation->nrWaves(); ++idxWave) {
			createPlotsWave(fitInformation,
			                fitData,
			                fitModel,
			                fitParameters,
			                cache,
			                outDirectory,
			                rangePlotting,
			                extraBinning,
			                idxWave,
			                idxBin);
		}

		for(size_t idxWave = 0; idxWave < fitInformation->nrWaves(); ++idxWave) {
			for(size_t jdxWave = idxWave+1; jdxWave < fitInformation->nrWaves(); ++jdxWave) {
				createPlotsWavePair(fitInformation,
				                    fitData,
				                    fitModel,
				                    fitParameters,
				                    cache,
				                    outDirectory,
				                    rangePlotting,
				                    extraBinning,
				                    idxWave,
				                    jdxWave,
				                    idxBin);
			}
		}

		if(fitModel->getFsmd() and (fitModel->getFsmd()->getNrBins() != 1 or not sameMassBinning)) {
			createPlotsFsmd(fitData,
			                fitModel,
			                fitParameters,
			                cache,
			                outDirectory,
			                extraBinning,
			                idxBin);
		}
	}

	if(fitInformation->nrBins() != 1 and sameMassBinning and fitModel->isMappingEqualInAllBins()) {
		for(size_t idxWave = 0; idxWave < fitInformation->nrWaves(); ++idxWave) {
			createPlotsWaveSum(fitInformation,
			                   fitData,
			                   fitModel,
			                   fitParameters,
			                   cache,
			                   outFile,
			                   rangePlotting,
			                   extraBinning,
			                   idxWave);
		}
	}

	if(fitModel->getFsmd() and (fitModel->getFsmd()->getNrBins() == 1 and sameMassBinning)) {
		createPlotsFsmd(fitData,
		                fitModel,
		                fitParameters,
		                cache,
		                outFile,
		                extraBinning,
		                0);
	}

	if(debug) {
		printDebug << "finished creating plots." << std::endl;
	}
}


void
rpwa::resonanceFit::setDebug(const bool newDebug)
{
	debug = newDebug;
}


bool rpwa::resonanceFit::massDepFit::_debug = false;


rpwa::resonanceFit::massDepFit::massDepFit()
{
}


void
rpwa::resonanceFit::massDepFit::setDebug(bool debug)
{
	_debug = debug;
	rpwa::resonanceFit::setDebug(_debug);
}


bool
rpwa::resonanceFit::massDepFit::readConfig(const YAML::Node& configRoot,
                                           rpwa::resonanceFit::informationConstPtr& fitInformation,
                                           rpwa::resonanceFit::dataConstPtr& fitData,
                                           rpwa::resonanceFit::modelConstPtr& fitModel,
                                           rpwa::resonanceFit::parameters& fitParameters,
                                           rpwa::resonanceFit::parameters& fitParametersError,
                                           std::map<std::string, double>& fitQuality,
                                           std::vector<std::string>& freeParameters,
                                           const bool useBranchings,
                                           const rpwa::resonanceFit::function::useCovarianceMatrix useCovariance,
                                           const std::string& valTreeName,
                                           const std::string& valBranchName)
{
	// get information of fit quality if a previous fit was stored in the
	// configuration file
	fitQuality = readFitQuality(configRoot);

	// get information for which parameters to release in which order
	freeParameters = readFreeParameters(configRoot);

	// get fit results and waves to use in the resonance fit
	fitInformation = readInformation(configRoot);
	if(not fitInformation) {
		printErr << "error while reading 'input' in configuration file." << std::endl;
		throw;
	}

	// extract information from fit results
	boost::multi_array<std::string, 2> waveNames;
	std::vector<size_t> nrMassBins;
	boost::multi_array<double, 2> massBinCenters;
	boost::multi_array<double, 3> phaseSpaceIntegrals;
	boost::multi_array<std::complex<double>, 3> inProductionAmplitudes;
	boost::multi_array<TMatrixT<double>, 2> inProductionAmplitudesCovariance;
	boost::multi_array<std::complex<double>, 4> inSpinDensityMatrices;
	boost::multi_array<TMatrixT<double>, 2> inSpinDensityMatricesCovariance;
	boost::multi_array<std::pair<double, double>, 3> inPlottingIntensities;
	boost::multi_array<std::pair<double, double>, 4> inPlottingSpinDensityMatrixElementsReal;
	boost::multi_array<std::pair<double, double>, 4> inPlottingSpinDensityMatrixElementsImag;
	boost::multi_array<std::pair<double, double>, 4> inPlottingPhases;
	boost::multi_array<std::pair<double, double>, 3> inSysPlottingIntensities;
	boost::multi_array<std::pair<double, double>, 4> inSysPlottingSpinDensityMatrixElementsReal;
	boost::multi_array<std::pair<double, double>, 4> inSysPlottingSpinDensityMatrixElementsImag;
	boost::multi_array<std::pair<double, double>, 4> inSysPlottingPhases;
	readInFiles(fitInformation,
	            waveNames,
	            nrMassBins,
	            massBinCenters,
	            phaseSpaceIntegrals,
	            inProductionAmplitudes,
	            inProductionAmplitudesCovariance,
	            inSpinDensityMatrices,
	            inSpinDensityMatricesCovariance,
	            inPlottingIntensities,
	            inPlottingSpinDensityMatrixElementsReal,
	            inPlottingSpinDensityMatrixElementsImag,
	            inPlottingPhases,
	            inSysPlottingIntensities,
	            inSysPlottingSpinDensityMatrixElementsReal,
	            inSysPlottingSpinDensityMatrixElementsImag,
	            inSysPlottingPhases,
	            valTreeName,
	            valBranchName);

	// prepare mass limits
	boost::multi_array<std::pair<size_t, size_t>, 3> wavePairMassBinLimits;
	prepareMassLimits(fitInformation,
	                  nrMassBins,
	                  massBinCenters,
	                  wavePairMassBinLimits);

	// set-up fit model (resonances, background, final-state mass-dependence)
	fitModel = readModel(configRoot,
	                     fitInformation,
	                     fitParameters,
	                     fitParametersError,
	                     waveNames,
	                     nrMassBins,
	                     massBinCenters,
	                     phaseSpaceIntegrals,
	                     useBranchings);

	std::ostringstream output;
	for(size_t idxComponent = 0; idxComponent < fitModel->getNrComponents(); ++idxComponent) {
		output << "    " << fitModel->getComponent(idxComponent)->getName() << std::endl;
	}
	printInfo << "fitting " << fitModel->getNrComponents() << " components to the data:" << std::endl
	          << output.str();
	if(fitModel->getFsmd()) {
		printInfo << "using final-state mass-dependence." << std::endl;
	} else {
		printInfo << "not using final-state mass-dependence." << std::endl;
	}

	// prepare production amplitudes and corresponding covariance matrices
	// for the fit
	prepareProductionAmplitudes(useCovariance,
	                            fitModel->getAnchorWave(),
	                            wavePairMassBinLimits,
	                            inProductionAmplitudes,
	                            inProductionAmplitudesCovariance);
	prepareSpinDensityMatrices(useCovariance,
	                           wavePairMassBinLimits,
	                           inSpinDensityMatrices,
	                           inSpinDensityMatricesCovariance);

	// create data object
	fitData.reset(new rpwa::resonanceFit::data(nrMassBins,
	                                           massBinCenters,
	                                           wavePairMassBinLimits,
	                                           phaseSpaceIntegrals,
	                                           inProductionAmplitudes,
	                                           inProductionAmplitudesCovariance,
	                                           inSpinDensityMatrices,
	                                           inSpinDensityMatricesCovariance,
	                                           useCovariance,
	                                           inPlottingIntensities,
	                                           inPlottingSpinDensityMatrixElementsReal,
	                                           inPlottingSpinDensityMatrixElementsImag,
	                                           inPlottingPhases,
	                                           inSysPlottingIntensities,
	                                           inSysPlottingSpinDensityMatrixElementsReal,
	                                           inSysPlottingSpinDensityMatrixElementsImag,
	                                           inSysPlottingPhases));

	return true;
}


bool
rpwa::resonanceFit::massDepFit::init(const rpwa::resonanceFit::dataConstPtr& fitData,
                                     const rpwa::resonanceFit::modelConstPtr& fitModel,
                                     const rpwa::resonanceFit::functionPtr& fitFunction)
{
	if(not fitFunction->init(fitData,
	                         fitModel)) {
		printErr << "error while initializing the function to minimize." << std::endl;
		return false;
	}

	return true;
}


bool
rpwa::resonanceFit::massDepFit::writeConfig(std::ostream& output,
                                            const rpwa::resonanceFit::informationConstPtr& fitInformation,
                                            const rpwa::resonanceFit::modelConstPtr& fitModel,
                                            const rpwa::resonanceFit::parameters& fitParameters,
                                            const rpwa::resonanceFit::parameters& fitParametersError,
                                            const std::map<std::string, double>& fitQuality,
                                            const std::vector<std::string>& freeParameters) const
{
	if(_debug) {
		printDebug << "writing configuration file." << std::endl;
	}

	YAML::Emitter yamlOutput(output);
	yamlOutput << YAML::BeginMap;

	writeFitQuality(yamlOutput, fitQuality);
	writeFreeParameters(yamlOutput, freeParameters);
	writeInformation(yamlOutput, fitInformation);
	writeModel(yamlOutput, fitModel, fitParameters, fitParametersError);

	yamlOutput << YAML::EndMap;

	// newline at end-of-file
	output << std::endl;

	return true;
}
