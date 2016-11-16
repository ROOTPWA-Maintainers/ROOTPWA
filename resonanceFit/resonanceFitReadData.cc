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


#include "resonanceFit.h"
#include "resonanceFitInternal.h"

#include <boost/multi_array.hpp>

#include <TFile.h>
#include <TMatrixT.h>
#include <TTree.h>

#include <fitResult.h>
#include <reportingUtils.hpp>

#include "data.h"
#include "information.h"
#include "resonanceFitHelper.h"


namespace {


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

		for(size_t idxBin = 0; idxBin < nrBins; ++idxBin) {
			// get range of used mass bins (in any wave) for each bin
			size_t idxMassMin = maxMassBins;
			size_t idxMassMax = 0;
			for(size_t idxWave = 0; idxWave < nrWaves; ++idxWave) {
				idxMassMin = std::min(idxMassMin, wavePairMassBinLimits[idxBin][idxWave][idxWave].first);
				idxMassMax = std::max(idxMassMax, wavePairMassBinLimits[idxBin][idxWave][idxWave].second);
			}

			for(size_t idxMass = idxMassMin; idxMass <= idxMassMax; ++idxMass) {
				// get a list of waves that are zero (those have to be excluded
				// from the inversion of the covariance matrix below) and test
				// that the anchor wave is non-zero over the complete fit range
				bool zeroAnchorWave = false;
				std::vector<size_t> zeroWaves;
				// test if the anchor wave is real valued
				bool realAnchorWave = true;
				// test that any non-anchor wave is not real valued
				bool realOtherWaves = false;
				for(size_t idxWave = 0; idxWave < nrWaves; ++idxWave) {
					bool zeroThisWave = true;
					zeroThisWave &= (productionAmplitudes[idxBin][idxMass][idxWave].real() == 0.);
					zeroThisWave &= (productionAmplitudes[idxBin][idxMass][idxWave].imag() == 0.);
					zeroThisWave &= (productionAmplitudesCovariance[idxBin][idxMass](2*idxWave,   2*idxWave  ) == 0.);
					zeroThisWave &= (productionAmplitudesCovariance[idxBin][idxMass](2*idxWave,   2*idxWave+1) == 0.);
					zeroThisWave &= (productionAmplitudesCovariance[idxBin][idxMass](2*idxWave+1, 2*idxWave  ) == 0.);
					zeroThisWave &= (productionAmplitudesCovariance[idxBin][idxMass](2*idxWave+1, 2*idxWave+1) == 0.);

					bool realThisWave = true;
					realThisWave &= (productionAmplitudes[idxBin][idxMass][idxWave].imag() == 0.);
					realThisWave &= (productionAmplitudesCovariance[idxBin][idxMass](2*idxWave,   2*idxWave+1) == 0.);
					realThisWave &= (productionAmplitudesCovariance[idxBin][idxMass](2*idxWave+1, 2*idxWave  ) == 0.);
					realThisWave &= (productionAmplitudesCovariance[idxBin][idxMass](2*idxWave+1, 2*idxWave+1) == 0.);

					if(zeroThisWave or idxMass < wavePairMassBinLimits[idxBin][idxWave][idxWave].first or idxMass > wavePairMassBinLimits[idxBin][idxWave][idxWave].second) {
						zeroWaves.push_back(idxWave);
					}

					if(idxWave == idxAnchorWave) {
						zeroAnchorWave |= zeroThisWave;
						realAnchorWave &= realThisWave;
					} else if(not zeroThisWave) {
						realOtherWaves |= realThisWave;
					}

					// check that a wave is not zero in its fit range
					if(zeroThisWave and idxMass >= wavePairMassBinLimits[idxBin][idxWave][idxWave].first and idxMass <= wavePairMassBinLimits[idxBin][idxWave][idxWave].second) {
						printErr << "production amplitudes of wave " << idxWave << " zero in its fit range (e.g. mass limit in mass-independent fit)." << std::endl;
						return false;
					}
				}

				// error if anchor wave is zero in one mass bin
				if(zeroAnchorWave) {
					printErr << "production amplitudes of anchor wave zero in some mass bins (mass limit in mass-independent fit)." << std::endl;
					return false;
				}

				// error if any non-anchor wave is real
				if(realOtherWaves) {
					printErr << "production amplitudes cannot be fitted if a non-anchor wave is real valued." << std::endl;
					return false;
				}

				// determine whether the anchor wave should be used in the current bin
				bool skipAnchor = false;
				for(std::vector<size_t>::const_iterator it = zeroWaves.begin(); it != zeroWaves.end(); ++it) {
					if(*it == idxAnchorWave) {
						skipAnchor = true;
					}
				}

				// import covariance matrix of production amplitudes
				const size_t matrixSize = 2 * (nrWaves - zeroWaves.size()) - (skipAnchor ? 0 : 1);
				TMatrixT<double> reducedCovMat(matrixSize, matrixSize);

				if(realAnchorWave) {
					for(size_t idxWave = 0, idxSkip = 0; idxWave < nrWaves; ++idxWave) {
						if(idxSkip < zeroWaves.size() and zeroWaves[idxSkip] == idxWave) {
							++idxSkip;
							continue;
						}

						for(size_t jdxWave = 0, jdxSkip = 0; jdxWave < nrWaves; ++jdxWave) {
							if(jdxSkip < zeroWaves.size() and zeroWaves[jdxSkip] == jdxWave) {
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
					// rotate production amplitudes and
					// covariance matrices such that the
					// anchor wave is real
					TMatrixT<double> covariance(matrixSize + (skipAnchor ? 0 : 1), matrixSize + (skipAnchor ? 0 : 1));
					TMatrixT<double> jacobian(matrixSize, matrixSize + (skipAnchor ? 0 : 1));

					for(size_t idxWave = 0, idxSkip = 0; idxWave < nrWaves; ++idxWave) {
						if(idxSkip < zeroWaves.size() and zeroWaves[idxSkip] == idxWave) {
							++idxSkip;
							continue;
						}

						for(size_t jdxWave = 0, jdxSkip = 0; jdxWave < nrWaves; ++jdxWave) {
							if(jdxSkip < zeroWaves.size() and zeroWaves[jdxSkip] == jdxWave) {
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

					// modify measured production amplitude such that the anchor wave is always real and positive
					const std::complex<double> anchorPhase = productionAmplitudes[idxBin][idxMass][idxAnchorWave] / abs(productionAmplitudes[idxBin][idxMass][idxAnchorWave]);
					for(size_t idxWave = 0; idxWave < nrWaves; ++idxWave) {
						productionAmplitudes[idxBin][idxMass][idxWave] /= anchorPhase;
					}
				}

				// set entries in covariance matrix to zero according to which parts are to be used
				if(useCovariance != rpwa::resonanceFit::function::useFullCovarianceMatrix) {
					for(size_t idxWave = 0, idxSkip = 0; idxWave < nrWaves; ++idxWave) {
						if(idxSkip < zeroWaves.size() and zeroWaves[idxSkip] == idxWave) {
							++idxSkip;
							continue;
						}

						for(size_t jdxWave = 0, jdxSkip = 0; jdxWave < nrWaves; ++jdxWave) {
							if(jdxSkip < zeroWaves.size() and zeroWaves[jdxSkip] == jdxWave) {
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
					if(idxSkip < zeroWaves.size() and zeroWaves[idxSkip] == idxWave) {
						++idxSkip;
						continue;
					}

					for(size_t jdxWave = 0, jdxSkip = 0; jdxWave < nrWaves; ++jdxWave) {
						if(jdxSkip < zeroWaves.size() and zeroWaves[jdxSkip] == jdxWave) {
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

			// remove (set to zero size) the unused covariance matrices
			for(size_t idxMass = 0; idxMass < idxMassMin; ++idxMass) {
				productionAmplitudesCovariance[idxBin][idxMass].ResizeTo(0, 0);
			}
			for(size_t idxMass = idxMassMax + 1; idxMass < maxMassBins; ++idxMass) {
				productionAmplitudesCovariance[idxBin][idxMass].ResizeTo(0, 0);
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

		for(size_t idxBin = 0; idxBin < nrBins; ++idxBin) {
			// get range of used mass bins (in any wave) for each bin
			size_t idxMassMin = maxMassBins;
			size_t idxMassMax = 0;
			for(size_t idxWave = 0; idxWave < nrWaves; ++idxWave) {
				idxMassMin = std::min(idxMassMin, wavePairMassBinLimits[idxBin][idxWave][idxWave].first);
				idxMassMax = std::max(idxMassMax, wavePairMassBinLimits[idxBin][idxWave][idxWave].second);
			}

			for(size_t idxMass = idxMassMin; idxMass <= idxMassMax; ++idxMass) {
				// get a list of waves that are zero (those have to be excluded
				// from the inversion of the covariance matrix below)
				std::vector<size_t> zeroWaves;
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
						zeroWaves.push_back(idxWave);
					}

					// check that a wave is not zero in its fit range
					if(zeroThisWave and idxMass >= wavePairMassBinLimits[idxBin][idxWave][idxWave].first and idxMass <= wavePairMassBinLimits[idxBin][idxWave][idxWave].second) {
						printErr << "spin-density matrix element of wave " << idxWave << " zero in its fit range (e.g. mass limit in mass-independent fit)." << std::endl;
						return false;
					}
				}

				// import covariance matrix of spin-density matrix elements
				const size_t reducedMatrixSize((nrWaves - zeroWaves.size())*(nrWaves - zeroWaves.size()));
				TMatrixT<double> reducedCovMat(reducedMatrixSize, reducedMatrixSize);

				{
					// i is for loop over rows
					size_t redIdx = 0;
					for(size_t iWave1 = 0, iSkip1 = 0; iWave1 < nrWaves; ++iWave1) {
						if(iSkip1 < zeroWaves.size() and zeroWaves[iSkip1] == iWave1) {
							++iSkip1;
							continue;
						}
						for(size_t iWave2 = iWave1, iSkip2 = iSkip1; iWave2 < nrWaves; ++iWave2) {
							if(iSkip2 < zeroWaves.size() and zeroWaves[iSkip2] == iWave2) {
								++iSkip2;
								continue;
							}
							const size_t idx = nrWaves*(nrWaves+1) - (nrWaves-iWave1)*(nrWaves-iWave1+1) + 2*(iWave2-iWave1);

							// j is for loop over columns
							size_t redJdx = 0;
							for(size_t jWave1 = 0, jSkip1 = 0; jWave1 < nrWaves; ++jWave1) {
								if(jSkip1 < zeroWaves.size() and zeroWaves[jSkip1] == jWave1) {
									++jSkip1;
									continue;
								}
								for(size_t jWave2 = jWave1, jSkip2 = jSkip1; jWave2 < nrWaves; ++jWave2) {
									if(jSkip2 < zeroWaves.size() and zeroWaves[jSkip2] == jWave2) {
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
					if(iSkip1 < zeroWaves.size() and zeroWaves[iSkip1] == iWave1) {
						++iSkip1;
						continue;
					}
					for(size_t iWave2 = iWave1, iSkip2 = iSkip1; iWave2 < nrWaves; ++iWave2) {
						if(iSkip2 < zeroWaves.size() and zeroWaves[iSkip2] == iWave2) {
							++iSkip2;
							continue;
						}
						const size_t idx = nrWaves*(nrWaves+1) - (nrWaves-iWave1)*(nrWaves-iWave1+1) + 2*(iWave2-iWave1);

						// j is for loop over columns
						size_t redJdx = 0;
						for(size_t jWave1 = 0, jSkip1 = 0; jWave1 < nrWaves; ++jWave1) {
							if(jSkip1 < zeroWaves.size() and zeroWaves[jSkip1] == jWave1) {
								++jSkip1;
								continue;
							}
							for(size_t jWave2 = jWave1, jSkip2 = jSkip1; jWave2 < nrWaves; ++jWave2) {
								if(jSkip2 < zeroWaves.size() and zeroWaves[jSkip2] == jWave2) {
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

			// remove (set to zero size) the unused covariance matrices
			for(size_t idxMass = 0; idxMass < idxMassMin; ++idxMass) {
				spinDensityMatricesCovariance[idxBin][idxMass].ResizeTo(0, 0);
			}
			for(size_t idxMass = idxMassMax + 1; idxMass < maxMassBins; ++idxMass) {
				spinDensityMatricesCovariance[idxBin][idxMass].ResizeTo(0, 0);
			}
		}

		return true;
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

		if(rpwa::resonanceFit::debug()) {
			printDebug << "getting center of mass bins from " << nrEntries << " entries in tree." << std::endl;
		}

		nrMassBins = 0;
		for(Long64_t idx = 0; idx < nrEntries; ++idx) {
			if(tree->GetEntry(idx) == 0) {
				printErr << "error while reading entry " << idx << " from tree." << std::endl;
				throw;
			}
			const double newMass = fit->massBinCenter();

			if(rpwa::resonanceFit::debug()) {
				printDebug << "entry " << idx << ": center of mass bin at " << newMass << " GeV/c^2" << std::endl;
			}

			bool found = false;
			for(size_t idxMass = 0; idxMass < nrMassBins; ++idxMass) {
				if(std::abs(massBinCenters[idxMass]-newMass) < 1000.*std::numeric_limits<double>::epsilon()) {
					found = true;
					if(rpwa::resonanceFit::debug()) {
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

		if(rpwa::resonanceFit::debug()) {
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

			if(rpwa::resonanceFit::debug()) {
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

			if(rpwa::resonanceFit::debug()) {
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

		if(rpwa::resonanceFit::debug()) {
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

		if(rpwa::resonanceFit::debug()) {
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

		if(rpwa::resonanceFit::debug()) {
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
			if(rpwa::resonanceFit::debug()) {
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

			if(rpwa::resonanceFit::debug()) {
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

		if(rpwa::resonanceFit::debug()) {
			printDebug << "reading phase-space integrals for " << waveNames.size() << " waves from fit result." << std::endl;
		}

		for(size_t idxMass = 0; idxMass < mapping.size(); ++idxMass) {
			if(rpwa::resonanceFit::debug()) {
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

		if(rpwa::resonanceFit::debug()) {
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
		if(rpwa::resonanceFit::debug()) {
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

		if(rpwa::resonanceFit::debug()) {
			printDebug << "searching for tree '" << valTreeName << "' in file '" << bin.fileName() << "'." << std::endl;
		}

		TTree* inTree;
		inFile->GetObject(valTreeName.c_str(), inTree);
		if(not inTree) {
			printErr << "input tree '" << valTreeName << "' not found in input file '" << bin.fileName() << "'."<< std::endl;
			throw;
		}

		if(rpwa::resonanceFit::debug()) {
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
		if(rpwa::resonanceFit::debug()) {
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

		if(rpwa::resonanceFit::debug()) {
			printDebug << "searching for tree '" << valTreeName << "' in file '" << bin.sysFileNames()[idxSystematics] << "'." << std::endl;
		}

		TTree* sysTree;
		sysFile->GetObject(valTreeName.c_str(), sysTree);
		if(not sysTree) {
			printErr << "input tree '" << valTreeName << "' not found in input file '" << bin.sysFileNames()[idxSystematics] << "'."<< std::endl;
			throw;
		}

		if(rpwa::resonanceFit::debug()) {
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
		if(rpwa::resonanceFit::debug()) {
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
	prepareMassLimit(const rpwa::resonanceFit::informationConstPtr& fitInformation,
	                 const size_t nrMassBins,
	                 const boost::multi_array<double, 1>& massBinCenters,
	                 boost::multi_array<std::pair<size_t, size_t>, 2>& wavePairMassBinLimits)
	{
		if(rpwa::resonanceFit::debug()) {
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
			if(rpwa::resonanceFit::debug()) {
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

		if(rpwa::resonanceFit::debug()) {
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
		if(rpwa::resonanceFit::debug()) {
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


}


rpwa::resonanceFit::dataConstPtr
rpwa::resonanceFit::readData(const rpwa::resonanceFit::informationConstPtr& fitInformation,
                             const std::string& anchorWaveName,
                             const rpwa::resonanceFit::function::useCovarianceMatrix useCovariance,
                             const std::string& valTreeName,
                             const std::string& valBranchName)
{
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

	// prepare production amplitudes and corresponding covariance matrices
	// for the fit
	std::vector<std::string> mainWaveNames(fitInformation->nrWaves());
	std::transform(fitInformation->waves().begin(), fitInformation->waves().end(), mainWaveNames.begin(), [](const rpwa::resonanceFit::information::wave& wave){ return wave.waveName(); });
	const size_t idxAnchorWave = std::find(mainWaveNames.begin(), mainWaveNames.end(), anchorWaveName) - mainWaveNames.begin();
	if(idxAnchorWave >= mainWaveNames.size()) {
		printErr << "anchor wave '" << anchorWaveName << "' not found in fit." << std::endl;
		throw;
	}
	prepareProductionAmplitudes(useCovariance,
	                            idxAnchorWave,
	                            wavePairMassBinLimits,
	                            inProductionAmplitudes,
	                            inProductionAmplitudesCovariance);
	prepareSpinDensityMatrices(useCovariance,
	                           wavePairMassBinLimits,
	                           inSpinDensityMatrices,
	                           inSpinDensityMatricesCovariance);

	// create data object
	return std::make_shared<rpwa::resonanceFit::data>(nrMassBins,
	                                                  massBinCenters,
	                                                  waveNames,
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
	                                                  inSysPlottingPhases);
}
