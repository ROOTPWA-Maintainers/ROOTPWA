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
#include "model.h"


namespace {


	template<typename T>
	void
	copyMultiArrayElement(const T& source,
	                      T& target)
	{
		target = source;
	}

	void
	copyMultiArrayElement(const TMatrixT<double>& source,
	                      TMatrixT<double>& target)
	{
		target.ResizeTo(source);
		target = source;
	}


	template<typename T, size_t dim1, size_t dim2>
	void
	copyMultiArray(const boost::multi_array<T, dim1>& source,
	               std::vector<size_t> sourceIndices,
	               boost::multi_array<T, dim2>& target,
	               std::vector<size_t> targetIndices)
	{
		assert(sourceIndices.size() <= dim1);
		assert(targetIndices.size() <= dim2);
		assert(dim1-sourceIndices.size() == dim2-targetIndices.size());
		assert(((sourceIndices.size() == dim1) and (targetIndices.size() == dim2))
		       or ((sourceIndices.size() != dim1) and (targetIndices.size() != dim2)));

		if(sourceIndices.size() == dim1) {
			copyMultiArrayElement(source(sourceIndices), target(targetIndices));
		} else {
			const size_t maxIdx = source.shape()[sourceIndices.size()];

			sourceIndices.push_back(0);
			targetIndices.push_back(0);
			for(size_t idx = 0; idx < maxIdx; ++idx) {
				sourceIndices.back() = idx;
				targetIndices.back() = idx;
				copyMultiArray(source, sourceIndices, target, targetIndices);
			}
		}
	}


	// increase the extensions of a boost::multi_array such that one part
	// of it can simply be reset
	template<typename T>
	void
	adjustSize(boost::multi_array<T, 1>& master,
	           const size_t minSizeFirst = 1)
	{
		// initialize the next size with the current size
		std::vector<size_t> newSize(master.shape(), master.shape()+master.num_dimensions());

		// resize if the minimal size required for the first dimension
		// is larger than its current size
		if(newSize[0] < minSizeFirst) {
			newSize[0] = minSizeFirst;
			master.resize(newSize);
		}
	}


	// increase the extensions of a boost::multi_array such that one part
	// of it can simply be reset
	template<typename T, size_t dim>
	void
	adjustSize(boost::multi_array<T, dim>& master,
	           const boost::multi_array<T, dim-1>& part,
	           const size_t minSizeFirst = 1)
	{
		// initialize the next size with the current size
		std::vector<size_t> newSize(master.shape(), master.shape()+master.num_dimensions());

		// resize if the minimal size required for the first dimension
		// is larger than its current size
		bool resize = newSize[0] < minSizeFirst;
		newSize[0] = std::max(newSize[0], minSizeFirst);

		// compare the other dimensions with the dimenstions of the
		// part to be set
		for(size_t i = 1; i < dim; ++i) {
			if(newSize[i] < part.shape()[i-1]) {
				resize = true;
				newSize[i] = part.shape()[i-1];
			}
		}

		if(resize) {
			master.resize(newSize);
		}
	}


	// increase the extensions of a boost::multi_array such that one part
	// of it can simply be reset
	// special case for TMatrixT because TMatrixT::operator= does not
	// adjust the size of the matrices
	template<size_t dim>
	void
	adjustSize(boost::multi_array<TMatrixT<double>, dim>& master,
	           const boost::multi_array<TMatrixT<double>, dim-1>& part,
	           const size_t minSizeFirst = 1)
	{
		// initialize the next size with the current size
		std::vector<size_t> newSize(master.shape(), master.shape()+master.num_dimensions());

		// resize if the minimal size required for the first dimension
		// is larger than its current size
		bool resize = newSize[0] < minSizeFirst;
		newSize[0] = std::max(newSize[0], minSizeFirst);

		// compare the other dimensions with the dimenstions of the
		// part to be set
		for(size_t i = 1; i < dim; ++i) {
			if(newSize[i] < part.shape()[i-1]) {
				resize = true;
				newSize[i] = part.shape()[i-1];
			}
		}

		if(resize) {
			boost::multi_array<TMatrixT<double>, dim> temp(std::vector<size_t>(master.shape(), master.shape()+master.num_dimensions()));

			copyMultiArray(master, std::vector<size_t>(), temp, std::vector<size_t>());

			// clear the current content, in particular make sure
			// that after the next 'resize' all TMatrixT have
			// dimension 0x0
			master.resize(std::vector<size_t>(master.num_dimensions(), 0));

			// resize to new size
			master.resize(newSize);

			copyMultiArray(temp, std::vector<size_t>(), master, std::vector<size_t>());
		}
	}


	// increase the extensions of a std::vector such that a new element can
	// simply be set
	template<typename T>
	void
	adjustSize(std::vector<T>& master,
	           const size_t minSize)
	{
		if(master.size() < minSize) {
			master.resize(minSize);
		}
	}


	// increase the extensions of a boost::multi_array such that one part
	// of it can simply be reset, and also set this part
	template<typename T>
	void
	adjustSizeAndSet(boost::multi_array<T, 1>& master,
	                 const size_t idx,
	                 const T& part)
	{
		adjustSize(master, idx+1);
		master[idx] = part;
	}


	// increase the extensions of a boost::multi_array such that one part
	// of it can simply be reset, and also set this part
	template<typename T, size_t dim>
	void
	adjustSizeAndSet(boost::multi_array<T, dim>& master,
	                 const size_t idx,
	                 const boost::multi_array<T, dim-1>& part)
	{
		adjustSize(master, part, idx+1);
		copyMultiArray(part, std::vector<size_t>(), master, std::vector<size_t>(1, idx));
	}


	// increase the extensions of a std::vector such that a new element can
	// simply be set, and also set this part
	template<typename T>
	void
	adjustSizeAndSet(std::vector<T>& master,
	                 const size_t idx,
	                 const T& part)
	{
		adjustSize(master, idx+1);
		master[idx] = part;
	}


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


}


bool rpwa::resonanceFit::massDepFit::_debug = false;


rpwa::resonanceFit::massDepFit::massDepFit()
	: _sameMassBinning(true),
	  _nrBins(0),
	  _maxMassBins(0),
	  _nrWaves(0)
{
}


bool
rpwa::resonanceFit::massDepFit::readConfig(const YAML::Node& configRoot,
                                           rpwa::resonanceFit::dataConstPtr& fitData,
                                           const rpwa::resonanceFit::modelPtr& fitModel,
                                           rpwa::resonanceFit::parameters& fitParameters,
                                           rpwa::resonanceFit::parameters& fitParametersError,
                                           std::map<std::string, double>& fitQuality,
                                           const bool useBranchings,
                                           const rpwa::resonanceFit::function::useCovarianceMatrix useCovariance,
                                           const std::string& valTreeName,
                                           const std::string& valBranchName)
{
	// fit result information
	const YAML::Node& configFitquality = configRoot["fitquality"];
	if(configFitquality) {
		if(not readConfigFitquality(configFitquality, fitQuality)) {
			printErr << "error while reading 'fitquality' in configuration file." << std::endl;
			return false;
		}
	} else {
		fitQuality.clear();
	}

	// input section
	const YAML::Node& configInput = configRoot["input"];
	if(not configInput) {
		printErr << "'input' does not exist in configuration file." << std::endl;
		return false;
	}
	if(not readConfigInput(configInput)) {
		printErr << "error while reading 'input' in configuration file." << std::endl;
		return false;
	}

	// extract information from fit results
	std::vector<size_t> nrMassBins;
	boost::multi_array<double, 2> massBinCenters;
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
	if(not readInFiles(nrMassBins,
	                   massBinCenters,
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
	                   valBranchName)) {
		printErr << "error while reading fit result." << std::endl;
		return false;
	}

	// prepare mass limits
	if(not prepareMassLimits(nrMassBins, massBinCenters)) {
		printErr << "error while determining which bins to use in the fit." << std::endl;
		return false;
	}

	// set-up fit model (resonances, background, final-state mass-dependence)
	const YAML::Node& configModel = configRoot["model"];
	if(not configModel) {
		printErr << "'model' does not exist in configuration file." << std::endl;
		return false;
	}
	if(not readConfigModel(configModel,
	                       fitModel,
	                       fitParameters,
	                       fitParametersError,
	                       nrMassBins,
	                       massBinCenters,
	                       useBranchings)) {
		printErr << "error while reading 'model' in configuration file." << std::endl;
		return false;
	}

	// prepare production amplitudes and corresponding covariance matrices
	// for the fit
	const size_t idxAnchorWave = std::find(_waveNames.begin(), _waveNames.end(), _anchorWaveName) - _waveNames.begin();
	if(idxAnchorWave >= _waveNames.size()) {
		printErr << "anchor wave '" << _anchorWaveName << "' not found in fit." << std::endl;
		return false;
	}
	prepareProductionAmplitudes(useCovariance,
	                            idxAnchorWave,
	                            _wavePairMassBinLimits,
	                            inProductionAmplitudes,
	                            inProductionAmplitudesCovariance);
	prepareSpinDensityMatrices(useCovariance,
	                           _wavePairMassBinLimits,
	                           inSpinDensityMatrices,
	                           inSpinDensityMatricesCovariance);

	// create data object
	fitData.reset(new rpwa::resonanceFit::data(nrMassBins,
	                                           massBinCenters,
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
rpwa::resonanceFit::massDepFit::readConfigFitquality(const YAML::Node& configFitquality,
                                                     std::map<std::string, double>& fitQuality) const
{
	// clear all information previously stored
	fitQuality.clear();

	if(not configFitquality) {
		printErr << "'configFitquality' is not a valid YAML node." << std::endl;
		return false;
	}
	if(not configFitquality.IsMap()) {
		printErr << "'fitquality' is not a YAML map." << std::endl;
		return false;
	}

	for(YAML::const_iterator it = configFitquality.begin(); it != configFitquality.end(); ++it) {
		if(not checkVariableType(it->first, YamlCppUtils::TypeString) or not checkVariableType(it->second, YamlCppUtils::TypeFloat)) {
			printErr << "entries in 'fitquality' must be pairs of 'string' and 'double'." << std::endl;
			return false;
		}

		const std::string key = it->first.as<std::string>();
		const double value = it->second.as<double>();

		if(fitQuality.count(key) != 0) {
			printErr << "variable '" << key << "' of 'fitquality' given multiple times." << std::endl;
			return false;
		}

		fitQuality[key] = value;
	}

	return true;
}


bool
rpwa::resonanceFit::massDepFit::readConfigInput(const YAML::Node& configInput)
{
	if(not configInput) {
		printErr << "'configInput' is not a valid YAML node." << std::endl;
		return false;
	}
	if(not configInput.IsMap()) {
		printErr << "'input' is not a YAML map." << std::endl;
		return false;
	}

	// get information about fit results from mass-independent
	const YAML::Node& configInputFitResults = configInput["fitresults"];
	if(not configInputFitResults) {
		printErr << "'fitresults' does not exist in 'input'." << std::endl;
		return false;
	}
	if(not readConfigInputFitResults(configInputFitResults)) {
		printErr << "error while reading 'fitresults' in 'input'." << std::endl;
		return false;
	}

	// get information about waves to be used in the fit
	const YAML::Node& configInputWaves = configInput["waves"];
	if(not configInputWaves) {
		printErr << "'waves' does not exist in 'input'." << std::endl;
		return false;
	}
	if(not readConfigInputWaves(configInputWaves)) {
		printErr << "error while reading 'waves' in 'input'." << std::endl;
		return false;
	}

	// get information for which parameters to release in which order
	const YAML::Node& configInputFreeParameters = configInput["freeparameters"];
	if(configInputFreeParameters) {
		if(not readConfigInputFreeParameters(configInputFreeParameters)) {
			printErr << "error while reading 'freeparameters' in 'input'." << std::endl;
			return false;
		}
	} else {
		_freeParameters.clear();
		_freeParameters.push_back("coupling branching");
		_freeParameters.push_back("coupling branching mass m0");
		_freeParameters.push_back("*");
	}

	return true;
}


bool
rpwa::resonanceFit::massDepFit::readConfigInputFitResults(const YAML::Node& configInputFitResults)
{
	if(not configInputFitResults) {
		printErr << "'configInputFitResults' is not a valid YAML node." << std::endl;
		return false;
	}
	if(not configInputFitResults.IsSequence()) {
		printErr << "'fitresults' is not a YAML sequence." << std::endl;
		return false;
	}

	_inFileName.clear();

	const size_t nrFitResults = configInputFitResults.size();
	for(size_t idxFitResult=0; idxFitResult<nrFitResults; ++idxFitResult) {
		if(_debug) {
			printDebug << "reading of entry " << idxFitResult << " in 'fitresults'." << std::endl;
		}

		const YAML::Node& configInputFitResult = configInputFitResults[idxFitResult];

		std::map<std::string, YamlCppUtils::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("name", YamlCppUtils::TypeString)
		                     ("tPrimeMean", YamlCppUtils::TypeFloat);
		if(not checkIfAllVariablesAreThere(configInputFitResult, mandatoryArguments)) {
			printErr << "'fitresults' entry at index " << idxFitResult << " does not contain all required variables." << std::endl;
			return false;
		}

		const std::string fileName = configInputFitResult["name"].as<std::string>();
		_inFileName.push_back(fileName);

		if(_debug) {
			printDebug << "read file name of fit results of mass-independent fit: '" << fileName << "'." << std::endl;
		}

		const double tPrimeMean = configInputFitResult["tPrimeMean"].as<double>();
		_tPrimeMeans.push_back(tPrimeMean);

		if(_debug) {
			printDebug << "read mean t' value: '" << tPrimeMean << "'." << std::endl;
		}

		if(configInputFitResult["rescaleErrors"]) {
			if(checkVariableType(configInputFitResult["rescaleErrors"], YamlCppUtils::TypeFloat)) {
				_rescaleErrors.push_back(configInputFitResult["rescaleErrors"].as<double>());
			} else {
				printErr << "variable 'rescaleErrors' of 'fitresults' entry at index " << idxFitResult << " is not a floating point number." << std::endl;
				return false;
			}
		} else {
			_rescaleErrors.push_back(1.);
		}

		if(_debug) {
			printDebug << "rescale errors by factor: '" << _rescaleErrors.back() << "'." << std::endl;
		}

		// get information for plotting of systematic error
		const YAML::Node& configInputFitResultSystematics = configInputFitResult["systematics"];
		if(configInputFitResultSystematics) {
			if(not readConfigInputFitResultSystematics(configInputFitResultSystematics)) {
				printErr << "error while reading 'systematics' in 'input'." << std::endl;
				return false;
			}
		} else {
			_nrSystematics.push_back(0);
			_sysFileNames.push_back(std::vector<std::string>());
		}
	}

	_nrBins = _inFileName.size();

	return true;
}


bool
rpwa::resonanceFit::massDepFit::readConfigInputFitResultSystematics(const YAML::Node& configInputFitResultSystematics)
{
	if(not configInputFitResultSystematics) {
		printErr << "'configInputSystematics' is not a valid YAML node." << std::endl;
		return false;
	}
	if(not configInputFitResultSystematics.IsSequence()) {
		printErr << "'systematics' is not a YAML sequence." << std::endl;
		return false;
	}

	_nrSystematics.push_back(configInputFitResultSystematics.size());
	if(_debug) {
		printDebug << "going to read information for " << _nrSystematics.back() << " files containing information for systematic errors." << std::endl;
	}

	std::vector<std::string> sysFileNames;
	for(size_t idxSystematics = 0; idxSystematics < _nrSystematics.back(); ++idxSystematics) {
		if(not checkVariableType(configInputFitResultSystematics[idxSystematics], YamlCppUtils::TypeString)) {
			printErr << "'systematics' entry at index " << idxSystematics << " is not a string." << std::endl;
			return false;
		}

		const std::string fileName = configInputFitResultSystematics[idxSystematics].as<std::string>();
		if(_debug) {
			printDebug << "'" << fileName << "' will be read to get information for systematic errors." << std::endl;
		}
		sysFileNames.push_back(fileName);
	}
	_sysFileNames.push_back(sysFileNames);

	if(_debug) {
		printDebug << "in total " << _nrSystematics.back() << " files to be read to get information for systematic errors." << std::endl;
	}

	return true;
}


bool
rpwa::resonanceFit::massDepFit::readConfigInputWaves(const YAML::Node& configInputWaves)
{
	if(not configInputWaves) {
		printErr << "'configInputWaves' is not a valid YAML node." << std::endl;
		return false;
	}
	if(not configInputWaves.IsSequence()) {
		printErr << "'waves' is not a YAML sequence." << std::endl;
		return false;
	}

	_nrWaves = configInputWaves.size();
	if(_debug) {
		printDebug << "going to read information of " << _nrWaves << " waves to be used in the fit." << std::endl;
	}

	for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
		const YAML::Node& configInputWave = configInputWaves[idxWave];

		std::map<std::string, YamlCppUtils::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("name", YamlCppUtils::TypeString);
		if(not checkIfAllVariablesAreThere(configInputWave, mandatoryArguments)) {
			printErr << "'waves' entry at index " << idxWave << " does not contain all required variables." << std::endl;
			return false;
		}

		const std::string name = configInputWave["name"].as<std::string>();

		double massLower = -1.;
		if(configInputWave["massLower"]) {
			if(checkVariableType(configInputWave["massLower"], YamlCppUtils::TypeFloat)) {
				massLower = configInputWave["massLower"].as<double>();
			} else {
				printErr << "variable 'massLower' of 'waves' entry at index " << idxWave << " is not a floating point number." << std::endl;
				return false;
			}
		}
		double massUpper = -1.;
		if(configInputWave["massUpper"]) {
			if(checkVariableType(configInputWave["massUpper"], YamlCppUtils::TypeFloat)) {
				massUpper = configInputWave["massUpper"].as<double>();
			} else {
				printErr << "variable 'massUpper' of 'waves' entry at index " << idxWave << " is not a floating point number." << std::endl;
				return false;
			}
		}

		std::vector<std::string> alternativeNames;
		if(configInputWave["alternativeNames"]) {
			if(checkVariableType(configInputWave["alternativeNames"], YamlCppUtils::TypeSequence)) {
				for(size_t idxAlt = 0; idxAlt < configInputWave["alternativeNames"].size(); ++idxAlt) {
					if(not checkVariableType(configInputWave["alternativeNames"][idxAlt], YamlCppUtils::TypeString)) {
						printErr << "element " << idxAlt << " of variable 'alternativeNames' of 'waves' entry at index " << idxWave << " is not a string." << std::endl;
						return false;
					}
					const std::string alternativeName = configInputWave["alternativeNames"][idxAlt].as<std::string>();

					// check that the alternative name does not yet exist
					if(alternativeName == name) {
						printErr << "alternative name '" << alternativeName << "' is equal to name of wave '" << name << "'." << std::endl;
						return false;
					}
					if(find(alternativeNames.begin(), alternativeNames.end(), alternativeName) != alternativeNames.end()) {
						printErr << "alternative name '" << alternativeName << "' of wave '" << name << "' defined twice." << std::endl;
						return false;
					}
					if(find(_waveNames.begin(), _waveNames.end(), alternativeName) != _waveNames.end()) {
						printErr << "alternative name '" << alternativeName << "' of wave '" << name << "' already defined as separate wave." << std::endl;
						return false;
					}
					for(size_t i = 0; i < _waveNameAlternatives.size(); ++i) {
						if(find(_waveNameAlternatives[i].begin(), _waveNameAlternatives[i].end(), alternativeName) != _waveNameAlternatives[i].end()) {
							printErr << "alternative name '" << alternativeName << "' of wave '" << name << "' already defined as alternative name of wave '" << _waveNames[i] << "'." << std::endl;
							return false;
						}
					}

					alternativeNames.push_back(alternativeName);
				}
			} else {
				printErr << "variable 'alternativeNames' of 'waves' entry at index " << idxWave << " is not a sequence." << std::endl;
				return false;
			}
		}

		// check that wave does not yet exist
		if(find(_waveNames.begin(), _waveNames.end(), name) != _waveNames.end()) {
			printErr << "wave '" << name << "' defined twice." << std::endl;
			return false;
		}
		for(size_t i = 0; i < _waveNameAlternatives.size(); ++i) {
			if(find(_waveNameAlternatives[i].begin(), _waveNameAlternatives[i].end(), name) != _waveNameAlternatives[i].end()) {
				printErr << "wave '" << name << "' already defined as alternative name of wave '" << _waveNames[i] << "'." << std::endl;
				return false;
			}
		}

		_waveNames.push_back(name);
		_waveNameAlternatives.push_back(alternativeNames);
		_waveIndices[name] = _waveNames.size() - 1;
		for(size_t idxAlt = 0; idxAlt < alternativeNames.size(); ++idxAlt) {
			_waveIndices[alternativeNames[idxAlt]] = _waveNames.size() - 1;
		}
		_waveMassLimits.push_back(std::make_pair(massLower, massUpper));

		if(_debug) {
			std::ostringstream alternatives;
			if(alternativeNames.size() > 0) {
				alternatives << ", alternative names: ";
				for(size_t idxAlt = 0; idxAlt < alternativeNames.size(); ++idxAlt) {
					if(idxAlt != 0) {
						alternatives << ", ";
					}
					alternatives << "'" << alternativeNames[idxAlt] << "'";
				}
			}

			printDebug << idxWave << ": " << name << " (mass range: " << massLower << "-" << massUpper << " GeV/c^2, index: " << _waveIndices[name] << alternatives.str() << ")" << std::endl;
		}
	}

	if(_debug) {
		printDebug << "read " << _nrWaves << " in total." << std::endl;
	}

	std::ostringstream output;
	for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
		output << "    " << _waveNames[idxWave] << std::endl;
	}
	printInfo << _nrWaves << " waves to be used in fit:" << std::endl
	          << output.str();

	return true;
}


bool
rpwa::resonanceFit::massDepFit::readConfigInputFreeParameters(const YAML::Node& configInputFreeParameters)
{
	if(not configInputFreeParameters) {
		printErr << "'configInputFreeParameters' is not a valid YAML node." << std::endl;
		return false;
	}
	if(not configInputFreeParameters.IsSequence()) {
		printErr << "'freeparameters' is not a YAML sequence." << std::endl;
		return false;
	}

	_freeParameters.clear();

	const size_t nrItems = configInputFreeParameters.size();
	if(nrItems == 0) {
		printErr << "'freeparameters' is an empty sequence, when defined it must at least contain one entry." << std::endl;
		return false;
	}

	if(_debug) {
		printDebug << "going to extract " << nrItems << " items from 'freeparameters'." << std::endl;
	}

	for(size_t idxItem=0; idxItem<nrItems; ++idxItem) {
		if(not checkVariableType(configInputFreeParameters[idxItem], YamlCppUtils::TypeString)) {
			printErr << "'freeparameters' entry at index " << idxItem << " is not a string." << std::endl;
			return false;
		}

		const std::string name = configInputFreeParameters[idxItem].as<std::string>();
		if(_debug) {
			printDebug << idxItem << ": '" << name << "'." << std::endl;
		}
		_freeParameters.push_back(name);
	}

	return true;
}


bool
rpwa::resonanceFit::massDepFit::readConfigModel(const YAML::Node& configModel,
                                                const rpwa::resonanceFit::modelPtr& fitModel,
                                                rpwa::resonanceFit::parameters& fitParameters,
                                                rpwa::resonanceFit::parameters& fitParametersError,
                                                const std::vector<size_t>& nrMassBins,
                                                const boost::multi_array<double, 2>& massBinCenters,
                                                const bool useBranchings)
{
	if(not configModel) {
		printErr << "'configModel' is not a valid YAML node." << std::endl;
		return false;
	}
	if(not configModel.IsMap()) {
		printErr << "'model' is not a YAML map." << std::endl;
		return false;
	}

	if(_debug) {
		printDebug << "reading fit model from configuration file." << std::endl;
	}

	// get information about anchor wave
	const YAML::Node& configAnchorWave = configModel["anchorwave"];
	if(not configAnchorWave) {
		printErr << "'anchorwave' does not exist in 'model'." << std::endl;
		return false;
	}
	if(not readConfigModelAnchorWave(configAnchorWave)) {
		printErr << "error while reading 'anchorwave' in 'model'." << std::endl;
		return false;
	}

	// read information for the individual components
	const YAML::Node& configComponents = configModel["components"];
	if(not configComponents) {
		printErr << "'components' does not exist in 'model'." << std::endl;
		return false;
	}
	if(not readConfigModelComponents(configComponents,
	                                 fitModel,
	                                 fitParameters,
	                                 fitParametersError,
	                                 nrMassBins,
	                                 massBinCenters,
	                                 useBranchings)) {
		printErr << "error while reading 'components' in 'model'." << std::endl;
		return false;
	}

	// get information for creating the final-state mass-dependence
	const YAML::Node& configFsmd = configModel["finalStateMassDependence"];
	if(configFsmd) {
		if(not readConfigModelFsmd(configFsmd, fitModel, fitParameters, fitParametersError)) {
			printErr << "error while reading 'finalStateMassDependence' in 'model'." << std::endl;
			return false;
		}
	} else {
		printInfo << "not using final-state mass-dependence." << std::endl;
	}

	return true;
}


bool
rpwa::resonanceFit::massDepFit::readConfigModelAnchorWave(const YAML::Node& configAnchorWave)
{
	if(not configAnchorWave) {
		printErr << "'configInputAnchorWave' is not a valid YAML node." << std::endl;
		return false;
	}
	if(not configAnchorWave.IsMap()) {
		printErr << "'anchorwave' is not a YAML map." << std::endl;
		return false;
	}

	std::map<std::string, YamlCppUtils::Type> mandatoryArguments;
	boost::assign::insert(mandatoryArguments)
	                     ("name", YamlCppUtils::TypeString)
	                     ("resonance", YamlCppUtils::TypeString);
	if(not checkIfAllVariablesAreThere(configAnchorWave, mandatoryArguments)) {
		printErr << "'anchorwave' does not contain all required variables." << std::endl;
		return false;
	}

	_anchorWaveName = configAnchorWave["name"].as<std::string>();
	_anchorComponentName = configAnchorWave["resonance"].as<std::string>();

	return true;
}


bool
rpwa::resonanceFit::massDepFit::readConfigModelComponents(const YAML::Node& configComponents,
                                                          const rpwa::resonanceFit::modelPtr& fitModel,
                                                          rpwa::resonanceFit::parameters& fitParameters,
                                                          rpwa::resonanceFit::parameters& fitParametersError,
                                                          const std::vector<size_t>& nrMassBins,
                                                          const boost::multi_array<double, 2>& massBinCenters,
                                                          const bool useBranchings) const
{
	if(not configComponents) {
		printErr << "'configComponents' is not a valid YAML node." << std::endl;
		return false;
	}
	if(not configComponents.IsSequence()) {
		printErr << "'components' is not a YAML sequence." << std::endl;
		return false;
	}

	const size_t nrComponents = configComponents.size();

	if(_debug) {
		printDebug << "reading " << nrComponents << " components from configuration file." << std::endl;
	}

	for(size_t idxComponent=0; idxComponent<nrComponents; ++idxComponent) {
		const YAML::Node& configComponent = configComponents[idxComponent];

		std::map<std::string, YamlCppUtils::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("name", YamlCppUtils::TypeString);
		if(not checkIfAllVariablesAreThere(configComponent, mandatoryArguments)) {
			printErr << "'components' entry at index " << idxComponent << " does not contain all required variables." << std::endl;
			return false;
		}

		const std::string name = configComponent["name"].as<std::string>();

		for(size_t idx = 0; idx < fitModel->getNrComponents(); ++idx) {
			if(fitModel->getComponent(idx)->getName() == name) {
				printErr << "component '" << name << "' defined twice." << std::endl;
				return false;
			}
		}

		std::string type;
		if(configComponent["type"]) {
			if(not checkVariableType(configComponent["type"], YamlCppUtils::TypeString)) {
				printErr << "component '" << name << "' has a type that is not a string." << std::endl;
				return false;
			}
			type = configComponent["type"].as<std::string>();
		} else {
			if(_debug) {
				printDebug << "component '" << name << "' has no type, use 'fixedWidthBreitWigner'." << std::endl;
			}
			type = "fixedWidthBreitWigner";
		}

		if(_debug) {
			printDebug << "found component '" << name << "' with type '" << type << "'." << std::endl;
		}

		rpwa::resonanceFit::componentPtr component;
		if(type == "fixedWidthBreitWigner") {
			component = std::make_shared<rpwa::resonanceFit::fixedWidthBreitWigner>(fitModel->getNrComponents(), name);
		} else if(type == "dynamicWidthBreitWigner") {
			component = std::make_shared<rpwa::resonanceFit::dynamicWidthBreitWigner>(fitModel->getNrComponents(), name);
		} else if(type == "integralWidthBreitWigner") {
			component = std::make_shared<rpwa::resonanceFit::integralWidthBreitWigner>(fitModel->getNrComponents(), name);
		} else if(type == "constantBackground") {
			component = std::make_shared<rpwa::resonanceFit::constantBackground>(fitModel->getNrComponents(), name);
		} else if(type == "exponentialBackground") {
			component = std::make_shared<rpwa::resonanceFit::exponentialBackground>(fitModel->getNrComponents(), name);
		} else if(type == "tPrimeDependentBackground") {
			component = std::make_shared<rpwa::resonanceFit::tPrimeDependentBackground>(fitModel->getNrComponents(), name);
			std::dynamic_pointer_cast<tPrimeDependentBackground>(component)->setTPrimeMeans(_tPrimeMeans);
		} else if(type == "exponentialBackgroundIntegral") {
			component = std::make_shared<rpwa::resonanceFit::exponentialBackgroundIntegral>(fitModel->getNrComponents(), name);
		} else if(type == "tPrimeDependentBackgroundIntegral") {
			component = std::make_shared<rpwa::resonanceFit::tPrimeDependentBackgroundIntegral>(fitModel->getNrComponents(), name);
			std::dynamic_pointer_cast<tPrimeDependentBackgroundIntegral>(component)->setTPrimeMeans(_tPrimeMeans);
		} else {
			printErr << "unknown type '" << type << "' for component '" << name << "'." << std::endl;
			return false;
		}

		if(not component->init(configComponent,
		                       fitParameters,
		                       fitParametersError,
		                       nrMassBins,
		                       massBinCenters,
		                       _waveIndices,
		                       _waveBins,
		                       _inPhaseSpaceIntegrals,
		                       useBranchings,
		                       _debug)) {
			printErr << "error while initializing component '" << name << "' of type '" << type << "'." << std::endl;
			return false;
		}

		fitModel->add(component);
	}

	std::ostringstream output;
	for(size_t idxComponent = 0; idxComponent < fitModel->getNrComponents(); ++idxComponent) {
		output << "    " << fitModel->getComponent(idxComponent)->getName() << std::endl;
	}
	printInfo << "fitting " << fitModel->getNrComponents() << " components to the data:" << std::endl
	          << output.str();

	return true;
}


bool
rpwa::resonanceFit::massDepFit::readConfigModelFsmd(const YAML::Node& configFsmd,
                                                    const rpwa::resonanceFit::modelPtr& fitModel,
                                                    rpwa::resonanceFit::parameters& fitParameters,
                                                    rpwa::resonanceFit::parameters& fitParametersError) const
{
	if(not configFsmd) {
		printErr << "'configFsmd' is not a valid YAML node." << std::endl;
		return false;
	}
	if(not configFsmd.IsMap() and not configFsmd.IsSequence()) {
		printErr << "'finalStateMassDependence' is not a YAML map or sequence." << std::endl;
		return false;
	}

	if(_debug) {
		printDebug << "reading final-state mass-dependence from configuration file." << std::endl;
	}

	rpwa::resonanceFit::fsmdPtr fsmd(new rpwa::resonanceFit::fsmd(fitModel->getNrComponents()));
	if(not fsmd->init(configFsmd, fitParameters, fitParametersError, _nrBins, _sameMassBinning, _debug)) {
		printErr << "error while initializing final-state mass-dependence." << std::endl;
		return false;
	}
	fitModel->setFsmd(fsmd);

	printInfo << "using final-state mass-dependence as defined in the configuration file." << std::endl;

	return true;
}


bool
rpwa::resonanceFit::massDepFit::init(const rpwa::resonanceFit::dataConstPtr& fitData,
                                     const rpwa::resonanceFit::modelPtr& fitModel,
                                     const rpwa::resonanceFit::functionPtr& fitFunction)
{
	if(not fitModel->init(_nrBins,
	                      _waveNames,
	                      _waveNameAlternatives,
	                      _anchorWaveName,
	                      _anchorComponentName)) {
		printErr << "error while initializing the fit model." << std::endl;
		return false;
	}

	if(not fitFunction->init(fitData,
	                         fitModel,
	                         _wavePairMassBinLimits)) {
		printErr << "error while initializing the function to minimize." << std::endl;
		return false;
	}

	return true;
}


bool
rpwa::resonanceFit::massDepFit::writeConfig(std::ostream& output,
                                            const rpwa::resonanceFit::modelConstPtr& fitModel,
                                            const rpwa::resonanceFit::parameters& fitParameters,
                                            const rpwa::resonanceFit::parameters& fitParametersError,
                                            const std::map<std::string, double>& fitQuality) const
{
	if(_debug) {
		printDebug << "writing configuration file." << std::endl;
	}

	YAML::Emitter yamlOutput(output);
	yamlOutput << YAML::BeginMap;

	yamlOutput << YAML::Key << "fitquality";
	yamlOutput << YAML::Value;
	if(not writeConfigFitquality(yamlOutput, fitQuality)) {
		printErr << "error while writing 'fitquality' to result file." << std::endl;
		return false;
	}

	yamlOutput << YAML::Key << "input";
	yamlOutput << YAML::Value;
	if(not writeConfigInput(yamlOutput)) {
		printErr << "error while writing 'input' to result file." << std::endl;
		return false;
	}

	yamlOutput << YAML::Key << "model";
	yamlOutput << YAML::Value;
	if(not writeConfigModel(yamlOutput, fitModel, fitParameters, fitParametersError)) {
		printErr << "error while writing 'model' to result file." << std::endl;
		return false;
	}

	yamlOutput << YAML::EndMap;

	// newline at end-of-file
	output << std::endl;

	return true;
}


bool
rpwa::resonanceFit::massDepFit::writeConfigFitquality(YAML::Emitter& yamlOutput,
                                                      const std::map<std::string, double>& fitQuality) const
{
	if(_debug) {
		printDebug << "writing 'fitquality'." << std::endl;
	}

	yamlOutput << YAML::BeginMap;

	for(std::map<std::string, double>::const_iterator it = fitQuality.begin(); it != fitQuality.end(); ++it) {
		yamlOutput << YAML::Key << it->first;
		yamlOutput << YAML::Value << it->second;
	}

	yamlOutput << YAML::EndMap;

	return true;
}


bool
rpwa::resonanceFit::massDepFit::writeConfigInput(YAML::Emitter& yamlOutput) const
{
	if(_debug) {
		printDebug << "writing 'input'." << std::endl;
	}

	yamlOutput << YAML::BeginMap;

	yamlOutput << YAML::Key << "fitresults";
	yamlOutput << YAML::Value;
	if(not writeConfigInputFitResults(yamlOutput)) {
		printErr << "error while writing 'fitresults' to result file." << std::endl;
		return false;
	}

	yamlOutput << YAML::Key << "waves";
	yamlOutput << YAML::Value;
	if(not writeConfigInputWaves(yamlOutput)) {
		printErr << "error while writing 'waves' to result file." << std::endl;
		return false;
	}

	yamlOutput << YAML::Key << "freeparameters";
	yamlOutput << YAML::Value;
	if(not writeConfigInputFreeParameters(yamlOutput)) {
		printErr << "error while writing 'freeparameters' to result file." << std::endl;
		return false;
	}

	yamlOutput << YAML::EndMap;

	return true;
}


bool
rpwa::resonanceFit::massDepFit::writeConfigInputFitResults(YAML::Emitter& yamlOutput) const
{
	if(_debug) {
		printDebug << "writing 'fitresults'." << std::endl;
	}

	yamlOutput << YAML::BeginSeq;

	for (size_t idxBin=0; idxBin<_nrBins; ++idxBin) {
		yamlOutput << YAML::BeginMap;

		yamlOutput << YAML::Key << "name";
		yamlOutput << YAML::Value << _inFileName[idxBin];

		yamlOutput << YAML::Key << "tPrimeMean";
		yamlOutput << YAML::Value << _tPrimeMeans[idxBin];

		if(_rescaleErrors[idxBin] != 1.) {
			yamlOutput << YAML::Key << "rescaleErrors";
			yamlOutput << YAML::Value << _rescaleErrors[idxBin];
		}

		if(_nrSystematics[idxBin] > 0) {
			yamlOutput << YAML::Key << "systematics";
			yamlOutput << YAML::Value;
			if(not writeConfigInputFitResultSystematics(yamlOutput, idxBin)) {
				printErr << "error while writing 'systematics' to result file." << std::endl;
				return false;
			}
		}

		yamlOutput << YAML::EndMap;
	}

	yamlOutput << YAML::EndSeq;

	return true;
}


bool
rpwa::resonanceFit::massDepFit::writeConfigInputFitResultSystematics(YAML::Emitter& yamlOutput,
                                                                     const size_t idxBin) const
{
	if(_debug) {
		printDebug << "writing 'systematics' for bin " << idxBin << "." << std::endl;
	}

	yamlOutput << _sysFileNames[idxBin];

	return true;
}


bool
rpwa::resonanceFit::massDepFit::writeConfigInputWaves(YAML::Emitter& yamlOutput) const
{
	if(_debug) {
		printDebug << "writing 'waves'." << std::endl;
	}

	yamlOutput << YAML::BeginSeq;

	for (size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
		yamlOutput << YAML::BeginMap;

		yamlOutput << YAML::Key << "name";
		yamlOutput << YAML::Value << _waveNames[idxWave];

		if (_waveMassLimits[idxWave].first >= 0) {
			yamlOutput << YAML::Key << "massLower";
			yamlOutput << YAML::Value << _waveMassLimits[idxWave].first;
		}
		if (_waveMassLimits[idxWave].second >= 0) {
			yamlOutput << YAML::Key << "massUpper";
			yamlOutput << YAML::Value << _waveMassLimits[idxWave].second;
		}

		if(_waveNameAlternatives[idxWave].size() > 0) {
			yamlOutput << YAML::Key << "alternativeNames";
			yamlOutput << YAML::Value;

			yamlOutput << YAML::Flow;
			yamlOutput << YAML::BeginSeq;

			for(size_t idxAlt = 0; idxAlt < _waveNameAlternatives[idxWave].size(); ++idxAlt) {
				yamlOutput << _waveNameAlternatives[idxWave][idxAlt];
			}

			yamlOutput << YAML::EndSeq;
			yamlOutput << YAML::Block;
		}

		yamlOutput << YAML::EndMap;
	}

	yamlOutput << YAML::EndSeq;

	return true;
}


bool
rpwa::resonanceFit::massDepFit::writeConfigInputFreeParameters(YAML::Emitter& yamlOutput) const
{
	if(_debug) {
		printDebug << "writing 'freeparameters'." << std::endl;
	}

	yamlOutput << _freeParameters;

	return true;
}


bool
rpwa::resonanceFit::massDepFit::writeConfigModel(YAML::Emitter& yamlOutput,
                                                 const rpwa::resonanceFit::modelConstPtr& fitModel,
                                                 const rpwa::resonanceFit::parameters& fitParameters,
                                                 const rpwa::resonanceFit::parameters& fitParametersError) const
{
	if(_debug) {
		printDebug << "writing 'model'." << std::endl;
	}

	yamlOutput << YAML::BeginMap;

	yamlOutput << YAML::Key << "anchorwave";
	yamlOutput << YAML::Value;
	if(not writeConfigModelAnchorWave(yamlOutput)) {
		printErr << "error while writing 'anchorwave' to result file." << std::endl;
		return false;
	}

	yamlOutput << YAML::Key << "components";
	yamlOutput << YAML::Value;
	if(not writeConfigModelComponents(yamlOutput, fitModel, fitParameters, fitParametersError)) {
		printErr << "error while writing 'components' to result file." << std::endl;
		return false;
	}

	yamlOutput << YAML::Key << "finalStateMassDependence";
	yamlOutput << YAML::Value;
	if(fitModel->getFsmd()) {
		if(not writeConfigModelFsmd(yamlOutput, fitModel, fitParameters, fitParametersError)) {
			printErr << "error while writing 'finalStateMassDependence' to result file." << std::endl;
			return false;
		}
	}

	yamlOutput << YAML::EndMap;

	return true;
}


bool
rpwa::resonanceFit::massDepFit::writeConfigModelAnchorWave(YAML::Emitter& yamlOutput) const
{
	if(_debug) {
		printDebug << "writing 'anchorwave'." << std::endl;
	}

	yamlOutput << YAML::BeginMap;

	yamlOutput << YAML::Key << "name";
	yamlOutput << YAML::Value << _anchorWaveName;

	yamlOutput << YAML::Key << "resonance";
	yamlOutput << YAML::Value << _anchorComponentName;

	yamlOutput << YAML::EndMap;

	return true;
}


bool
rpwa::resonanceFit::massDepFit::writeConfigModelComponents(YAML::Emitter& yamlOutput,
                                                           const rpwa::resonanceFit::modelConstPtr& fitModel,
                                                           const rpwa::resonanceFit::parameters& fitParameters,
                                                           const rpwa::resonanceFit::parameters& fitParametersError) const
{
	if(_debug) {
		printDebug << "writing 'components'." << std::endl;
	}

	yamlOutput << YAML::BeginSeq;

	const size_t nrComponents = fitModel->getNrComponents();
	for(size_t idxComponent=0; idxComponent<nrComponents; ++idxComponent) {
		if(not fitModel->getComponent(idxComponent)->write(yamlOutput, fitParameters, fitParametersError, _debug)) {
			printErr << "error while writing component at index " << idxComponent << " to result file." << std::endl;
			return false;
		}
	}

	yamlOutput << YAML::EndSeq;

	return true;
}


bool
rpwa::resonanceFit::massDepFit::writeConfigModelFsmd(YAML::Emitter& yamlOutput,
                                                     const rpwa::resonanceFit::modelConstPtr& fitModel,
                                                     const rpwa::resonanceFit::parameters& fitParameters,
                                                     const rpwa::resonanceFit::parameters& fitParametersError) const
{
	if(_debug) {
		printDebug << "writing 'finalStateMassDependence'." << std::endl;
	}

	if(not fitModel->getFsmd()) {
		printErr << "writing final-state mass-dependence requested, but there is no final-state mass-dependence." << std::endl;
		return false;
	}

	if(not fitModel->getFsmd()->write(yamlOutput, fitParameters, fitParametersError, _debug)) {
		printErr << "error while writing final-state mass-dependence to result file." << std::endl;
		return false;
	}

	return true;
}


bool
rpwa::resonanceFit::massDepFit::readInFiles(std::vector<size_t>& nrMassBins,
                                            boost::multi_array<double, 2>& massBinCenters,
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
                                            const std::string& valTreeName,
                                            const std::string& valBranchName)
{
	for(size_t idxBin = 0; idxBin < _nrBins; ++idxBin) {
		boost::multi_array<std::complex<double>, 2> tempProductionAmplitudes;
		boost::multi_array<TMatrixT<double>, 1> tempProductionAmplitudesCovariance;
		boost::multi_array<std::complex<double>, 3> tempSpinDensityMatrices;
		boost::multi_array<TMatrixT<double>, 1> tempSpinDensityMatricesCovariance;
		boost::multi_array<std::pair<double, double>, 2> tempPlottingIntensities;
		boost::multi_array<std::pair<double, double>, 3> tempPlottingSpinDensityMatrixElementsReal;
		boost::multi_array<std::pair<double, double>, 3> tempPlottingSpinDensityMatrixElementsImag;
		boost::multi_array<std::pair<double, double>, 3> tempPlottingPhases;

		if(not readInFile(idxBin,
		                  nrMassBins,
		                  massBinCenters,
		                  tempProductionAmplitudes,
		                  tempProductionAmplitudesCovariance,
		                  tempSpinDensityMatrices,
		                  tempSpinDensityMatricesCovariance,
		                  tempPlottingIntensities,
		                  tempPlottingSpinDensityMatrixElementsReal,
		                  tempPlottingSpinDensityMatrixElementsImag,
		                  tempPlottingPhases,
		                  valTreeName,
		                  valBranchName)) {
			printErr << "error while reading file entry " << idxBin << "." << std::endl;
			return false;
		}

		adjustSizeAndSet(productionAmplitudes, idxBin, tempProductionAmplitudes);
		adjustSizeAndSet(productionAmplitudesCovariance, idxBin, tempProductionAmplitudesCovariance);
		adjustSizeAndSet(spinDensityMatrices, idxBin, tempSpinDensityMatrices);
		adjustSizeAndSet(spinDensityMatricesCovariance, idxBin, tempSpinDensityMatricesCovariance);
		adjustSizeAndSet(plottingIntensities, idxBin, tempPlottingIntensities);
		adjustSizeAndSet(plottingSpinDensityMatrixElementsReal, idxBin, tempPlottingSpinDensityMatrixElementsReal);
		adjustSizeAndSet(plottingSpinDensityMatrixElementsImag, idxBin, tempPlottingSpinDensityMatrixElementsImag);
		adjustSizeAndSet(plottingPhases, idxBin, tempPlottingPhases);

		// extract information for systematic errors
		// initialize with real fit result
		boost::multi_array<std::pair<double, double>, 2> tempSysPlottingIntensities(std::vector<size_t>(tempPlottingIntensities.shape(), tempPlottingIntensities.shape()+tempPlottingIntensities.num_dimensions()));
		boost::multi_array<std::pair<double, double>, 3> tempSysPlottingSpinDensityMatrixElementsReal(std::vector<size_t>(tempPlottingSpinDensityMatrixElementsReal.shape(), tempPlottingSpinDensityMatrixElementsReal.shape()+tempPlottingSpinDensityMatrixElementsReal.num_dimensions()));
		boost::multi_array<std::pair<double, double>, 3> tempSysPlottingSpinDensityMatrixElementsImag(std::vector<size_t>(tempPlottingSpinDensityMatrixElementsImag.shape(), tempPlottingSpinDensityMatrixElementsImag.shape()+tempPlottingSpinDensityMatrixElementsImag.num_dimensions()));
		boost::multi_array<std::pair<double, double>, 3> tempSysPlottingPhases(std::vector<size_t>(tempPlottingPhases.shape(), tempPlottingPhases.shape()+tempPlottingPhases.num_dimensions()));

		for(size_t idxMass = 0; idxMass < nrMassBins[idxBin]; ++idxMass) {
			for(size_t idxWave = 0; idxWave < _nrWaves; ++idxWave) {
				tempSysPlottingIntensities[idxMass][idxWave] = std::make_pair(tempPlottingIntensities[idxMass][idxWave].first,
				                                                              tempPlottingIntensities[idxMass][idxWave].first);

				for(size_t jdxWave = 0; jdxWave < _nrWaves; ++jdxWave) {
					tempSysPlottingSpinDensityMatrixElementsReal[idxMass][idxWave][jdxWave] = std::make_pair(tempPlottingSpinDensityMatrixElementsReal[idxMass][idxWave][jdxWave].first,
					                                                                                         tempPlottingSpinDensityMatrixElementsReal[idxMass][idxWave][jdxWave].first);
					tempSysPlottingSpinDensityMatrixElementsImag[idxMass][idxWave][jdxWave] = std::make_pair(tempPlottingSpinDensityMatrixElementsImag[idxMass][idxWave][jdxWave].first,
					                                                                                         tempPlottingSpinDensityMatrixElementsImag[idxMass][idxWave][jdxWave].first);
					tempSysPlottingPhases[idxMass][idxWave][jdxWave] = std::make_pair(tempPlottingPhases[idxMass][idxWave][jdxWave].first,
					                                                                  tempPlottingPhases[idxMass][idxWave][jdxWave].first);
				}
			}
		}

		if(_nrSystematics[idxBin] > 0) {
			if(not readSystematicsFiles(idxBin,
			                            nrMassBins[idxBin],
			                            massBinCenters[idxBin],
			                            tempPlottingPhases,
			                            tempSysPlottingIntensities,
			                            tempSysPlottingSpinDensityMatrixElementsReal,
			                            tempSysPlottingSpinDensityMatrixElementsImag,
			                            tempSysPlottingPhases,
			                            valTreeName,
			                            valBranchName)) {
				printErr << "error while reading fit results for systematic errors in bin " << idxBin << "." << std::endl;
				return false;
			}
		}

		adjustSizeAndSet(sysPlottingIntensities, idxBin, tempSysPlottingIntensities);
		adjustSizeAndSet(sysPlottingSpinDensityMatrixElementsReal, idxBin, tempSysPlottingSpinDensityMatrixElementsReal);
		adjustSizeAndSet(sysPlottingSpinDensityMatrixElementsImag, idxBin, tempSysPlottingSpinDensityMatrixElementsImag);
		adjustSizeAndSet(sysPlottingPhases, idxBin, tempSysPlottingPhases);
	}

	return true;
}


bool
rpwa::resonanceFit::massDepFit::readInFile(const size_t idxBin,
                                           std::vector<size_t>& nrMassBins,
                                           boost::multi_array<double, 2>& massBinCenters,
                                           boost::multi_array<std::complex<double>, 2>& productionAmplitudes,
                                           boost::multi_array<TMatrixT<double>, 1>& productionAmplitudesCovariance,
                                           boost::multi_array<std::complex<double>, 3>& spinDensityMatrices,
                                           boost::multi_array<TMatrixT<double>, 1>& spinDensityMatricesCovariance,
                                           boost::multi_array<std::pair<double, double>, 2>& plottingIntensities,
                                           boost::multi_array<std::pair<double, double>, 3>& plottingSpinDensityMatrixElementsReal,
                                           boost::multi_array<std::pair<double, double>, 3>& plottingSpinDensityMatrixElementsImag,
                                           boost::multi_array<std::pair<double, double>, 3>& plottingPhases,
                                           const std::string& valTreeName,
                                           const std::string& valBranchName)
{
	if(_debug) {
		printDebug << "reading fit result from file '" << _inFileName[idxBin] << "'." << std::endl;
	}

	TFile* inFile = TFile::Open(_inFileName[idxBin].c_str());
	if(not inFile) {
		printErr << "input file '" << _inFileName[idxBin] << "' not found."<< std::endl;
		return false;
	}
	if(inFile->IsZombie()) {
		printErr << "error while reading input file '" << _inFileName[idxBin] << "'."<< std::endl;
		delete inFile;
		return false;
	}

	if(_debug) {
		printDebug << "searching for tree '" << valTreeName << "' in file '" << _inFileName[idxBin] << "'." << std::endl;
	}

	TTree* inTree;
	inFile->GetObject(valTreeName.c_str(), inTree);
	if(not inTree) {
		printErr << "input tree '" << valTreeName << "' not found in input file '" << _inFileName[idxBin] << "'."<< std::endl;
		delete inFile;
		return false;
	}

	if(_debug) {
		printDebug << "searching for branch '" << valBranchName << "' in tree '" << valTreeName << "'." << std::endl;
	}

	fitResult* inFit = NULL;
	if(inTree->SetBranchAddress(valBranchName.c_str(), &inFit)) {
		printErr << "branch '" << valBranchName << "' not found in input tree '" << valTreeName << "'." << std::endl;
		delete inFile;
		return false;
	}

	size_t tempNrMassBins;
	boost::multi_array<double, 1> tempMassBinCenters;
	if(not readFitResultMassBins(inTree,
	                             inFit,
	                             tempNrMassBins,
	                             tempMassBinCenters)) {
		printErr << "could not extract mass bins from fit result tree in '" << _inFileName[idxBin] << "'." << std::endl;
		delete inFile;
		return false;
	}

	adjustSizeAndSet(nrMassBins, idxBin, tempNrMassBins);
	adjustSizeAndSet(massBinCenters, idxBin, tempMassBinCenters);

	if(_maxMassBins < tempNrMassBins) {
		_maxMassBins = tempNrMassBins;

		// resize all array to store the information
		_inPhaseSpaceIntegrals.resize(boost::extents[_nrBins][_maxMassBins][_nrWaves]);
	}

	bool readMapping(false);
	std::vector<Long64_t> inMapping;
	if(_sameMassBinning and idxBin > 0) {
		if(checkFitResultMassBins(inTree,
		                          inFit,
		                          nrMassBins[idxBin-1],
		                          massBinCenters[idxBin-1],
		                          inMapping)) {
			if(_debug) {
				printDebug << "bin " << idxBin << " has the same mass binning as bin " << idxBin-1 << "." << std::endl;
			}

			readMapping = true;
		} else {
			printInfo << "bin " << idxBin << " does not have the same mass binning as bin " << idxBin-1 << ", use individual mass binning for each bin." << std::endl;

			_sameMassBinning = false;
		}
	}
	if(not readMapping) {
		if(not checkFitResultMassBins(inTree,
		                              inFit,
		                              nrMassBins[idxBin],
		                              massBinCenters[idxBin],
		                              inMapping)) {
			printErr << "error while checking and mapping mass bins from fit result tree in '" << _inFileName[idxBin] << "'." << std::endl;
			delete inFile;
			return false;
		}
	}

	std::vector<std::string> waveNames;
	if(not readFitResultMatrices(inTree,
	                             inFit,
	                             inMapping,
	                             _rescaleErrors[idxBin],
	                             waveNames,
	                             productionAmplitudes,
	                             productionAmplitudesCovariance,
	                             spinDensityMatrices,
	                             spinDensityMatricesCovariance,
	                             plottingIntensities,
	                             plottingSpinDensityMatrixElementsReal,
	                             plottingSpinDensityMatrixElementsImag,
	                             plottingPhases)) {
		printErr << "error while reading spin-density matrix from fit result tree in '" << _inFileName[idxBin] << "'." << std::endl;
		delete inFile;
		return false;
	}
	for(std::vector<std::string>::const_iterator it = waveNames.begin(); it != waveNames.end(); ++it) {
		_waveBins[*it].push_back(idxBin);
	}

	boost::multi_array<double, 2> tempPhaseSpaceIntegrals;
	if(not readFitResultIntegrals(inTree, inFit, inMapping, waveNames, tempPhaseSpaceIntegrals)) {
		printErr << "error while reading phase-space integrals from fit result tree in '" << _inFileName[idxBin] << "'." << std::endl;
		delete inFile;
		return false;
	}
	_inPhaseSpaceIntegrals[idxBin] = tempPhaseSpaceIntegrals;

	delete inFile;
	return true;
}


bool
rpwa::resonanceFit::massDepFit::readSystematicsFiles(const size_t idxBin,
                                                     const size_t nrMassBins,
                                                     const boost::multi_array<double, 1>& massBinCenters,
                                                     const boost::multi_array<std::pair<double, double>, 3>& plottingPhases,
                                                     boost::multi_array<std::pair<double, double>, 2>& sysPlottingIntensities,
                                                     boost::multi_array<std::pair<double, double>, 3>& sysPlottingSpinDensityMatrixElementsReal,
                                                     boost::multi_array<std::pair<double, double>, 3>& sysPlottingSpinDensityMatrixElementsImag,
                                                     boost::multi_array<std::pair<double, double>, 3>& sysPlottingPhases,
                                                     const std::string& valTreeName,
                                                     const std::string& valBranchName)
{
	if(_nrSystematics[idxBin] == 0) {
		return true;
	}

	if(_debug) {
		printDebug << "reading fit results for systematic errors for bin " << idxBin << " from " << _nrSystematics[idxBin] << " files." << std::endl;
	}

	for(size_t idxSystematics = 0; idxSystematics < _nrSystematics[idxBin]; ++idxSystematics) {
		if(not readSystematicsFile(idxBin,
		                           idxSystematics,
		                           nrMassBins,
		                           massBinCenters,
		                           plottingPhases,
		                           sysPlottingIntensities,
		                           sysPlottingSpinDensityMatrixElementsReal,
		                           sysPlottingSpinDensityMatrixElementsImag,
		                           sysPlottingPhases,
		                           valTreeName,
		                           valBranchName)) {
			printErr << "error while reading fit results for systematic errors." << std::endl;
			return false;
		}
	}

	return true;
}


bool
rpwa::resonanceFit::massDepFit::readSystematicsFile(const size_t idxBin,
                                                    const size_t idxSystematics,
                                                    const size_t nrMassBins,
                                                    const boost::multi_array<double, 1>& massBinCenters,
                                                    const boost::multi_array<std::pair<double, double>, 3>& plottingPhases,
                                                    boost::multi_array<std::pair<double, double>, 2>& sysPlottingIntensities,
                                                    boost::multi_array<std::pair<double, double>, 3>& sysPlottingSpinDensityMatrixElementsReal,
                                                    boost::multi_array<std::pair<double, double>, 3>& sysPlottingSpinDensityMatrixElementsImag,
                                                    boost::multi_array<std::pair<double, double>, 3>& sysPlottingPhases,
                                                    const std::string& valTreeName,
                                                    const std::string& valBranchName)
{
	if(_debug) {
		printDebug << "reading fit result for systematics for bin " << idxBin << " from file at index " << idxSystematics << ": '" << _sysFileNames[idxBin][idxSystematics] << "'." << std::endl;
	}

	TFile* sysFile = TFile::Open(_sysFileNames[idxBin][idxSystematics].c_str());
	if(not sysFile) {
		printErr << "input file '" << _sysFileNames[idxBin][idxSystematics] << "' not found."<< std::endl;
		return false;
	}
	if(sysFile->IsZombie()) {
		printErr << "error while reading input file '" << _sysFileNames[idxBin][idxSystematics] << "'."<< std::endl;
		delete sysFile;
		return false;
	}

	if(_debug) {
		printDebug << "searching for tree '" << valTreeName << "' in file '" << _sysFileNames[idxBin][idxSystematics] << "'." << std::endl;
	}

	TTree* sysTree;
	sysFile->GetObject(valTreeName.c_str(), sysTree);
	if(not sysTree) {
		printErr << "input tree '" << valTreeName << "' not found in input file '" << _sysFileNames[idxBin][idxSystematics] << "'."<< std::endl;
		delete sysFile;
		return false;
	}

	if(_debug) {
		printDebug << "searching for branch '" << valBranchName << "' in tree '" << valTreeName << "'." << std::endl;
	}

	fitResult* sysFit = NULL;
	if(sysTree->SetBranchAddress(valBranchName.c_str(), &sysFit)) {
		printErr << "branch '" << valBranchName << "' not found in input tree '" << valTreeName << "'." << std::endl;
		delete sysFile;
		return false;
	}

	std::vector<Long64_t> sysMapping;
	if(not checkFitResultMassBins(sysTree,
	                              sysFit,
	                              nrMassBins,
	                              massBinCenters,
	                              sysMapping)) {
		printErr << "error while checking and mapping mass bins from fit result tree in '" << _sysFileNames[idxBin][idxSystematics] << "'." << std::endl;
		delete sysFile;
		return false;
	}

	std::vector<std::string> waveNames;
	boost::multi_array<std::complex<double>, 2> tempProductionAmplitudes;
	boost::multi_array<TMatrixT<double>, 1> tempProductionAmplitudesCovariance;
	boost::multi_array<std::complex<double>, 3> tempSpinDensityMatrices;
	boost::multi_array<TMatrixT<double>, 1> tempSpinDensityCovarianceMatrices;
	boost::multi_array<std::pair<double, double>, 2> tempSysPlottingIntensities;
	boost::multi_array<std::pair<double, double>, 3> tempSysPlottingSpinDensityMatrixElementsReal;
	boost::multi_array<std::pair<double, double>, 3> tempSysPlottingSpinDensityMatrixElementsImag;
	boost::multi_array<std::pair<double, double>, 3> tempSysPlottingPhases;
	if(not readFitResultMatrices(sysTree,
	                             sysFit,
	                             sysMapping,
	                             _rescaleErrors[idxBin],
	                             waveNames,
	                             tempProductionAmplitudes,
	                             tempProductionAmplitudesCovariance,
	                             tempSpinDensityMatrices,
	                             tempSpinDensityCovarianceMatrices,
	                             tempSysPlottingIntensities,
	                             tempSysPlottingSpinDensityMatrixElementsReal,
	                             tempSysPlottingSpinDensityMatrixElementsImag,
	                             tempSysPlottingPhases)) {
		printErr << "error while reading spin-density matrix from fit result tree in '" << _sysFileNames[idxBin][idxSystematics] << "'." << std::endl;
		delete sysFile;
		return false;
	}

	for(size_t idxMass = 0; idxMass < nrMassBins; ++idxMass) {
		for(size_t idxWave = 0; idxWave < _nrWaves; ++idxWave) {
			sysPlottingIntensities[idxMass][idxWave].first = std::min(sysPlottingIntensities[idxMass][idxWave].first,
			                                                          tempSysPlottingIntensities[idxMass][idxWave].first);
			sysPlottingIntensities[idxMass][idxWave].second = std::max(sysPlottingIntensities[idxMass][idxWave].second,
			                                                           tempSysPlottingIntensities[idxMass][idxWave].first);

			for(size_t jdxWave = 0; jdxWave < _nrWaves; ++jdxWave) {
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

	delete sysFile;
	return true;
}


bool
rpwa::resonanceFit::massDepFit::checkFitResultMassBins(TTree* tree,
                                                       rpwa::fitResult* fit,
                                                       const size_t nrMassBins,
                                                       const boost::multi_array<double, 1>& massBinCenters,
                                                       std::vector<Long64_t>& mapping) const
{
	if(not tree or not fit) {
		printErr << "'tree' or 'fit' is not a pointer to a valid object." << std::endl;
		return false;
	}

	// reset mapping
	mapping.assign(nrMassBins, std::numeric_limits<Long64_t>::max());

	// extract data from tree
	const Long64_t nrEntries = tree->GetEntries();

	if(_debug) {
		printDebug << "check that the centers of mass bins of " << nrEntries << " entries in tree are at a known place, "
		           << "and map the " << nrMassBins << " mass bins to those entries." << std::endl;
	}

	for(Long64_t idx=0; idx<nrEntries; ++idx) {
		if(tree->GetEntry(idx) == 0) {
			printErr << "error while reading entry " << idx << " from tree." << std::endl;
			return false;
		}
		//FIXME: this would also be the place to select the best fit in case one file contains more than one fit result per mass bin
		const double mass = fit->massBinCenter();

		if(_debug) {
			printDebug << "entry " << idx << ": center of mass bin at " << mass << " GeV/c^2" << std::endl;
		}

		bool found = false;
		size_t idxMass=0;
		while(idxMass<nrMassBins) {
			if(std::abs(massBinCenters[idxMass]-mass) < 1000.*std::numeric_limits<double>::epsilon()) {
				found = true;
				break;
			}
			++idxMass;
		}

		if(not found) {
			printErr << "could not map mass bin centered at " << mass << " GeV/c^2 to a known mass bin." << std::endl;
			return false;
		}

		if(mapping[idxMass] != std::numeric_limits<Long64_t>::max()) {
			printErr << "cannot map tree entry " << idx << " to mass bin " << idxMass << " (" << massBinCenters[idxMass] << " GeV/c^2)  "
			         << "which is already mapped to tree entry " << mapping[idxMass] << "." << std::endl;
			return false;
		}

		if(_debug) {
			printDebug << "mapping mass bin " << idxMass << " (" << massBinCenters[idxMass] << " GeV/c^2) to tree entry " << idx << "." << std::endl;
		}
		mapping[idxMass] = idx;
	} // end loop over entries in tree

	// check that all mass bins are mapped
	for(size_t idx=0; idx<mapping.size(); ++idx) {
		if(mapping[idx] == std::numeric_limits<Long64_t>::max()) {
			printErr << "mass bin " << idx << " (" << massBinCenters[idx] << " GeV/c^2) not mapped." << std::endl;
			return false;
		}
	}

	if(_debug) {
		std::ostringstream output;
		for(size_t idx=0; idx<mapping.size(); ++idx) {
			output << " " << idx << "->" << mapping[idx];
		}
		printDebug << "established mapping:" << output.str() << std::endl;
	}

	return true;
}


bool
rpwa::resonanceFit::massDepFit::readFitResultMassBins(TTree* tree,
                                                      rpwa::fitResult* fit,
                                                      size_t& nrMassBins,
                                                      boost::multi_array<double, 1>& massBinCenters) const
{
	if(not tree or not fit) {
		printErr << "'tree' or 'fit' is not a pointer to a valid object." << std::endl;
		return false;
	}

	// extract data from tree
	const Long64_t nrEntries = tree->GetEntries();

	if(_debug) {
		printDebug << "getting center of mass bins from " << nrEntries << " entries in tree." << std::endl;
	}

	nrMassBins = 0;
	for(Long64_t idx=0; idx<nrEntries; ++idx) {
		if(tree->GetEntry(idx) == 0) {
			printErr << "error while reading entry " << idx << " from tree." << std::endl;
			return false;
		}
		const double newMass = fit->massBinCenter();

		if(_debug) {
			printDebug << "entry " << idx << ": center of mass bin at " << newMass << " GeV/c^2" << std::endl;
		}

		bool found = false;
		for(size_t idxMass = 0; idxMass < nrMassBins; ++idxMass) {
			if(std::abs(massBinCenters[idxMass]-newMass) < 1000.*std::numeric_limits<double>::epsilon()) {
				found = true;
				if(_debug) {
					printDebug << "this center of mass bin already was encountered before." << std::endl;
				}
				break;
			}
		}

		if(not found) {
			adjustSizeAndSet(massBinCenters, nrMassBins++, newMass);
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
			return false;
		}
	}

	return true;
}


bool
rpwa::resonanceFit::massDepFit::readFitResultMatrices(TTree* tree,
                                                      rpwa::fitResult* fit,
                                                      const std::vector<Long64_t>& mapping,
                                                      const double rescaleErrors,
                                                      std::vector<std::string>& waveNames,
                                                      boost::multi_array<std::complex<double>, 2>& productionAmplitudes,
                                                      boost::multi_array<TMatrixT<double>, 1>& productionAmplitudesCovariance,
                                                      boost::multi_array<std::complex<double>, 3>& spinDensityMatrices,
                                                      boost::multi_array<TMatrixT<double>, 1>& spinDensityCovarianceMatrices,
                                                      boost::multi_array<std::pair<double, double>, 2>& plottingIntensities,
                                                      boost::multi_array<std::pair<double, double>, 3>& plottingSpinDensityMatrixElementsReal,
                                                      boost::multi_array<std::pair<double, double>, 3>& plottingSpinDensityMatrixElementsImag,
                                                      boost::multi_array<std::pair<double, double>, 3>& plottingPhases) const
{
	if(not tree or not fit) {
		printErr << "'tree' or 'fit' is not a pointer to a valid object." << std::endl;
		return false;
	}

	if(_debug) {
		printDebug << "reading spin-density matrices for " << _nrWaves << " waves from fit result." << std::endl;
	}

	// read wave names from first fit result in tree
	waveNames.resize(_nrWaves);
	if(tree->GetEntry(0) == 0) {
		printErr << "error while reading entry " << 0 << " from tree." << std::endl;
		return false;
	}
	for(size_t idxWave = 0; idxWave < _nrWaves; ++idxWave) {
		int idx = fit->waveIndex(_waveNames[idxWave]);
		// try alternative wave names
		for(size_t idxAlt = 0; idxAlt < _waveNameAlternatives[idxWave].size(); ++idxAlt) {
			const int altIdx = fit->waveIndex(_waveNameAlternatives[idxWave][idxAlt]);
			if(altIdx != -1) {
				if(idx != -1) {
					printErr << "more than one wave name or alternative wave name is matching wave in fit result for wave '" << _waveNames[idxWave] << "'." << std::endl;
					return false;
				}
				idx = altIdx;
			}
		}
		if(idx == -1) {
			printErr << "wave '" << _waveNames[idxWave] << "' not in fit result." << std::endl;
			return false;
		}
		waveNames[idxWave] = fit->waveName(idx);
	}

	productionAmplitudes.resize(boost::extents[_maxMassBins][_nrWaves]);
	productionAmplitudesCovariance.resize(boost::extents[_maxMassBins]);

	spinDensityMatrices.resize(boost::extents[_maxMassBins][_nrWaves][_nrWaves]);
	spinDensityCovarianceMatrices.resize(boost::extents[_maxMassBins]);

	plottingIntensities.resize(boost::extents[_maxMassBins][_nrWaves]);
	plottingSpinDensityMatrixElementsReal.resize(boost::extents[_maxMassBins][_nrWaves][_nrWaves]);
	plottingSpinDensityMatrixElementsImag.resize(boost::extents[_maxMassBins][_nrWaves][_nrWaves]);
	plottingPhases.resize(boost::extents[_maxMassBins][_nrWaves][_nrWaves]);

	for(size_t idxMass = 0; idxMass < mapping.size(); ++idxMass) {
		if(_debug) {
			printDebug << "reading entry " << mapping[idxMass] << " for mass bin " << idxMass << " from tree." << std::endl;
		}
		// FIXME: in case of reading the fit result for a systematic tree this might happen, so this should be allowed in certain cases
		if(tree->GetEntry(mapping[idxMass]) == 0) {
			printErr << "error while reading entry " << mapping[idxMass] << " from tree." << std::endl;
			return false;
		}

		spinDensityCovarianceMatrices[idxMass].ResizeTo(_nrWaves * (_nrWaves+1), _nrWaves * (_nrWaves+1));
		for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
			const int idx = fit->waveIndex(waveNames[idxWave]);
			if(idx == -1) {
				printErr << "wave '" << _waveNames[idxWave] << "' not in fit result." << std::endl;
				return false;
			}

			plottingIntensities[idxMass][idxWave] = std::make_pair(fit->intensity(idx),
			                                                       fit->intensityErr(idx) * sqrt(rescaleErrors));

			for(size_t jdxWave=0; jdxWave<_nrWaves; ++jdxWave) {
				const int jdx = fit->waveIndex(waveNames[jdxWave]);
				if(jdx == -1) {
					printErr << "wave '" << _waveNames[jdxWave] << "' not in fit result." << std::endl;
					return false;
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

					const size_t idxCov = _nrWaves*(_nrWaves+1) - (_nrWaves-idxWave)*(_nrWaves-idxWave+1) + 2*(jdxWave-idxWave);
					spinDensityCovarianceMatrices[idxMass].SetSub(idxCov, idxCov, spinDensityMatrixElemCov);
				}
			}
		}

		// for the production amplitudes loop over the production
		// amplitudes of the fit result
		std::vector<unsigned int> prodAmpIndicesForCov(_nrWaves);
		for(unsigned int idxProdAmp=0; idxProdAmp < fit->nmbProdAmps(); ++idxProdAmp) {
			const std::string waveName = fit->waveNameForProdAmp(idxProdAmp);

			const std::map<std::string, size_t>::const_iterator it = _waveIndices.find(waveName);
			// most of the waves are ignored
			if(it == _waveIndices.end()) {
				continue;
			}
			size_t idxWave = it->second;

			int rank = fit->rankOfProdAmp(idxProdAmp);
			// TODO: multiple ranks, in that case also check that rank is not -1
			if(rank != 0) {
				printErr << "can only handle rank-1 fit (production amplitude '" << fit->prodAmpName(idxProdAmp)
				         << "' of wave '" << waveName << "' has rank " << rank << ")." << std::endl;
				return false;
			}

			productionAmplitudes[idxMass][idxWave] = fit->prodAmp(idxProdAmp);

			prodAmpIndicesForCov[idxWave] = idxProdAmp;
		}
		const TMatrixT<double> prodAmpCov = fit->prodAmpCov(prodAmpIndicesForCov) * rescaleErrors;
		productionAmplitudesCovariance[idxMass].ResizeTo(prodAmpCov);
		productionAmplitudesCovariance[idxMass] = prodAmpCov;

		if(_debug) {
			std::ostringstream outputProdAmp;
			std::ostringstream outputProdAmpCovariance;
			std::ostringstream output;
			std::ostringstream outputCovariance;

			outputProdAmp << " (";
			for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
				outputProdAmp << " " << productionAmplitudes[idxMass][idxWave];

				outputProdAmpCovariance << " (";
				output << " (";
				for(size_t jdxWave=0; jdxWave<_nrWaves; ++jdxWave) {
					output << " " << spinDensityMatrices[idxMass][idxWave][jdxWave];

					outputProdAmpCovariance << " (";
					for(size_t idx=0; idx<2; ++idx) {
						outputProdAmpCovariance << " (";
						for(size_t jdx=0; jdx<2; ++jdx) {
							outputProdAmpCovariance << " " << productionAmplitudesCovariance[idxMass](idxWave+idx, jdxWave+jdx);
						}
						outputProdAmpCovariance << " )";
					}
					outputProdAmpCovariance << " )";

					if(jdxWave >= idxWave) {
						const size_t idxCov = _nrWaves*(_nrWaves+1) - (_nrWaves-idxWave)*(_nrWaves-idxWave+1) + 2*(jdxWave-idxWave);
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

	return true;
}


bool
rpwa::resonanceFit::massDepFit::readFitResultIntegrals(TTree* tree,
                                                       rpwa::fitResult* fit,
                                                       const std::vector<Long64_t>& mapping,
                                                       const std::vector<std::string>& waveNames,
                                                       boost::multi_array<double, 2>& phaseSpaceIntegrals) const
{
	if(not tree or not fit) {
		printErr << "'tree' or 'fit' is not a pointer to a valid object." << std::endl;
		return false;
	}

	phaseSpaceIntegrals.resize(boost::extents[_maxMassBins][_nrWaves]);

	if(_debug) {
		printDebug << "reading phase-space integrals for " << _nrWaves << " waves from fit result." << std::endl;
	}

	for(size_t idxMass = 0; idxMass < mapping.size(); ++idxMass) {
		if(_debug) {
			printDebug << "reading entry " << mapping[idxMass] << " for mass bin " << idxMass << " from tree." << std::endl;
		}
		if(tree->GetEntry(mapping[idxMass]) == 0) {
			printErr << "error while reading entry " << mapping[idxMass] << " from tree." << std::endl;
			return false;
		}

		for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
			const double ps = fit->phaseSpaceIntegral(waveNames[idxWave]);
			phaseSpaceIntegrals[idxMass][idxWave] = ps;
		}
	}

	if(_debug) {
		for(size_t idxWave = 0; idxWave < _nrWaves; ++idxWave) {
			std::ostringstream output;
			for(size_t idxMass = 0; idxMass < mapping.size(); ++idxMass) {
				output << " " << phaseSpaceIntegrals[idxMass][idxWave];
			}
			printDebug << "phase-space integrals for wave '" << _waveNames[idxWave] << "' (" << idxWave << "):" << output.str() << std::endl;
		}
	}

	return true;
}


bool
rpwa::resonanceFit::massDepFit::prepareMassLimits(const std::vector<size_t>& nrMassBins,
                                                  const boost::multi_array<double, 2>& massBinCenters)
{
	_waveMassBinLimits.resize(boost::extents[_nrBins][_nrWaves]);
	_wavePairMassBinLimits.resize(boost::extents[_nrBins][_nrWaves][_nrWaves]);

	for(size_t idxBin = 0; idxBin < _nrBins; ++idxBin) {
		if(not prepareMassLimit(nrMassBins[idxBin], massBinCenters[idxBin], idxBin)) {
			printErr << "error while determine bins to use in fit for bin " << idxBin << "." << std::endl;
			return false;
		}
	}

	return true;
}


bool
rpwa::resonanceFit::massDepFit::prepareMassLimit(const size_t nrMassBins,
                                                 const boost::multi_array<double, 1>& massBinCenters,
                                                 const size_t idxBin)
{
	if(_debug) {
		printDebug << "determine which mass bins to use in the fit for " << nrMassBins << " mass bins, center of first and last mass bins: "
		           << massBinCenters[0] << " and " << massBinCenters[nrMassBins - 1] << " GeV/c^2." << std::endl;
	}

	for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
		size_t binFirst = 0;
		size_t binLast = nrMassBins-1;
		for(size_t idxMass = 0; idxMass < nrMassBins; ++idxMass) {
			if(massBinCenters[idxMass] < _waveMassLimits[idxWave].first) {
				binFirst = idxMass+1;
			}
			if(massBinCenters[idxMass] == _waveMassLimits[idxWave].first) {
				binFirst = idxMass;
			}
			if(massBinCenters[idxMass] <= _waveMassLimits[idxWave].second) {
				binLast = idxMass;
			}
		}
		if(_waveMassLimits[idxWave].first < 0) {
			binFirst = 0;
		}
		if(_waveMassLimits[idxWave].second < 0) {
			binLast = nrMassBins-1;
		}

		const double massStep = (massBinCenters[nrMassBins-1] - massBinCenters[0]) / (nrMassBins - 1);
		const double massFirst = massBinCenters[binFirst] - massStep/2.;
		const double massLast = massBinCenters[binLast] + massStep/2.;
		if(_debug) {
			printDebug << idxWave << ": " << _waveNames[idxWave] << ": "
			           << "mass range: " << massFirst << "-" << massLast << " GeV/c^2, "
			           << "bin range: " << binFirst << "-" << binLast << std::endl;
		}
		_waveMassBinLimits[idxBin][idxWave] = std::make_pair(binFirst, binLast);
	}

	for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
		for(size_t jdxWave=0; jdxWave<_nrWaves; ++jdxWave) {
			_wavePairMassBinLimits[idxBin][idxWave][jdxWave] = std::make_pair(std::max(_waveMassBinLimits[idxBin][idxWave].first,  _waveMassBinLimits[idxBin][jdxWave].first),
			                                                                  std::min(_waveMassBinLimits[idxBin][idxWave].second, _waveMassBinLimits[idxBin][jdxWave].second));
		}
	}

	if(_debug) {
		printDebug << "waves and mass limits:" << std::endl;
		for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
			std::ostringstream output;
			for(size_t jdxWave=0; jdxWave<_nrWaves; ++jdxWave) {
				output << _wavePairMassBinLimits[idxBin][idxWave][jdxWave].first << "-" << _wavePairMassBinLimits[idxBin][idxWave][jdxWave].second << " ";
			}
			printDebug << _waveNames[idxWave] << " " << _waveMassBinLimits[idxBin][idxWave].first << "-" << _waveMassBinLimits[idxBin][idxWave].second
			           << ": " << output.str() << std::endl;
		}
	}

	return true;
}


bool
rpwa::resonanceFit::massDepFit::createPlots(const rpwa::resonanceFit::dataConstPtr& fitData,
                                            const rpwa::resonanceFit::modelConstPtr& fitModel,
                                            const rpwa::resonanceFit::parameters& fitParameters,
                                            rpwa::resonanceFit::cache& cache,
                                            TFile* outFile,
                                            const bool rangePlotting,
                                            const size_t extraBinning) const
{
	if(_debug) {
		printDebug << "start creating plots." << std::endl;
	}

	for(size_t idxBin=0; idxBin<_nrBins; ++idxBin) {
		TDirectory* outDirectory = NULL;
		if(_nrBins == 1) {
			outDirectory = outFile;
		} else {
			std::ostringstream name;
			name << "bin" << idxBin;
			outDirectory = outFile->mkdir(name.str().c_str());
		}

		for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
			if(not createPlotsWave(fitData, fitModel, fitParameters, cache, outDirectory, rangePlotting, extraBinning, idxWave, idxBin)) {
				printErr << "error while creating intensity plots for wave '" << _waveNames[idxWave] << "' in bin " << idxBin << "." << std::endl;
				return false;
			}
		}

		for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
			for(size_t jdxWave=idxWave+1; jdxWave<_nrWaves; ++jdxWave) {
				if(not createPlotsWavePair(fitData, fitModel, fitParameters, cache, outDirectory, rangePlotting, extraBinning, idxWave, jdxWave, idxBin)) {
					printErr << "error while creating intensity plots for wave pair '" << _waveNames[idxWave] << "' and '" << _waveNames[jdxWave] << "' in bin " << idxBin << "." << std::endl;
					return false;
				}
			}
		}

		if(fitModel->getFsmd() and (fitModel->getFsmd()->getNrBins() != 1 or not _sameMassBinning)) {
			if(not createPlotsFsmd(fitData, fitModel, fitParameters, cache, outDirectory, rangePlotting, extraBinning, idxBin)) {
				printErr << "error while creating plots for final-state mass-dependence in bin " << idxBin << "." << std::endl;
				return false;
			}
		}
	}

	if(_nrBins != 1 and _sameMassBinning and fitModel->isMappingEqualInAllBins()) {
		for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
			if(not createPlotsWaveSum(fitData, fitModel, fitParameters, cache, outFile, rangePlotting, extraBinning, idxWave)) {
				printErr << "error while creating intensity plots for wave '" << _waveNames[idxWave] << "' for sum over all bins." << std::endl;
				return false;
			}
		}
	}

	if(fitModel->getFsmd() and (fitModel->getFsmd()->getNrBins() == 1 and _sameMassBinning)) {
		if(not createPlotsFsmd(fitData, fitModel, fitParameters, cache, outFile, rangePlotting, extraBinning, 0)) {
			printErr << "error while creating plots for final-state mass-dependence." << std::endl;
			return false;
		}
	}

	if(_debug) {
		printDebug << "finished creating plots." << std::endl;
	}

	return true;
}


bool
rpwa::resonanceFit::massDepFit::createPlotsWave(const rpwa::resonanceFit::dataConstPtr& fitData,
                                                const rpwa::resonanceFit::modelConstPtr& fitModel,
                                                const rpwa::resonanceFit::parameters& fitParameters,
                                                rpwa::resonanceFit::cache& cache,
                                                TDirectory* outDirectory,
                                                const bool rangePlotting,
                                                const size_t extraBinning,
                                                const size_t idxWave,
                                                const size_t idxBin) const
{
	if(_debug) {
		printDebug << "start creating plots for wave '" << _waveNames[idxWave] << "' in bin " << idxBin << "." << std::endl;
	}

	TMultiGraph graphs;
	graphs.SetName(_waveNames[idxWave].c_str());
	graphs.SetTitle(_waveNames[idxWave].c_str());

	TGraphErrors* systematics = NULL;
	if(_nrSystematics[idxBin] > 0) {
		systematics = new TGraphErrors;
		systematics->SetName((_waveNames[idxWave] + "__sys").c_str());
		systematics->SetTitle((_waveNames[idxWave] + "__sys").c_str());
		systematics->SetLineColor(kAzure-9);
		systematics->SetFillColor(kAzure-9);
		graphs.Add(systematics, "2");
	}

	TGraphErrors* data = new TGraphErrors;
	data->SetName((_waveNames[idxWave] + "__data").c_str());
	data->SetTitle((_waveNames[idxWave] + "__data").c_str());
	graphs.Add(data, "P");

	TGraph* fit = new TGraph;
	fit->SetName((_waveNames[idxWave] + "__fit").c_str());
	fit->SetTitle((_waveNames[idxWave] + "__fit").c_str());
	fit->SetLineColor(kRed);
	fit->SetLineWidth(2);
	fit->SetMarkerColor(kRed);
	graphs.Add(fit, "L");

	TGraph* phaseSpace = new TGraph;
	phaseSpace->SetName((_waveNames[idxWave] + "__ps").c_str());
	phaseSpace->SetTitle((_waveNames[idxWave] + "__ps").c_str());
	graphs.Add(phaseSpace, "L");

	const std::vector<std::pair<size_t, size_t> >& compChannel = fitModel->getComponentChannel(idxBin, idxWave);
	std::vector<TGraph*> components;
	for(size_t idxComponents=0; idxComponents<compChannel.size(); ++idxComponents) {
		const size_t idxComponent = compChannel[idxComponents].first;
		TGraph* component = new TGraph;
		component->SetName((_waveNames[idxWave] + "__" + fitModel->getComponent(idxComponent)->getName()).c_str());
		component->SetTitle((_waveNames[idxWave] + "__" + fitModel->getComponent(idxComponent)->getName()).c_str());

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

		if(_nrSystematics[idxBin] > 0) {
			const double minSI = fitData->sysPlottingIntensities()[idxBin][idxMass][idxWave].first;
			const double maxSI = fitData->sysPlottingIntensities()[idxBin][idxMass][idxWave].second;
			systematics->SetPoint(point, mass, (maxSI+minSI)/2.);
			systematics->SetPointError(point, halfBin, (maxSI-minSI)/2.);
			maxIE = std::max(maxIE, maxSI);
		}
	}

	// plot fit, either over full or limited mass range
	const size_t firstPoint = rangePlotting ? (extraBinning*_waveMassBinLimits[idxBin][idxWave].first) : 0;
	const size_t lastPoint = rangePlotting ? (extraBinning*_waveMassBinLimits[idxBin][idxWave].second) : (extraBinning*(fitData->nrMassBins()[idxBin]-1));
	for(size_t point=firstPoint; point<=lastPoint; ++point) {
		const size_t idxMass = (point%extraBinning == 0) ? (point/extraBinning) : std::numeric_limits<size_t>::max();
		const double massStep = (fitData->massBinCenters()[idxBin][fitData->nrMassBins()[idxBin]-1] - fitData->massBinCenters()[idxBin][0]) / (fitData->nrMassBins()[idxBin] - 1) / extraBinning;
		const double mass = (idxMass != std::numeric_limits<size_t>::max()) ? fitData->massBinCenters()[idxBin][idxMass] : (fitData->massBinCenters()[idxBin][point/extraBinning] + (point%extraBinning) * massStep);

		const double intensity = fitModel->intensity(fitParameters, cache, idxWave, idxBin, mass, idxMass);
		fit->SetPoint(point-firstPoint, mass, intensity);
		maxIE = std::max(maxIE, intensity);

		for(size_t idxComponents=0; idxComponents<compChannel.size(); ++idxComponents) {
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
	boost::multi_array<double, 3>::const_array_view<1>::type viewInt = _inPhaseSpaceIntegrals[boost::indices[idxBin][boost::multi_array<double, 3>::index_range(0, fitData->nrMassBins()[idxBin])][idxWave]];
	ROOT::Math::Interpolator phaseSpaceInterpolator(std::vector<double>(viewM.begin(), viewM.end()), std::vector<double>(viewInt.begin(), viewInt.end()), ROOT::Math::Interpolation::kLINEAR);

	// plot phase-space
	double maxP = -std::numeric_limits<double>::max();
	for(size_t point = 0; point <= (extraBinning*(fitData->nrMassBins()[idxBin]-1)); ++point) {
		const size_t idxMass = (point%extraBinning == 0) ? (point/extraBinning) : std::numeric_limits<size_t>::max();
		const double massStep = (fitData->massBinCenters()[idxBin][fitData->nrMassBins()[idxBin]-1] - fitData->massBinCenters()[idxBin][0]) / (fitData->nrMassBins()[idxBin] - 1) / extraBinning;
		const double mass = (idxMass != std::numeric_limits<size_t>::max()) ? fitData->massBinCenters()[idxBin][idxMass] : (fitData->massBinCenters()[idxBin][point/extraBinning] + (point%extraBinning) * massStep);

		double ps = pow((idxMass != std::numeric_limits<size_t>::max()) ? _inPhaseSpaceIntegrals[idxBin][idxMass][idxWave] : phaseSpaceInterpolator.Eval(mass), 2);
		if(fitModel->getFsmd()) {
			ps *= std::norm(fitModel->getFsmd()->val(fitParameters, cache, idxBin, mass, idxMass));
		}
		phaseSpace->SetPoint(point, mass, ps);
		maxP = std::max(maxP, ps);
	}

	// scale phase-space graph to half-height of intensity graphs
	for(Int_t idx=0; idx<phaseSpace->GetN(); ++idx) {
		double x, y;
		phaseSpace->GetPoint(idx, x, y);
		phaseSpace->SetPoint(idx, x, y * 0.5 * maxIE/maxP);
	}

	outDirectory->cd();
	graphs.Write();

	return true;
}


bool
rpwa::resonanceFit::massDepFit::createPlotsWaveSum(const rpwa::resonanceFit::dataConstPtr& fitData,
                                                   const rpwa::resonanceFit::modelConstPtr& fitModel,
                                                   const rpwa::resonanceFit::parameters& fitParameters,
                                                   rpwa::resonanceFit::cache& cache,
                                                   TDirectory* outDirectory,
                                                   const bool rangePlotting,
                                                   const size_t extraBinning,
                                                   const size_t idxWave) const
{
	if(_debug) {
		printDebug << "start creating plots for wave '" << _waveNames[idxWave] << "' for sum over all bins." << std::endl;
	}

	// all mass binnings must be the same to be able to create the sum plots
	if(not _sameMassBinning or not fitModel->isMappingEqualInAllBins()) {
		printErr << "cannot create plots for wave '" << _waveNames[idxWave] << "' for sum over all bins if the bins used different mass binnings." << std::endl;
		return false;
	}
	const size_t idxBin = 0;

	TMultiGraph graphs;
	graphs.SetName(_waveNames[idxWave].c_str());
	graphs.SetTitle(_waveNames[idxWave].c_str());

	TGraphErrors* data = new TGraphErrors;
	data->SetName((_waveNames[idxWave] + "__data").c_str());
	data->SetTitle((_waveNames[idxWave] + "__data").c_str());
	graphs.Add(data, "P");

	TGraph* fit = new TGraph;
	fit->SetName((_waveNames[idxWave] + "__fit").c_str());
	fit->SetTitle((_waveNames[idxWave] + "__fit").c_str());
	fit->SetLineColor(kRed);
	fit->SetLineWidth(2);
	fit->SetMarkerColor(kRed);
	graphs.Add(fit, "L");

	const std::vector<std::pair<size_t, size_t> >& compChannel = fitModel->getComponentChannel(idxBin, idxWave);
	std::vector<TGraph*> components;
	for(size_t idxComponents=0; idxComponents<compChannel.size(); ++idxComponents) {
		const size_t idxComponent = compChannel[idxComponents].first;
		TGraph* component = new TGraph;
		component->SetName((_waveNames[idxWave] + "__" + fitModel->getComponent(idxComponent)->getName()).c_str());
		component->SetTitle((_waveNames[idxWave] + "__" + fitModel->getComponent(idxComponent)->getName()).c_str());

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
		for(size_t idxBin=0; idxBin<_nrBins; ++idxBin) {
			sum += fitData->plottingIntensities()[idxBin][idxMass][idxWave].first;
			error2 += std::pow(fitData->plottingIntensities()[idxBin][idxMass][idxWave].second, 2);
		}
		data->SetPoint(point, mass, sum);
		data->SetPointError(point, halfBin, sqrt(error2));
	}

	// plot fit, either over full or limited mass range
	const size_t firstPoint = rangePlotting ? (extraBinning*_waveMassBinLimits[idxBin][idxWave].first) : 0;
	const size_t lastPoint = rangePlotting ? (extraBinning*_waveMassBinLimits[idxBin][idxWave].second) : (extraBinning*(fitData->nrMassBins()[idxBin]-1));
	for(size_t point=firstPoint; point<=lastPoint; ++point) {
		const size_t idxMass = (point%extraBinning == 0) ? (point/extraBinning) : std::numeric_limits<size_t>::max();
		const double massStep = (fitData->massBinCenters()[idxBin][fitData->nrMassBins()[idxBin]-1] - fitData->massBinCenters()[idxBin][0]) / (fitData->nrMassBins()[idxBin] - 1) / extraBinning;
		const double mass = (idxMass != std::numeric_limits<size_t>::max()) ? fitData->massBinCenters()[idxBin][idxMass] : (fitData->massBinCenters()[idxBin][point/extraBinning] + (point%extraBinning) * massStep);

		double sum = 0.;
		for(size_t idxBin=0; idxBin<_nrBins; ++idxBin) {
			sum += fitModel->intensity(fitParameters, cache, idxWave, idxBin, mass, idxMass);
		}
		fit->SetPoint(point-firstPoint, mass, sum);

		for(size_t idxComponents=0; idxComponents<compChannel.size(); ++idxComponents) {
			const size_t idxComponent = compChannel[idxComponents].first;
			const size_t idxChannel = compChannel[idxComponents].second;

			double sum = 0.;
			for(size_t idxBin=0; idxBin<_nrBins; ++idxBin) {
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

	return true;
}


bool
rpwa::resonanceFit::massDepFit::createPlotsWavePair(const rpwa::resonanceFit::dataConstPtr& fitData,
                                                    const rpwa::resonanceFit::modelConstPtr& fitModel,
                                                    const rpwa::resonanceFit::parameters& fitParameters,
                                                    rpwa::resonanceFit::cache& cache,
                                                    TDirectory* outDirectory,
                                                    const bool rangePlotting,
                                                    const size_t extraBinning,
                                                    const size_t idxWave,
                                                    const size_t jdxWave,
                                                    const size_t idxBin) const
{
	if(_debug) {
		printDebug << "start creating plots for wave pair '" << _waveNames[idxWave] << "' and '" << _waveNames[jdxWave] << "' in bin " << idxBin << "." << std::endl;
	}

	const std::string realName = _waveNames[idxWave] + "__" + _waveNames[jdxWave] + "__real";
	const std::string imagName = _waveNames[idxWave] + "__" + _waveNames[jdxWave] + "__imag";
	const std::string phaseName = _waveNames[idxWave] + "__" + _waveNames[jdxWave] + "__phase";

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
	if(_nrSystematics[idxBin] > 0) {
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

		if(_nrSystematics[idxBin] > 0) {
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
	const size_t firstPoint = rangePlotting ? (extraBinning*_wavePairMassBinLimits[idxBin][idxWave][jdxWave].first) : 0;
	const size_t lastPoint = rangePlotting ? (extraBinning*_wavePairMassBinLimits[idxBin][idxWave][jdxWave].second) : (extraBinning*(fitData->nrMassBins()[idxBin]-1));
	for(size_t point=firstPoint; point<=lastPoint; ++point) {
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
			for(int offs=-5; offs<6; ++offs) {
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

		if (idxMass != std::numeric_limits<size_t>::max()) {
			int bestOffs = 0;
			double bestDiff = std::numeric_limits<double>::max();

			double data;
			phaseData->GetPoint(idxMass, x, data);
			for(int offs=-5; offs<6; ++offs) {
				if(std::abs(data + offs*360. - valueFit) < bestDiff) {
					bestDiff = std::abs(data + offs*360. - valueFit);
					bestOffs = offs;
				}
			}

			phaseData->SetPoint(idxMass, x, data + bestOffs*360.);
			if(_nrSystematics[idxBin] > 0) {
				phaseSystematics->GetPoint(idxMass, x, data);
				phaseSystematics->SetPoint(idxMass, x, data + bestOffs*360.);
			}
		}

		// check that this mass bin should be taken into account for this
		// combination of waves
		if(point < firstPoint || point > lastPoint) {
			continue;
		}

		phaseFit->SetPoint(point-firstPoint, mass, valueFit);
	}

	outDirectory->cd();
	real.Write();
	imag.Write();
	phase.Write();

	return true;
}


bool
rpwa::resonanceFit::massDepFit::createPlotsFsmd(const rpwa::resonanceFit::dataConstPtr& fitData,
                                                const rpwa::resonanceFit::modelConstPtr& fitModel,
                                                const rpwa::resonanceFit::parameters& fitParameters,
                                                rpwa::resonanceFit::cache& cache,
                                                TDirectory* outDirectory,
                                                const bool /*rangePlotting*/,
                                                const size_t extraBinning,
                                                const size_t idxBin) const
{
	if(_debug) {
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

	return true;
}
