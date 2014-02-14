
#include "massDepFitLikeli.h"
#include "pwacomponent.h"


using namespace rpwa;
using namespace std;


massDepFitLikeli*
massDepFitLikeli::Clone() const {
	return new massDepFitLikeli(*this);
}


unsigned int
massDepFitLikeli::NDim() const {
	return _compset->numPar();
}


unsigned int
massDepFitLikeli::NDataPoints() const {
	// calculate data points:
	// * diagonal elements are real numbers
	// * non-diagonal elements are complex numbers
	// * remember (Re,Im) => factor 2
	// * diagonal elements are only checked once, of diagonal elements with
	//   the two different combinations (i,j) and (j,i)
	unsigned int nrPts(0);

	for(size_t idxWave=0; idxWave<_wavePairMassBinLimits.size(); ++idxWave) {
		for(size_t jdxWave=0; jdxWave<_wavePairMassBinLimits[idxWave].size(); ++jdxWave) {
			nrPts += _wavePairMassBinLimits[idxWave][jdxWave].second - _wavePairMassBinLimits[idxWave][jdxWave].first + 1;
		}
	}
  
	return nrPts;
}


void
massDepFitLikeli::init(pwacompset* compset,
                       const std::vector<double>& massBinCenters,
                       const std::vector<spinDensityMatrixType>& spinDensityMatrices,
                       const std::vector<spinDensityCovarianceMatrixType>& spinDensityCovarianceMatrices,
                       const std::vector<std::vector<std::pair<size_t, size_t> > >& wavePairMassBinLimits,
                       bool useCovariance)
{
	_compset = compset;

	_massBinCenters = massBinCenters;

	_spinDensityMatrices = spinDensityMatrices;
	_spinDensityCovarianceMatrices = spinDensityCovarianceMatrices;

	_wavePairMassBinLimits = wavePairMassBinLimits;

	_useCovariance = useCovariance;
}


double
massDepFitLikeli::DoEval(const double* par) const {
	const size_t nrWaves = _wavePairMassBinLimits.size();
	const size_t nrMassBins = _massBinCenters.size();
     
	// set parameters for resonances, background and phase space
	_compset->setPar(par);

	double chi2=0;
 
	// loop over mass-bins
	for(unsigned idxMass=0; idxMass<nrMassBins; ++idxMass) {
		const double mass = _massBinCenters[idxMass];
		const spinDensityMatrixType& spinDensityMatrix = _spinDensityMatrices[idxMass];
		const spinDensityCovarianceMatrixType& spinDensityCovarianceMatrix = _spinDensityCovarianceMatrices[idxMass];

		// sum over the contributions to chi2 -> rho_ij
		for(size_t idxWave=0; idxWave<nrWaves; ++idxWave) {
			for(size_t jdxWave=idxWave; jdxWave<nrWaves; ++jdxWave) {
				// check that this mass bin should be taken into account for this
				// combination of waves
				if(idxMass < _wavePairMassBinLimits[idxWave][jdxWave].first || idxMass > _wavePairMassBinLimits[idxWave][jdxWave].second) {
					continue;
				}

				// calculate target spin density matrix element
				const complex<double> rhoFit = _compset->overlap(idxWave, jdxWave, mass, idxMass);

				const complex<double> rhoDiff = rhoFit - spinDensityMatrix(idxWave, jdxWave);

				const complexCovarianceMatrixType& covariance = spinDensityCovarianceMatrix(idxWave, jdxWave);

				double dchi;
				if(idxWave==jdxWave) {
					dchi = norm(rhoDiff) / covariance(0,0);
				} else {
					if(_useCovariance) {
						dchi  = rhoDiff.real()*rhoDiff.real() * covariance(1,1);
						dchi -= rhoDiff.real()*rhoDiff.imag() * (covariance(0,1) + covariance(1,0));
						dchi += rhoDiff.imag()*rhoDiff.imag() * covariance(0,0);

						dchi /= covariance(0,0)*covariance(1,1) - covariance(0,1)*covariance(1,0);
					} else {
						dchi = rhoDiff.real()*rhoDiff.real()/covariance(0,0) + rhoDiff.imag()*rhoDiff.imag()/covariance(1,1);
					}
				}
				chi2 += dchi;
			} // end loop over jdxWave
		} // end loop over idxWave
	} // end loop over mass-bins
  
	return chi2;
} // end DoEval
