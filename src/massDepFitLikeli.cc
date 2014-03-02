#include "massDepFitLikeli.h"

#include "massDepFitModel.h"

using namespace std;


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
rpwa::massDepFit::likelihood::init(rpwa::massDepFit::model* compset,
                                   const std::vector<double>& massBinCenters,
                                   const boost::multi_array<std::complex<double>, 3>& spinDensityMatrices,
                                   const boost::multi_array<double, 5>& spinDensityCovarianceMatrices,
                                   const boost::multi_array<std::pair<size_t, size_t>, 2>& wavePairMassBinLimits,
                                   bool useCovariance)
{
	_compset = compset;

	_massBinCenters = massBinCenters;

	_spinDensityMatrices.resize(std::vector<size_t>(spinDensityMatrices.shape(), spinDensityMatrices.shape()+spinDensityMatrices.num_dimensions()));
	_spinDensityMatrices = spinDensityMatrices;
	_spinDensityCovarianceMatrices.resize(std::vector<size_t>(spinDensityCovarianceMatrices.shape(), spinDensityCovarianceMatrices.shape()+spinDensityCovarianceMatrices.num_dimensions()));
	_spinDensityCovarianceMatrices = spinDensityCovarianceMatrices;

	_wavePairMassBinLimits.resize(std::vector<size_t>(wavePairMassBinLimits.shape(), wavePairMassBinLimits.shape()+wavePairMassBinLimits.num_dimensions()));
	_wavePairMassBinLimits = wavePairMassBinLimits;

	_useCovariance = useCovariance;

	_nrMassBins = _massBinCenters.size();
	_nrWaves = _wavePairMassBinLimits.size();
}


double
rpwa::massDepFit::likelihood::DoEval(const double* par) const {
	// set parameters for resonances, background and phase space
	_compset->setParameters(par);

	double chi2=0;

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
				// FIXME: replace 0 by idxBin
				const complex<double> rhoFit = _compset->spinDensityMatrix(idxWave, jdxWave, 0, mass, idxMass);

				const complex<double> rhoDiff = rhoFit - _spinDensityMatrices[idxMass][idxWave][jdxWave];

				double dchi;
				if(idxWave==jdxWave) {
					dchi = norm(rhoDiff) / _spinDensityCovarianceMatrices[idxMass][idxWave][jdxWave][0][0];
				} else {
					if(_useCovariance) {
						dchi  = rhoDiff.real()*rhoDiff.real() * _spinDensityCovarianceMatrices[idxMass][idxWave][jdxWave][1][1];
						dchi -= rhoDiff.real()*rhoDiff.imag() * _spinDensityCovarianceMatrices[idxMass][idxWave][jdxWave][0][1];
						dchi -= rhoDiff.real()*rhoDiff.imag() * _spinDensityCovarianceMatrices[idxMass][idxWave][jdxWave][1][0];
						dchi += rhoDiff.imag()*rhoDiff.imag() * _spinDensityCovarianceMatrices[idxMass][idxWave][jdxWave][0][0];

						dchi /= _spinDensityCovarianceMatrices[idxMass][idxWave][jdxWave][0][0]*_spinDensityCovarianceMatrices[idxMass][idxWave][jdxWave][1][1]
						        - _spinDensityCovarianceMatrices[idxMass][idxWave][jdxWave][0][1]*_spinDensityCovarianceMatrices[idxMass][idxWave][jdxWave][1][0];
					} else {
						dchi  = rhoDiff.real()*rhoDiff.real() / _spinDensityCovarianceMatrices[idxMass][idxWave][jdxWave][0][0];
						dchi += rhoDiff.imag()*rhoDiff.imag() / _spinDensityCovarianceMatrices[idxMass][idxWave][jdxWave][1][1];
					}
				}
				chi2 += dchi;
			} // end loop over jdxWave
		} // end loop over idxWave
	} // end loop over mass-bins

	return chi2;
}
