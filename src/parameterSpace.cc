#include "parameterSpace.h"

#include "ampIntegralMatrix.h"
#include "reportingUtils.hpp"

#include <TMath.h>

using namespace std;
using namespace rpwa;


void parameterSpace::pickAngles() {

	printErr << "This method is not finished yet." << endl;

	for(unsigned int i = 0; i < _nDimensions; ++i) {
		_phi[i] = (3. / _nDimensions) * i;
		_sigma[i] = calculateSigma(i);
	}

	initialize();

}


void parameterSpace::setAngles(const std::vector<double>& phi) {

	if(phi.size() != _phi.size()) {
		printErr << "Size mismatch between number of phis provided ("
		         << phi.size() << ") and number of phis needed ("
		         << _phi.size() << "). Aborting..." << endl;
		throw;
	}
	for(unsigned int i = 0; i < _phi.size(); ++i) {
		_phi[i] = phi[i];
		_sigma[i] = calculateSigma(i);
	}

	initialize();

}


double parameterSpace::getDS(const unsigned int& phiIndex) const
{
	double dS = 0.;
	for(unsigned int i = 0; i < (2*_nDimensions); ++i) {
		double x = (getDRDPhi(phiIndex) * _sigma[i]) + (getR() * getDSigmaDPhi(i, phiIndex));
		dS += x * x;
	}
	return sqrt(dS);
}


double parameterSpace::getDA() const
{
	double dA = 1.;
	for(unsigned int i = 0; i < ((2*_nDimensions)-1); ++i) {
		dA *= getDS(i);
	}
	return dA;
}


void parameterSpace::calculateAllDRDPhi() {
	for(unsigned int i = 0; i < ((2*_nDimensions)-1); ++i) {
		_DRDPhi[i] = calculateDRDPhi(i);
	}
}


double parameterSpace::calculateDRDPhi(const unsigned int& phiIndex) const
{
	const double& rho = getRho();
	return -(((double)_nEvents) / (2. * sqrt(rho) * rho * rho)) * calculateDRhoDPhi(phiIndex);
}


void parameterSpace::calculateRho()
{

	complex<double> retval(0., 0.);

	for(unsigned int alpha = 0; alpha < (2*_nDimensions); ++alpha) {
		for(unsigned int beta = 0; beta < (2*_nDimensions); ++beta) {
			const complex<double>& IAlphaBeta = _integralMatrix.matrix()[alpha][beta];
			retval = _sigma[alpha] * _sigma[beta] * IAlphaBeta;
		}
	}

	for(unsigned int alpha = _nDimensions; alpha < (2*_nDimensions); ++alpha) {
		for(unsigned int beta = 0; beta < _nDimensions; ++beta) {
			double lAlphaBeta = _integralMatrix.matrix()[alpha][beta].imag();
			retval -= 2. * _sigma[alpha] * _sigma[beta] * lAlphaBeta;
		}
	}

	if(fabs(retval.imag()) > 1e-15) {
		printErr << "Calculated rho with finite imaginary part. Aborting..." << endl;
		throw;
	}

	_rho = retval.real();

}


double parameterSpace::calculateDRhoDPhi(const unsigned int& phiIndex) const
{

	double retval = 0.;

	for(unsigned int beta = 0; beta < (2*_nDimensions); ++beta) {
		for(unsigned int alpha = phiIndex + 1; alpha < (2*_nDimensions); ++alpha) {
			const double& kAlphaBeta = _integralMatrix.matrix()[alpha][beta].real();
			retval += (2 * kAlphaBeta * _sigma[alpha] * _sigma[beta]) / TMath::Tan(_phi[phiIndex]);
		}
	}

	for(unsigned int alpha = 0; alpha < (2*_nDimensions); ++alpha) {
		const double& kAlphaBeta = _integralMatrix.matrix()[alpha][phiIndex].real();
		retval -= 2 * kAlphaBeta * _sigma[alpha] * _sigma[phiIndex] * TMath::Tan(_phi[phiIndex]);
	}

	if(phiIndex < _nDimensions) {

		for(unsigned int alpha = _nDimensions; alpha < (2*_nDimensions); ++alpha) {
			for(unsigned int beta = 0; beta < _nDimensions; ++beta) {
				double lAlphaBeta = _integralMatrix.matrix()[alpha][beta].imag();
				double factor = 1.;
				if(beta > phiIndex) {
					factor = 2.;
				}
				retval -= (factor * 2. * _sigma[alpha] * _sigma[beta] * lAlphaBeta) / TMath::Tan(_phi[phiIndex]);
			}
		}

		for(unsigned int alpha = _nDimensions; alpha < (2*_nDimensions); ++alpha) {
			const double& lAlphaBeta = _integralMatrix.matrix()[alpha][phiIndex].imag();
			retval += 2 * lAlphaBeta * _sigma[alpha] * _sigma[phiIndex] * TMath::Tan(_phi[phiIndex]);
		}

	} else if(phiIndex == _nDimensions) {

		for(unsigned int alpha = _nDimensions; alpha < (2*_nDimensions); ++alpha) {
			for(unsigned int beta = 0; beta < _nDimensions; ++beta) {
				double lAlphaBeta = _integralMatrix.matrix()[alpha][beta].imag();
				retval -= (2. * _sigma[alpha] * _sigma[beta] * lAlphaBeta) / TMath::Tan(_phi[phiIndex]);
			}
		}

		for(unsigned int alpha = _nDimensions; alpha < (2*_nDimensions); ++alpha) {
			const double& lAlphaBeta = _integralMatrix.matrix()[alpha][phiIndex].imag();
			retval += 2 * lAlphaBeta * _sigma[alpha] * _sigma[phiIndex] * TMath::Tan(_phi[phiIndex]);
		}

	} else {

		for(unsigned int beta = 0; beta < _nDimensions; ++beta) {
			for(unsigned int alpha = phiIndex+1; alpha < (2*_nDimensions); ++alpha) {
				double lAlphaBeta = _integralMatrix.matrix()[alpha][beta].imag();
				retval -= (2. * _sigma[alpha] * _sigma[beta] * lAlphaBeta) / TMath::Tan(_phi[phiIndex]);
			}
		}

		for(unsigned int beta = 0; beta < _nDimensions; ++beta) {
			const double& lAlphaBeta = _integralMatrix.matrix()[phiIndex][beta].imag();
			retval += 2 * lAlphaBeta * _sigma[phiIndex] * _sigma[beta] * TMath::Tan(_phi[phiIndex]);
		}

	}

	return retval;

}


double parameterSpace::calculateSigma(const unsigned int& index) const
{
	double sigma = TMath::Cos(_phi[index]);
	for(unsigned int i = 0; i < index; ++i) {
		sigma *= TMath::Sin(_phi[i]);
	}
	return sigma;
}


double parameterSpace::getDSigmaDPhi(const unsigned int& sigmaIndex, const unsigned int& phiIndex) const
{
	if(phiIndex > sigmaIndex) {
		return 0;
	} else if(phiIndex == sigmaIndex) {
		return (-TMath::Cos(_phi[phiIndex]) * _sigma[sigmaIndex]);
	} else {
		return (_sigma[sigmaIndex] / TMath::Cos(_phi[phiIndex]));
	}
}
