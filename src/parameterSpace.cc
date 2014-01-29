#include "parameterSpace.h"

#include <TMath.h>
#include <TRandom3.h>

#include "ampIntegralMatrix.h"
#include "reportingUtils.hpp"


using namespace boost;
using namespace std;
using namespace rpwa;


std::vector<double> parameterSpace::convertToSphericalCoordinates(std::vector<double> x)
{

	std::vector<double> sphericalCoords(x.size(), 0.);

	double r = 0.;
	vector<double> xSums(x.size(), 0.);
	for(int i = x.size()-1; i >= 0; --i) {
		unsigned int j = (unsigned int)i;
		r += x[j] * x[j];
		xSums[j] = sqrt(r);
	}
	r = sqrt(r);
	sphericalCoords[0] = r;

	for(unsigned int i = 1; i < (x.size()-1); ++i) {
		sphericalCoords[i] = TMath::ACos(x[i-1] / xSums[i-1]);
	}

	const double lastCosPhi = x[x.size()-2] / xSums[x.size()-2];
	if(x[x.size()-1] >= 0.) {
		sphericalCoords[x.size()-1] = TMath::ACos(lastCosPhi);
	} else {
		sphericalCoords[x.size()-1] = TMath::TwoPi() - TMath::ACos(lastCosPhi);
	}

	return sphericalCoords;

}


parameterSpace::parameterSpace(const ampIntegralMatrix& integralMatrix)
	: _nEvents(integralMatrix.nmbEvents()),
	  _nDimensions(integralMatrix.nmbWaves()),
	  _phi((2 * _nDimensions) - 1, 0.),
	  _cosPhi((2 * _nDimensions) - 1, 0.),
	  _sinPhi((2 * _nDimensions) - 1, 0.),
	  _sigma(2 * _nDimensions, 0.),
	  _DRDPhi((2 * _nDimensions) - 1, 0.),
	  _DSigmaDPhi(2 * _nDimensions, vector<double>((2 * _nDimensions) - 1, 0.)),
	  _randomNumberGen(new TRandom3(123456789))
{
	_integralMatrix.resize(_nDimensions, vector<complex<double> >(_nDimensions, complex<double>(0., 0.)));
	for(unsigned int i = 0; i < _nDimensions; ++i) {
		for(unsigned int j = 0; j < _nDimensions; ++j) {
			_integralMatrix[i][j] = integralMatrix.matrix()[i][j] / (double)_nEvents;
		}
	}
}


parameterSpace::~parameterSpace() {
	delete _randomNumberGen;
}


void parameterSpace::pickUniformAngles() {

	vector<double> x((2*_nDimensions), 0.);
	vector<double> xSums(2*_nDimensions, 0.);

	double r = 0.;
	for(int i = (2*_nDimensions)-1; i >= 0; --i) {
		unsigned int j = (unsigned int)i;
		double xi = _randomNumberGen->Gaus();
		x[j] = xi;
		r += xi * xi;
		xSums[j] = sqrt(r);
	}
	r = sqrt(r);

	for(unsigned int i = 0; i < ((2*_nDimensions)-2); ++i) {
		_cosPhi[i] = x[i] / xSums[i];
		_sinPhi[i] = sqrt(1. - _cosPhi[i]*_cosPhi[i]);
		_phi[i] = TMath::ACos(_cosPhi[i]);
	}

	unsigned int lastPhiIndex = (2*_nDimensions)-2;
	_cosPhi[lastPhiIndex] = x[lastPhiIndex] / xSums[lastPhiIndex];
	if(x[(2*_nDimensions) - 1] >= 0.) {
		_phi[lastPhiIndex] = TMath::ACos(_cosPhi[lastPhiIndex]);
		_sinPhi[lastPhiIndex] = sqrt(1. - _cosPhi[lastPhiIndex]*_cosPhi[lastPhiIndex]);
	} else {
		_phi[lastPhiIndex] = TMath::TwoPi() - TMath::ACos(_cosPhi[lastPhiIndex]);
		_sinPhi[lastPhiIndex] = -sqrt(1. - _cosPhi[lastPhiIndex]*_cosPhi[lastPhiIndex]);
	}

	for(unsigned int i = 0; i < (2*_nDimensions); ++i) {
		_sigma[i] = x[i] / r;
	}

	initialize();

}


void parameterSpace::setAngles(const vector<double>& phi) {

	if(phi.size() != _phi.size()) {
		printErr << "Size mismatch between number of phis provided ("
		         << phi.size() << ") and number of phis needed ("
		         << _phi.size() << "). Aborting..." << endl;
		throw;
	}
	for(unsigned int i = 0; i < _phi.size(); ++i) {
		_phi[i] = phi[i];
		_cosPhi[i] = TMath::Cos(_phi[i]);
		_sinPhi[i] = TMath::Sin(_phi[i]);
		_sigma[i] = calculateSigma(i);
	}

	initialize();

}


double parameterSpace::getDS(const unsigned int& phiIndex) const
{
	double dS = 0.;
	const double R = getR();
	const double dRDPhi = getDRDPhi(phiIndex);
	for(unsigned int i = 0; i < (2*_nDimensions); ++i) {
		double x = (dRDPhi* _sigma[i]) + (R * getDSigmaDPhi(i, phiIndex));
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


double parameterSpace::getDAHat() const
{
	double dA = 1.;
	for(unsigned int i = 0; i < ((2*_nDimensions)-1); ++i) {
		dA *= getDSHat(i);
	}
	return dA;
}


double parameterSpace::getDASphereHat() const
{
	double dA = 1.;
	const unsigned int pow = (2*_nDimensions) - 2;
	for(unsigned int i = 0; i < pow; ++i) {
		for(unsigned int j = 0; j < (pow-i); ++j) {
			dA *= _sinPhi[i];
		}
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
	return -((sqrt(_nEvents * rho)) / (rho * rho)) * calculateDRhoDPhi(phiIndex);
}


void parameterSpace::calculateRho()
{

	complex<double> complexFirstSum(0., 0.);

	for(unsigned int alpha = 0; alpha < (2*_nDimensions); ++alpha) {
		unsigned int betaLowerLimit = (alpha < _nDimensions) ? 0 : _nDimensions;
		unsigned int betaUpperLimit = (betaLowerLimit == 0) ? _nDimensions : 2*_nDimensions;
		for(unsigned int beta = betaLowerLimit; beta < betaUpperLimit; ++beta) {
			const complex<double>& IAlphaBeta = _integralMatrix[alpha % _nDimensions][beta % _nDimensions];
			complexFirstSum += _sigma[alpha] * _sigma[beta] * IAlphaBeta;
		}
	}

	if(fabs(complexFirstSum.imag()) > 1e-10) {
		printErr << "Calculated rho with finite imaginary part (" << complexFirstSum.imag() << "). Aborting..." << endl;
		throw;
	}
	double firstSum = complexFirstSum.real();

	double secondSum = 0.;
	for(unsigned int alpha = 0; alpha < _nDimensions; ++alpha) {
		for(unsigned int beta = 0; beta < _nDimensions; ++beta) {
			const double& lAlphaBeta = _integralMatrix[alpha][beta].imag();
			secondSum += _sigma[alpha+_nDimensions] * _sigma[beta] * lAlphaBeta;
		}
	}

	_rho = firstSum - 2. * secondSum;

}


double parameterSpace::calculateDRhoDPhi(const unsigned int& phiIndex) const
{

	double retval = 0.;

	for(unsigned int alpha = phiIndex; alpha < (2*_nDimensions); ++alpha) {
		unsigned int betaLowerLimit = (alpha < _nDimensions) ? 0 : _nDimensions;
		unsigned int betaUpperLimit = (betaLowerLimit == 0) ? _nDimensions : 2*_nDimensions;
		for(unsigned int beta = betaLowerLimit; beta < betaUpperLimit; ++beta) {
			const double& kAlphaBeta = _integralMatrix[alpha % _nDimensions][beta % _nDimensions].real();
			retval += kAlphaBeta * _sigma[beta] * getDSigmaDPhi(alpha, phiIndex);
		}
	}
	retval *= 2.;

	double dSigma2DPhi = 0.;
	for(unsigned int beta = 0; beta < _nDimensions; ++beta) {
		for(unsigned int alpha = 0; alpha < _nDimensions; ++alpha) {
			const double& lAlphaBeta = _integralMatrix[alpha][beta].imag();
			dSigma2DPhi += lAlphaBeta * ((_sigma[beta] * getDSigmaDPhi(alpha+_nDimensions, phiIndex)) +
			                              (_sigma[alpha+_nDimensions] * getDSigmaDPhi(beta, phiIndex)));
		}
	}
	retval -= 2. * dSigma2DPhi;

	return retval;

}


double parameterSpace::calculateSigma(const unsigned int& index) const
{
	double sigma = (index == ((2*_nDimensions)-1)) ? 1. : _cosPhi[index];
	for(unsigned int i = 0; i < index; ++i) {
		sigma *= _sinPhi[i];
	}
	return sigma;
}


void parameterSpace::calculateAllDSigmaDPhi() {
	for(unsigned int sigmaIndex = 0; sigmaIndex < (2*_nDimensions); ++sigmaIndex) {
		for(unsigned int phiIndex = 0; phiIndex < ((2*_nDimensions)-1); ++phiIndex) {
			_DSigmaDPhi[sigmaIndex][phiIndex] = calculateDSigmaDPhi(sigmaIndex, phiIndex);
		}
	}
}


double parameterSpace::calculateDSigmaDPhi(const unsigned int& sigmaIndex, const unsigned int& phiIndex) const
{
	double retval = 0.;
	if(phiIndex > sigmaIndex) {
		retval = 0.;
	} else if(phiIndex == sigmaIndex) {
		retval = -1.;
		for(unsigned int i = 0; i <= sigmaIndex; ++i) {
			retval *= _sinPhi[i];
		}
	} else {
		retval = (sigmaIndex == ((2*_nDimensions)-1)) ? 1. : _cosPhi[sigmaIndex];
		for(unsigned int i = 0; i < sigmaIndex; ++i) {
			if(i == phiIndex) {
				retval *= _cosPhi[i];
			} else {
				retval *= _sinPhi[i];
			}
		}
	}
	return retval;
}
