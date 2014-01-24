
#include <boost/progress.hpp>

#include <TMath.h>

#include "ampIntegralMatrix.h"
#include "parameterSpace.h"
#include "reportingUtils.hpp"


using namespace boost;
using namespace rpwa;
using namespace std;


typedef boost::multi_array<std::complex<double>, 2> integralMatrixType;

double surSphere(double r, double n) {
	return ((2*(pow(TMath::Pi(), n / 2.)) / TMath::Gamma(n / 2.)) * pow(r, n-1));
}

double logSurSphere(double r, double n) {
	return log(2.) + (n / 2.)*log(TMath::Pi()) - TMath::LnGamma(n/2.) + (n-1.)*log(r);
}


int main() {
/*
	const unsigned int N_DIMS = 30;
	const unsigned int N_EVENTS = 10000;
	ampIntegralMatrix integralMatrix;
	integralMatrixType& matrix = integralMatrix.matrix();
	matrix.resize(extents[N_DIMS][N_DIMS]);

	for(unsigned int i = 0; i < N_DIMS; ++i) {
		for(unsigned int j = 0; j < N_DIMS; ++j) {
			if(i == j) {
				matrix[i][j] = complex<double>(1., 0.);
			} else {
				matrix[i][j] = complex<double>(0., 0.);
			}
		}
	}
*/
	ampIntegralMatrix integralMatrix;
	integralMatrix.readAscii("/home/kbicker/analysis/monteCarlo/run7/pwaBins7.2/2680.2710/ACCAMPS/norm.int");
	const unsigned int N_DIMS = integralMatrix.nmbWaves();
	const unsigned int N_EVENTS = integralMatrix.nmbEvents();

	std::cout << "N_DIMS=" << N_DIMS << std::endl;
	std::cout << "N_EVENTS=" << N_EVENTS << std::endl;

	parameterSpace paramSpace(N_EVENTS, N_DIMS, integralMatrix);

	std::vector<double> dAs;
	double integral = 0.;
	const unsigned int NSAMPLES = 1000;
	progress_display progressIndicator(NSAMPLES, cout, "");
	for(unsigned int i = 0; i < NSAMPLES; ++i) {
		paramSpace.pickUniformAngles();
/*		std::cout<<"----------------"<<std::endl;
		std::cout<<"R="<<paramSpace.getR()<<" rho="<<paramSpace.getRho()<<" dA="<<paramSpace.getDARatio()<<std::endl;
		std::cout<<"dA="<<paramSpace.getDA()<<" dA^="<<paramSpace.getDAHat()<<" dA_sphere^="<<paramSpace.getDASphereHat()<<std::endl;
		std::cout<<"----------------"<<std::endl;
*/		integral += paramSpace.getDARatio();
		++progressIndicator;
	}

	integral /= (double)NSAMPLES;
	std::cout<<"integral="<<integral<<std::endl;
//	double V_S = surSphere(sqrt((double)N_EVENTS), 2. * (double)N_DIMS);
	double logV_S = logSurSphere(sqrt((double)N_EVENTS), 2. * (double)N_DIMS);
	std::cout<<"log(V_P)="<<log(integral) + logV_S<<std::endl;
	std::cout<<"log(V_S)="<<logV_S<<std::endl;

	return 0;

}
