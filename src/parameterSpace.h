#ifndef EVIDENCECALCULATOR_H
#define EVIDENCECALCULATOR_H

#include<cmath>
#include<vector>


namespace rpwa {

	class ampIntegralMatrix;

	class parameterSpace {

	  public:

		parameterSpace(const unsigned int& nEvents,
		                   const unsigned int& nDimensions,
		                   const ampIntegralMatrix& integralMatrix)
			: _nEvents(nEvents),
			  _nDimensions(nDimensions),
			  _integralMatrix(integralMatrix) { _phi.resize(2 * _nDimensions); }

		void pickAngles();
		void setAngles(const std::vector<double>& phi);

		double getR() const { return std::sqrt( _nEvents / getRho() ); }
		double getDRDPhi(const unsigned int& phiIndex) const { return _DRDPhi[phiIndex]; }

		double getDS(const unsigned int& phiIndex) const;
		double getDA() const;


	  private:

		void initialize() { calculateRho(); calculateAllDRDPhi(); }

		void calculateAllDRDPhi();
		double calculateDRDPhi(const unsigned int& phiIndex) const;

		void calculateRho();
		const double& getRho() const { return _rho; }

		double calculateDRhoDPhi(const unsigned int& phiIndex) const;

		double calculateSigma(const unsigned int& index) const;
		double getDSigmaDPhi(const unsigned int& sigmaIndex, const unsigned int& phiIndex) const;

		const unsigned int _nEvents;
		const unsigned int _nDimensions;
		const ampIntegralMatrix& _integralMatrix;

		std::vector<double> _phi;
		std::vector<double> _sigma;

		double _rho;
		std::vector<double> _DRDPhi;


	};

}

#endif
