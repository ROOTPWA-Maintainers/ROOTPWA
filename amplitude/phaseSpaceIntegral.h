
#include<complex>

#include"massDependence.h"
#include"isobarDecayTopology.h"
#include"isobarDecayVertex.h"
#include"particle.h"

#ifndef PHASESPACEINTEGRAL_H
#define PHASESPACEINTEGRAL_H

class TFile;
class TTree;

namespace rpwa {

	class phaseSpaceIntegral {

	  public:

		static phaseSpaceIntegral* instance();

		void setDecay(const isobarDecayTopology& decay) { _decay = &decay; }

		std::complex<double> operator()(const isobarDecayVertex& vertex);

	  private:

		phaseSpaceIntegral() { }

		double dyn();
		double readIntegralValueFromTree(const double& M, TTree* tree) const;

		std::vector<std::string> getFilenames() const;
		void createIntegralFile(std::string filename) const;

		double phaseSpace1D(double* x, double* p) const;

		double evalInt(const double& M) const;

		static phaseSpaceIntegral* _instance;

		const isobarDecayVertex* _vertex;

		const isobarDecayTopology* _decay;

		const static std::string TREE_NAME;
		const static std::string DIRECTORY;
		const static int N_POINTS = 10000;
		const static double LOWER_BOUND = 0.41871054;
		const static double UPPER_BOUND = 4.5;

};

}

#endif
