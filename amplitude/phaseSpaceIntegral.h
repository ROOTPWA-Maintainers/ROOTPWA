
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

		std::complex<double> operator()(const isobarDecayVertex& vertex);

	  private:

		phaseSpaceIntegral() { }

		double dyn();
		double readIntegralValueFromTree(const double& M, TTree* tree) const;

		void createIntegralFile() const;

		double evalInt(const double& M, const unsigned int& nEvents) const;

		static phaseSpaceIntegral* _instance;

		isobarDecayVertexPtr _vertex;
		std::string _filename;
		isobarDecayTopologyPtr _subDecay;

		const static std::string TREE_NAME;
		const static std::string DIRECTORY;
		const static int N_POINTS = 50;
		const static int N_MC_EVENTS = 1000000;
		const static int N_MC_EVENTS_FOR_M0 = 10000000;
		const static int MC_SEED = 987654321;
		const static double LOWER_BOUND = 0.41871054;
		const static double UPPER_BOUND = 4.5;
		const static bool NEW_FILENAME_CONVENTION = false;

};

}

#endif
