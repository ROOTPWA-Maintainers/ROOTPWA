#ifndef PHASESPACEINTEGRAL_H
#define PHASESPACEINTEGRAL_H

#include<complex>

#include"isobarDecayVertex.h"
#include"isobarDecayTopology.h"
#include"particle.h"

class TFile;
class TTree;

namespace rpwa {

	class integralTablePoint {

	  public:

		integralTablePoint(double MM = 0., double intV = 0., double intErr = 0.)
			: M(MM),
			  integralValue(intV),
			  integralError(intErr) { };

		double M;
		double integralValue;
		double integralError;

	};

	class integralTableContainer {

	  public:

		integralTableContainer() : _init(false) { }
		integralTableContainer(const isobarDecayVertex& vertex);
		~integralTableContainer() { }

		std::complex<double> operator()(double M, double M0, double Gamma0);

		static std::string getSubWaveNameFromVertex(const isobarDecayVertex& vertex);
		static std::string getSubWaveNameFromVertex(const isobarDecayVertex& vertex,
		                                            isobarDecayVertexPtr& vertexPtr,
		                                            isobarDecayTopologyPtr& subDecay);

	  private:

		double dyn(double M, double M0);
		double interpolate(const double& M) const;
		double getInt0(const double& M0);

		void fillIntegralTable();
		void addToIntegralTable(const integralTablePoint& newPoint);

		void readIntegralFile();
		void writeIntegralTableToDisk(bool overwriteFile = false) const;

		integralTablePoint evalInt(const double& M, const unsigned int& nEvents) const;

		std::vector<integralTablePoint> _integralTable;
		std::vector<double> _M0s;

		isobarDecayVertexPtr _vertex;
		std::string _subWaveName;
		std::string _fullPathToFile;
		isobarDecayTopologyPtr _subDecay;

		bool _init;

		const static int N_POINTS = 50;
		const static int N_MC_EVENTS = 1000000;
		const static int N_MC_EVENTS_FOR_M0 = 10000000;
		const static int MC_SEED = 987654321;
		const static double UPPER_BOUND = 4.5;
		const static bool NEW_FILENAME_CONVENTION = false;
		const static bool CALCULATE_ERRORS = true;

		const static std::string DIRECTORY;
		const static std::string TREE_NAME;

	};


	class phaseSpaceIntegral {

	  public:

		static phaseSpaceIntegral* instance();
		std::complex<double> operator()(const isobarDecayVertex& vertex);

	  private:

		phaseSpaceIntegral() { };

		static phaseSpaceIntegral* _instance;
		std::map<const isobarDecayVertex*, std::string> _vertexToSubwaveName;
		std::map<std::string, integralTableContainer> _subwaveNameToIntegral;

	};

}

#endif
