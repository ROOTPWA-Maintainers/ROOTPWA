#ifndef GENERATORMANAGER_HH_
#define GENERATORMANAGER_HH_

#include<string>

#include "generatorParameters.hpp"

class TVector3;

namespace rpwa {

	class generatorManager {

	  public:

		generatorManager();
		~generatorManager();

		bool readReactionFile(const std::string& fileName);

		static bool debug() { return _debug; };
		static void setDebug(bool debug = true) { _debug = debug; };

	  private:

		rpwa::Beam _beam;
		rpwa::Target _target;
		rpwa::FinalState _finalState;

		static bool _debug;

	};

}

#endif

