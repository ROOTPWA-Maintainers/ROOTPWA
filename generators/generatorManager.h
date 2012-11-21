#ifndef GENERATORMANAGER_HH_
#define GENERATORMANAGER_HH_


#include<string>

#include "generatorParameters.hpp"


class TVector3;

namespace rpwa {

	class generator;
	class massAndTPrimePicker;
	class primaryVertexGen;

	class generatorManager {

	  public:

		generatorManager();
		~generatorManager();

		bool readReactionFile(const std::string& fileName);

		bool initializeGenerator();

		static bool debug() { return _debug; };
		static void setDebug(bool debug = true) { _debug = debug; };

	  private:

		rpwa::Beam _beam;
		rpwa::Target _target;
		rpwa::FinalState _finalState;

		rpwa::primaryVertexGen* _primaryVertexGen;

		rpwa::massAndTPrimePicker* _pickerFunction;

		bool _reactionFileRead;

		rpwa::generator* _generator;

		static bool _debug;

	};

}

#endif

