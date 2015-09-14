
#include "randomNumberGenerator.h"


using namespace rpwa;


randomNumberGenerator* randomNumberGenerator::_randomNumberGenerator = 0;


randomNumberGenerator* randomNumberGenerator::instance() {
	if(not _randomNumberGenerator) {
		_randomNumberGenerator = new randomNumberGenerator();
	}
	return _randomNumberGenerator;
}


double randomNumberGenerator::rndm() {
	return _rndGen.Rndm();
}
