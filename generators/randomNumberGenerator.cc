
#include "randomNumberGenerator.h"


using namespace rpwa;


randomNumberGenerator* randomNumberGenerator::_randomNumberGenerator = NULL;


randomNumberGenerator* randomNumberGenerator::instance() {
	if(_randomNumberGenerator == NULL) {
		_randomNumberGenerator = new randomNumberGenerator();
	}
	return _randomNumberGenerator;
}


unsigned int randomNumberGenerator::seed() {
	return _rndGen.GetSeed();
}


void randomNumberGenerator::setSeed(unsigned int seed) {
	_rndGen.SetSeed(seed);
}


double randomNumberGenerator::rndm() {
	return _rndGen.Rndm();
}
