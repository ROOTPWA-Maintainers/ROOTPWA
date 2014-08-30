
#include "dataMetadata.h"
#include "reportingUtils.hpp"


using namespace std;
using namespace rpwa;


rpwa::dataMetadata::dataMetadata()
	: _userString(""),
	  _contentHash(""),
	  _initialStateParticleNames(),
	  _finalStateParticleNames(),
	  _binningMap()
{ }


rpwa::dataMetadata::~dataMetadata() { };


ostream& rpwa::dataMetadata::print(ostream& out) const
{
	out << "dataMetadata: " << endl
	    << "    userString ...................... '" << _userString << "'"         << endl
	    << "    contentHash ..................... '" << _contentHash << "'"        << endl
	    << "    initial state particle names: ... "  << _initialStateParticleNames << endl
        << "    final state particle names: ..... "  << _finalStateParticleNames   << endl
	    << "    binning map";
	if(_binningMap.empty()) {
		out << " ..................... " << "<empty>" << endl;
	} else {
		out << ": " << endl;
		for(binningMapType::const_iterator it = _binningMap.begin(); it != _binningMap.end(); ++it) {
			out << "        variable '" << it->first << "' range " << it->second << endl;
		}
	}
	return out;
}


void rpwa::dataMetadata::setInitialStateParticleNames(const vector<string>& initialStateParticleNames)
{
	_initialStateParticleNames = initialStateParticleNames;
}


void rpwa::dataMetadata::setFinalStateParticleNames(const vector<string>& finalStateParticleNames)
{
	_finalStateParticleNames = finalStateParticleNames;
}


void rpwa::dataMetadata::setBinningVariableLabels(const vector<string>& labels)
{
	for(unsigned int i = 0; i < labels.size(); ++i) {
		_binningMap[labels[i]] = rangePairType(0., 0.);
	}
}


void rpwa::dataMetadata::setBinningVariableRange(const string& label, const rangePairType& range)
{
	_binningMap[label] = range;
}


void rpwa::dataMetadata::setBinningMap(const binningMapType& binningMap)
{
	_binningMap = binningMap;
}
