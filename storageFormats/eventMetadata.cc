
#include "eventMetadata.h"
#include "reportingUtils.hpp"


using namespace std;
using namespace rpwa;


rpwa::eventMetadata::eventMetadata()
	: _userString(""),
	  _contentHash(""),
	  _initialStateParticleNames(),
	  _finalStateParticleNames(),
	  _binningMap()
{ }


rpwa::eventMetadata::~eventMetadata() { };


ostream& rpwa::eventMetadata::print(ostream& out) const
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


void rpwa::eventMetadata::setInitialStateParticleNames(const vector<string>& initialStateParticleNames)
{
	_initialStateParticleNames = initialStateParticleNames;
}


void rpwa::eventMetadata::setFinalStateParticleNames(const vector<string>& finalStateParticleNames)
{
	_finalStateParticleNames = finalStateParticleNames;
}


void rpwa::eventMetadata::setBinningVariableLabels(const vector<string>& labels)
{
	for(unsigned int i = 0; i < labels.size(); ++i) {
		_binningMap[labels[i]] = rangePairType(0., 0.);
	}
}


void rpwa::eventMetadata::setBinningVariableRange(const string& label, const rangePairType& range)
{
	_binningMap[label] = range;
}


void rpwa::eventMetadata::setBinningMap(const binningMapType& binningMap)
{
	_binningMap = binningMap;
}
