
#include <TFile.h>

#include "randomNumberGenerator.h"
#include "reportingUtils.hpp"
#include "primaryVertexGen.h"


using namespace std;
using namespace rpwa;

beamAndVertexGenerator::beamAndVertexGenerator(string rootFileName,
                                   double massBeamParticle)
	: _rootFile(0),
	  _beamTree(0),
	  _vertex(TVector3(0., 0., 0.)),
	  _beam(TLorentzVector(0., 0., 0., 0.)),
	  _massBeamParticle(massBeamParticle),
	  _sigmasPresent(false)
{
	_rootFile = TFile::Open(rootFileName.c_str(), "READ");
	if(not _rootFile) {
		printErr << "Could not open root file '" << rootFileName << "' when initializing beam simulation." << endl;
	}
}


beamAndVertexGenerator::~beamAndVertexGenerator() {
	_rootFile->Close();
}


bool beamAndVertexGenerator::check() {
	if(_beamTree) {
		return true;
	} else {
		return false;
	}
}


bool beamAndVertexGenerator::event() {
	return false;
}
