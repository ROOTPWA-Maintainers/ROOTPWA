
#include <TFile.h>
#include <TTree.h>

#include "generatorParameters.hpp"
#include "randomNumberGenerator.h"
#include "reportingUtils.hpp"
#include "beamAndVertexGenerator.h"


using namespace std;
using namespace rpwa;

beamAndVertexGenerator::beamAndVertexGenerator(string rootFileName,
                                   double massBeamParticle,
                                   Target target)
	: _rootFile(NULL),
	  _beamTree(NULL),
	  _vertex(TVector3(0., 0., 0.)),
	  _beam(TLorentzVector(0., 0., 0., 0.)),
	  _vertexX(pair<double, double>(0., 0.)),
	  _vertexY(pair<double, double>(0., 0.)),
	  _vertexZLow(target.position.Z() - 0.5 * target.length),
	  _targetLength(target.length),
	  _targetInteractionLength(target.interactionLength),
	  _beamMomentumX(pair<double, double>(0., 0.)),
	  _beamMomentumY(pair<double, double>(0., 0.)),
	  _beamMomentumZ(pair<double, double>(0., 0.)),
	  _massBeamParticle(massBeamParticle),
	  _sigmasPresent(false)
{
	_rootFile = TFile::Open(rootFileName.c_str(), "READ");
	if(not _rootFile) {
		printErr << "Could not open root file '" << rootFileName
		         << "' when initializing beam simulation." << endl;
		return;
	}
	const string beamTreeName = "beamTree";
	_beamTree = dynamic_cast<TTree*>(_rootFile->Get(beamTreeName.c_str()));
	if(not _beamTree) {
		printErr << "Could not find '" << beamTreeName << "' in root file '" << rootFileName
		         << "' when initializing beam simulation." << endl;
		return;
	}
	TBranch* vertexXPosition = _beamTree->FindBranch("vertex_x_position");
	TBranch* vertexYPosition = _beamTree->FindBranch("vertex_y_position");
	TBranch* beamMomentumX = _beamTree->FindBranch("beam_momentum_x");
	TBranch* beamMomentumY = _beamTree->FindBranch("beam_momentum_y");
	TBranch* beamMomentumZ = _beamTree->FindBranch("beam_momentum_z");
	if(not (vertexXPosition and vertexYPosition and
	        beamMomentumX and beamMomentumY and beamMomentumZ))
	{
		printErr << "One of the required branches is missing in tree '" << beamTreeName
		         << "' in root file '" << rootFileName << ". Required are "
		         << " 'vertex_x_position', 'vertex_y_position', 'beam_momentum_x'"
		         << ", 'beam_momentum_y' 'and beam_momentum_z'." << endl;
		_rootFile->Close();
		delete _rootFile;
		_rootFile = NULL;
		delete _beamTree;
		_beamTree = NULL;
		return;
	}
	_beamTree->SetBranchAddress("vertex_x_position", &_vertexX.first);
	_beamTree->SetBranchAddress("vertex_y_position", &_vertexY.first);
	_beamTree->SetBranchAddress("beam_momentum_x", &_beamMomentumX.first);
	_beamTree->SetBranchAddress("beam_momentum_y", &_beamMomentumY.first);
	_beamTree->SetBranchAddress("beam_momentum_z", &_beamMomentumZ.first);
	TBranch* vertexXPositionSigma = _beamTree->FindBranch("vertex_x_position_sigma");
	TBranch* vertexYPositionSigma = _beamTree->FindBranch("vertex_y_position_sigma");
	TBranch* beamMomentumXSigma =  _beamTree->FindBranch("beam_momentum_x_sigma");
	TBranch* beamMomentumYSigma =  _beamTree->FindBranch("beam_momentum_y_sigma");
	TBranch* beamMomentumZSigma =  _beamTree->FindBranch("beam_momentum_z_sigma");
	if(not (vertexXPositionSigma and vertexYPositionSigma and
	        beamMomentumXSigma and beamMomentumYSigma and beamMomentumZSigma))
	{
		printWarn << "One of the optional branches is missing in tree '" << beamTreeName
		          << "' in root file '" << rootFileName << ". Required are"
		          << " 'vertex_x_position_sigma', 'vertex_y_position_sigma', 'beam_momentum_x_sigma'"
		          << ", 'beam_momentum_y_sigma' 'and beam_momentum_z_sigma'. The input events"
		          << " will be used as-is." << endl;
		_sigmasPresent = false;
	} else {
		_beamTree->SetBranchAddress("vertex_x_position_sigma", &_vertexX.second);
		_beamTree->SetBranchAddress("vertex_y_position_sigma", &_vertexY.second);
		_beamTree->SetBranchAddress("beam_momentum_x_sigma", &_beamMomentumX.second);
		_beamTree->SetBranchAddress("beam_momentum_y_sigma", &_beamMomentumY.second);
		_beamTree->SetBranchAddress("beam_momentum_z_sigma", &_beamMomentumZ.second);
	}
}


beamAndVertexGenerator::~beamAndVertexGenerator() {
	if(_rootFile) {
		_rootFile->Close();
		delete _rootFile;
	}
	if(_beamTree) {
		delete _beamTree;
	}
}


bool beamAndVertexGenerator::check() {
	if(_beamTree) {
		return true;
	} else {
		return false;
	}
}


bool beamAndVertexGenerator::event() {
	TRandom3* randomGen = randomNumberGenerator::instance()->getGenerator();
	long nEntries = _beamTree->GetEntries();
	long entry = (long)-randomGen->Uniform(-nEntries, 0); // because Uniform(a, b) is in ]a, b]
	_beamTree->GetEntry(entry);
	double z;
	do {
		z = randomGen->Uniform();
	} while (randomGen->Uniform() < z*_targetInteractionLength);
	z = _vertexZLow + z * _targetLength;
	if(_sigmasPresent) {
		_vertex.SetXYZ(randomGen->Gaus(_vertexX.first, _vertexX.second),
		               randomGen->Gaus(_vertexY.first, _vertexY.second),
		               z);
		_beam.SetXYZM(randomGen->Gaus(_beamMomentumX.first, _beamMomentumX.second),
		              randomGen->Gaus(_beamMomentumY.first, _beamMomentumY.second),
		              randomGen->Gaus(_beamMomentumZ.first, _beamMomentumZ.second),
		              _massBeamParticle);
	} else {
		_vertex.SetXYZ(_vertexX.first, _vertexY.first, z);
		_beam.SetXYZM(_beamMomentumX.first,
		              _beamMomentumY.first,
		              _beamMomentumZ.first,
		              _massBeamParticle);
	}
	return false;
}
