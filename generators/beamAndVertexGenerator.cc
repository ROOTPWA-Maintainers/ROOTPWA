
#include <TFile.h>
#include <TTree.h>

#include "generator.h"
#include "generatorParameters.hpp"
#include "randomNumberGenerator.h"
#include "reportingUtils.hpp"
#include "beamAndVertexGenerator.h"


using namespace std;
using namespace rpwa;


beamAndVertexGenerator::beamAndVertexGenerator()
	: _beamFileName(""),
	  _vertex(TVector3(0., 0., 0.)),
	  _beam(TLorentzVector(0., 0., 0., 0.)),
	  _readBeamfileSequentially(false),
	  _currentBeamfileEntry(0),
	  _simpleSimulation(true),
	  _rootFile(NULL),
	  _beamTree(NULL),
	  _vertexX(pair<double, double>(0., 0.)),
	  _vertexY(pair<double, double>(0., 0.)),
	  _beamMomentumX(pair<double, double>(0., 0.)),
	  _beamMomentumY(pair<double, double>(0., 0.)),
	  _beamMomentumZ(pair<double, double>(0., 0.)),
	  _sigmasPresent(false),
	  _sigmaScalingFactor(1.) { }


bool beamAndVertexGenerator::loadBeamFile(const string& beamFileName)
{
	_beamFileName = beamFileName;
	_rootFile = TFile::Open(_beamFileName.c_str(), "READ");
	if(not _rootFile) {
		printErr << "Could not open root file '" << _beamFileName
		         << "' when initializing beam simulation." << endl;
		return false;
	}
	const string beamTreeName = "beamTree";
	_beamTree = dynamic_cast<TTree*>(_rootFile->Get(beamTreeName.c_str()));
	if(not _beamTree) {
		printErr << "Could not find '" << beamTreeName << "' in root file '" << _beamFileName
		         << "' when initializing beam simulation." << endl;
		return false;
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
		         << "' in root file '" << _beamFileName << ". Required are "
		         << " 'vertex_x_position', 'vertex_y_position', 'beam_momentum_x'"
		         << ", 'beam_momentum_y' 'and beam_momentum_z'." << endl;
		_rootFile->Close();
		delete _rootFile;
		_rootFile = NULL;
		delete _beamTree;
		_beamTree = NULL;
		return false;
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
		          << "' in root file '" << _beamFileName << ". Required are"
		          << " 'vertex_x_position_sigma', 'vertex_y_position_sigma', 'beam_momentum_x_sigma'"
		          << ", 'beam_momentum_y_sigma' 'and beam_momentum_z_sigma'. The input events"
		          << " will be used as-is." << endl;
		_sigmasPresent = false;
	} else {
		_sigmasPresent = true;
		_beamTree->SetBranchAddress("vertex_x_position_sigma", &_vertexX.second);
		_beamTree->SetBranchAddress("vertex_y_position_sigma", &_vertexY.second);
		_beamTree->SetBranchAddress("beam_momentum_x_sigma", &_beamMomentumX.second);
		_beamTree->SetBranchAddress("beam_momentum_y_sigma", &_beamMomentumY.second);
		_beamTree->SetBranchAddress("beam_momentum_z_sigma", &_beamMomentumZ.second);
	}
	_simpleSimulation = false;
	return check();
}


beamAndVertexGenerator::~beamAndVertexGenerator() {
	if(_rootFile) {
		_rootFile->Close();
	}
}


bool beamAndVertexGenerator::check() {
	if(_beamTree or _simpleSimulation) {
		return true;
	} else {
		return false;
	}
}


bool beamAndVertexGenerator::event(const Target& target, const Beam& beam) {
	double z = getVertexZ(target);
	TRandom3* randomGen = randomNumberGenerator::instance()->getGenerator();
	if(_simpleSimulation) {
		double x;
		double y;
		double radius = std::sqrt(randomGen->Uniform(0, target.radius * target.radius));
		randomGen->Circle(x, y, radius);
		_vertex.SetXYZ(target.position.X() + x,
		               target.position.Y() + y,
		               z);
		// throw magnituide of beam momentum
		const double pBeam = randomGen->Gaus(beam.momentum, beam.momentumSigma);
		// throw beam inclination
		const double dxdz = randomGen->Gaus(beam.DxDz, beam.DxDzSigma);
		const double dydz = randomGen->Gaus(beam.DyDz, beam.DyDzSigma);
		// construct tilted beam momentum Lorentz vector
		const double pz        = pBeam / sqrt(1 + dxdz * dxdz + dydz * dydz);
		const double px        = dxdz * pz;
		const double py        = dydz * pz;
		const double EBeam     = sqrt(pBeam * pBeam + beam.particle.mass2());
		_beam.SetXYZT(px, py, pz, EBeam);
	} else {
		long nEntries = _beamTree->GetEntries();
		if(not _readBeamfileSequentially) {
			_currentBeamfileEntry = (long)-randomGen->Uniform(-nEntries, 0); // because Uniform(a, b) is in ]a, b]
			_beamTree->GetEntry(_currentBeamfileEntry);
		} else {
			if(_currentBeamfileEntry >= nEntries) {
				printInfo << "reached end of beamfile, looping back to first event." << endl;
				_currentBeamfileEntry = 0;
			}
			_beamTree->GetEntry(_currentBeamfileEntry++);
		}
		if(_sigmasPresent and _sigmaScalingFactor != 0.) {
			_vertex.SetXYZ(randomGen->Gaus(_vertexX.first, _sigmaScalingFactor * _vertexX.second),
			               randomGen->Gaus(_vertexY.first, _sigmaScalingFactor * _vertexY.second),
			               z);
			_beam.SetXYZM(randomGen->Gaus(_beamMomentumX.first, _sigmaScalingFactor * _beamMomentumX.second),
			              randomGen->Gaus(_beamMomentumY.first, _sigmaScalingFactor * _beamMomentumY.second),
			              randomGen->Gaus(_beamMomentumZ.first, _sigmaScalingFactor * _beamMomentumZ.second),
						  beam.particle.mass());
		} else {
			_vertex.SetXYZ(_vertexX.first, _vertexY.first, z);
			_beam.SetXYZM(_beamMomentumX.first,
			              _beamMomentumY.first,
			              _beamMomentumZ.first,
			              beam.particle.mass());
		}
	}
	return true;
}


double beamAndVertexGenerator::getVertexZ(const Target& target) const {
	TRandom3* randomGen = randomNumberGenerator::instance()->getGenerator();
	double z;
	do {
		z = randomGen->Exp(target.interactionLength);
	} while(z > target.length);
	z = (target.position.Z() - target.length * 0.5) + z;
	return z;
}


ostream& beamAndVertexGenerator::print(ostream& out) {

	out << "Beam package information:" << endl;
	out << "    Using beam file ................ ";
	if(_simpleSimulation) {
		out << "No" << endl;
	} else {
		out << "Yes" << endl;
	}
	out << "    Beam file path ................. " << _beamFileName << endl;
	out << "    Beam file loaded ............... ";
	if(_rootFile) {
		out << "Yes" << endl;
	} else {
		out << "No" << endl;
	}
	out << "    Beam tree loaded ............... ";
	if(_beamTree) {
		out << "Yes" << endl;
	} else {
		out << "No" << endl;
	}
	out << "    Sigmas found ................... ";
	if(_sigmasPresent) {
		out << "Yes" << endl;
	} else {
		out << "No" << endl;
	}
	out << "    Sigma scaling factor ........... " << _sigmaScalingFactor << endl;
	out << "    Read beam file sequentially .... ";
	if(_readBeamfileSequentially) {
		out << "Yes" << endl;
	} else {
		out << "No" << endl;
	}
	return out;

}
