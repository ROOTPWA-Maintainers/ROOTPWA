
#include <TFile.h>
#include <TTree.h>
#include <TVectorD.h>
#include <TMatrixDSymEigen.h>

#include "generator.h"
#include "generatorParameters.hpp"
#include "randomNumberGenerator.h"
#include "reportingUtils.hpp"
#include "beamAndVertexGenerator.h"

namespace {
	TVectorD multivariateGauss(const TVectorD& mean, const TMatrixDSym& cov, const double sigmaScalingFactor=1.0);
}

using namespace std;
using namespace rpwa;


beamAndVertexGenerator::beamAndVertexGenerator()
	: _beamFileName(""),
	  _vertex(TVector3(0., 0., 0.)),
	  _beam(TLorentzVector(0., 0., 0., 0.)),
	  _readBeamfileSequentially(false),
	  _currentBeamfileEntry(0),
	  _simulationMode(simpleSimulation),
	  _rootFile(NULL),
	  _beamTree(NULL),
	  _vertexX(0.0),
	  _vertexY(0.0),
	  _vertexZ(0.0),
	  _beamMomentumX(0.0),
	  _beamMomentumY(0.0),
	  _beamMomentumZ(0.0),
	  _beamdXdZ(0.0),
	  _beamdYdZ(0.0),
	  _beamMomentum(0.0),
	  _vertexXSigma(0.0),
	  _vertexYSigma(0.0),
	  _beamMomentumXSigma(0.0),
	  _beamMomentumYSigma(0.0),
	  _beamMomentumZSigma(0.0),
	  _beamMomentumSigma(0.0),
	  _vertexCovarianceMatrix(nullptr),
	  _sigmasPresent(false),
	  _sigmaScalingFactor(1.),
	  _takeZpositionFromData(false),
	  _momentumResolution(0.0),
	  _momentumPDFSigma(0.0),
	  _momentumPDFMean(0.0){ }


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
	TBranch* beamMomentumX = _beamTree->FindBranch("beam_momentum_x");
	TBranch* beamdXdZ = _beamTree->FindBranch("beam_dXdZ");
	if (beamMomentumX != nullptr) {
		return loadBeamFileWithMomenta();
	} else if (beamdXdZ != nullptr) {
		return loadBeamFileWithInclinations();
	} else {
		printErr << "The beam file tree '" << _beamTree->GetName()
		         << "' contains neither momentum branches nor inclination branches." << endl;
	}
	return false;
}

bool beamAndVertexGenerator::loadBeamFileWithMomenta()
{

	// set required branches
	TBranch* vertexXPosition = _beamTree->FindBranch("vertex_x_position");
	TBranch* vertexYPosition = _beamTree->FindBranch("vertex_y_position");
	TBranch* beamMomentumX = _beamTree->FindBranch("beam_momentum_x");
	TBranch* beamMomentumY = _beamTree->FindBranch("beam_momentum_y");
	TBranch* beamMomentumZ = _beamTree->FindBranch("beam_momentum_z");
	if(not (vertexXPosition and vertexYPosition and
	        beamMomentumX and beamMomentumY and beamMomentumZ))
	{
		printErr << "One of the required branches is missing in tree '" << _beamTree->GetName()
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
	_beamTree->SetBranchAddress("vertex_x_position", &_vertexX);
	_beamTree->SetBranchAddress("vertex_y_position", &_vertexY);
	_beamTree->SetBranchAddress("beam_momentum_x", &_beamMomentumX);
	_beamTree->SetBranchAddress("beam_momentum_y", &_beamMomentumY);
	_beamTree->SetBranchAddress("beam_momentum_z", &_beamMomentumZ);

	// set optional branch with vertex z-position
	TBranch* vertexZPosition = _beamTree->FindBranch("vertex_z_position");
	if(vertexZPosition){
		_beamTree->SetBranchAddress("vertex_z_position", &_vertexZ);
	} else {
		_vertexZ = 0.0; // per default the beam parameters are given at Z=0
	}

	// set optional branches for Gaussian smearing
	TBranch* vertexXPositionSigma = _beamTree->FindBranch("vertex_x_position_sigma");
	TBranch* vertexYPositionSigma = _beamTree->FindBranch("vertex_y_position_sigma");
	TBranch* beamMomentumXSigma =  _beamTree->FindBranch("beam_momentum_x_sigma");
	TBranch* beamMomentumYSigma =  _beamTree->FindBranch("beam_momentum_y_sigma");
	TBranch* beamMomentumZSigma =  _beamTree->FindBranch("beam_momentum_z_sigma");
	if(not (vertexXPositionSigma and vertexYPositionSigma and
	        beamMomentumXSigma and beamMomentumYSigma and beamMomentumZSigma))
	{
		printWarn << "One of the optional branches for smearing the generated events with Gaussian uncertainties is missing in tree '"
		          << _beamTree->GetName() << "' in root file '" << _beamFileName << ". Required branches for smearing are"
		          << " 'vertex_x_position_sigma', 'vertex_y_position_sigma', 'beam_momentum_x_sigma'"
		          << ", 'beam_momentum_y_sigma' 'and beam_momentum_z_sigma'. The input events"
		          << " will be used without smearing." << endl;
		_sigmasPresent = false;
	} else {
		_sigmasPresent = true;
		_beamTree->SetBranchAddress("vertex_x_position_sigma", &_vertexXSigma);
		_beamTree->SetBranchAddress("vertex_y_position_sigma", &_vertexYSigma);
		_beamTree->SetBranchAddress("beam_momentum_x_sigma", &_beamMomentumXSigma);
		_beamTree->SetBranchAddress("beam_momentum_y_sigma", &_beamMomentumYSigma);
		_beamTree->SetBranchAddress("beam_momentum_z_sigma", &_beamMomentumZSigma);
	}

	if(_momentumPDFMean == 0.0){ // no momentum resolution correction
		_simulationMode = simulationFromMomenta;
	} else {
		_simulationMode = simulationFromMomentaCorrectMomentumResolution;
	}

	return check();
}

bool beamAndVertexGenerator::loadBeamFileWithInclinations()
{
	TBranch* vertexXPosition = _beamTree->FindBranch("vertex_x_position");
	TBranch* vertexYPosition = _beamTree->FindBranch("vertex_y_position");
	TBranch* vertexZPosition = _beamTree->FindBranch("vertex_z_position");
	TBranch* beamdXdZ = _beamTree->FindBranch("beam_dXdZ");
	TBranch* beamdYdZ = _beamTree->FindBranch("beam_dYdZ");
	TBranch* beamMomentum = _beamTree->FindBranch("beam_momentum");
	TBranch* vertexCov = _beamTree->FindBranch("vertex_XYdXdY_cov");
	TBranch* beamMomentumSigma = _beamTree->FindBranch("beam_momentum_sigma");
	if(not (vertexXPosition and vertexYPosition and vertexZPosition and
	        beamdXdZ and beamdYdZ and beamMomentum and vertexCov and beamMomentumSigma))
	{
		printErr << "One of the required branches is missing in tree '" << _beamTree->GetName()
		         << "' in root file '" << _beamFileName << ". Required are "
		         << " 'vertex_x_position', 'vertex_y_position', 'vertex_z_position', 'beam_dYdZ'"
		         << ", 'beam_dYdZ', 'beam_momentum' 'and vertex_XYdXdY_cov'." << endl;
		_rootFile->Close();
		delete _rootFile;
		_rootFile = nullptr;
		delete _beamTree;
		_beamTree = nullptr;
		return false;
	}
	_vertexCovarianceMatrix = nullptr;
	_beamTree->SetBranchAddress("vertex_x_position", &_vertexX);
	_beamTree->SetBranchAddress("vertex_y_position", &_vertexY);
	_beamTree->SetBranchAddress("vertex_z_position", &_vertexZ);
	_beamTree->SetBranchAddress("beam_dXdZ", &_beamdXdZ);
	_beamTree->SetBranchAddress("beam_dYdZ", &_beamdYdZ);
	_beamTree->SetBranchAddress("beam_momentum", &_beamMomentum);
	_beamTree->SetBranchAddress("vertex_XYdXdY_cov", &_vertexCovarianceMatrix);
	_beamTree->SetBranchAddress("beam_momentum_sigma", &_beamMomentumSigma);
	_simulationMode = simulationFromInclinations;
	_sigmasPresent = true;
	return check();
}


beamAndVertexGenerator::~beamAndVertexGenerator() {
	if(_rootFile) {
		_rootFile->Close();
	}
}


void beamAndVertexGenerator::randomizeBeamfileStartingPosition() {

	if(not _readBeamfileSequentially) {
		printWarn << "randomizing beamfile starting position without "
		          << "sequential beamfile reading has no effect." << endl;
	}
	TRandom3* randomGen = randomNumberGenerator::instance()->getGenerator();
	_currentBeamfileEntry = (long)-randomGen->Uniform(-_beamTree->GetEntries(), 0);

}


bool beamAndVertexGenerator::check() const {
	if( not (_beamTree or _simulationMode == simpleSimulation))  return false;
	if (_momentumPDFMean != 0.0){// do momentum resolution correction
		if (_simulationMode != simulationFromMomentaCorrectMomentumResolution){
			printErr << "Wrong simulation mode for for correcting the resolution in the beam momentum" << endl;
			return false;
		}
		if (_momentumPDFSigma == 0.0 or _momentumResolution == 0.0){
			printErr << "momentumResolution, momentumPDFSigma, and momentumPDFMean must be defined to use resolution correction" << endl;
			return false;
		}
	}
	return true;
}


bool beamAndVertexGenerator::event(const Target& target, const Beam& beam) {
	switch(_simulationMode){
		case simpleSimulation:
			return eventSimple(target, beam);
			break;
		case simulationFromMomenta:
			return eventFromMomenta(target, beam);
			break;
		case simulationFromMomentaCorrectMomentumResolution:
			return eventFromMomentaCorrectMomentumResolution(target, beam);
			break;
		case simulationFromInclinations:
			return eventFromInclinations(target, beam);
			break;
		default:
			printErr << "Unknown simulation mode" << endl;
			return false;
			break;
	}
}


bool beamAndVertexGenerator::eventSimple(const Target& target, const Beam& beam) {
	double z = getVertexZ(target);
	TRandom3* randomGen = randomNumberGenerator::instance()->getGenerator();
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
	return true;
}

bool
beamAndVertexGenerator::eventFromMomenta(const Target& target, const Beam& beam) {
	TRandom3* randomGen = randomNumberGenerator::instance()->getGenerator();
	if(not loadNextEventFromBeamfile()) return false;

	// either use z-position from data or draw random z-position according to the target definition
	const double z = getVertexZ(target); // has to be done AFTER loading the next event

	const double dxdz = _beamMomentumX / _beamMomentumZ;
	const double dydz = _beamMomentumY / _beamMomentumZ;
	const double projectedVertexX = _vertexX + dxdz * (z-_vertexZ);
	const double projectedVertexY = _vertexY + dydz * (z-_vertexZ);
	if(_sigmasPresent and _sigmaScalingFactor != 0.) {
		_vertex.SetXYZ(randomGen->Gaus(projectedVertexX, _sigmaScalingFactor * _vertexXSigma),
						randomGen->Gaus(projectedVertexY, _sigmaScalingFactor * _vertexYSigma),
						z);
		_beam.SetXYZM(randomGen->Gaus(_beamMomentumX, _sigmaScalingFactor * _beamMomentumXSigma),
						randomGen->Gaus(_beamMomentumY, _sigmaScalingFactor * _beamMomentumYSigma),
						randomGen->Gaus(_beamMomentumZ, _sigmaScalingFactor * _beamMomentumZSigma),
						beam.particle.mass());
	} else {
		_vertex.SetXYZ(projectedVertexX, projectedVertexY, z);
		_beam.SetXYZM(_beamMomentumX,
						_beamMomentumY,
						_beamMomentumZ,
						beam.particle.mass());
	}
	return true;
}


bool
beamAndVertexGenerator::eventFromMomentaCorrectMomentumResolution(const Target& target, const Beam& beam) {
	TRandom3* randomGen = randomNumberGenerator::instance()->getGenerator();
	if(not loadNextEventFromBeamfile()) return false;

	// either use z-position from data or draw random z-position according to the target definition
	const double z = getVertexZ(target); // has to be done AFTER loading the next event

	const double dxdz = _beamMomentumX / _beamMomentumZ;
	const double dydz = _beamMomentumY / _beamMomentumZ;
	const double projectedVertexX = _vertexX + dxdz * (z-_vertexZ);
	const double projectedVertexY = _vertexY + dydz * (z-_vertexZ);

	_beamMomentum = sqrt(_beamMomentumX*_beamMomentumX + _beamMomentumY*_beamMomentumY + _beamMomentumZ*_beamMomentumZ);
	const double truthMomentumPDFSigma = 1.0/sqrt(1.0/_momentumResolution/_momentumResolution+ 1.0/_momentumPDFSigma/_momentumPDFSigma);
	const double truthMomentumPDFMean = truthMomentumPDFSigma * truthMomentumPDFSigma * (_beamMomentum / (_momentumResolution * _momentumResolution)
	                                                                                     + _momentumPDFMean / (_momentumPDFSigma * _momentumPDFSigma));
	const double truthMomentum = randomGen->Gaus(truthMomentumPDFMean, truthMomentumPDFSigma);

	_vertex.SetXYZ(projectedVertexX, projectedVertexY, z);
	const double norm = sqrt(dxdz*dxdz + dydz*dydz + 1);
	_beam.SetXYZM(truthMomentum * dxdz / norm,
	              truthMomentum * dydz / norm,
	              truthMomentum * 1 / norm,
	              beam.particle.mass());
	return true;
}


bool
beamAndVertexGenerator::eventFromInclinations(const Target& target, const Beam& beam) {
	TRandom3* randomGen = randomNumberGenerator::instance()->getGenerator();
	if(not loadNextEventFromBeamfile()) return false;

	// either use z-position from data or draw random z-position according to the target definition
	const double z = getVertexZ(target); // has to be done AFTER loading the next event

	const bool smearing = _sigmasPresent and _sigmaScalingFactor != 0.0;

	TVectorD mean(4);
	mean[0] = _vertexX;
	mean[1] = _vertexY;
	mean[2] = _beamdXdZ;
	mean[3] = _beamdYdZ;
	const TVectorD sample = (smearing)? multivariateGauss(mean, *_vertexCovarianceMatrix, _sigmaScalingFactor) : mean;
	const double x = sample[0];
	const double y = sample[1];
	const double dxdz = sample[2];
	const double dydz = sample[3];
	const double norm = sqrt(dxdz*dxdz + dydz*dydz + 1);

	const double projectedVertexX = x + dxdz * (z-_vertexZ);
	const double projectedVertexY = y + dydz * (z-_vertexZ);
	const double beamMomentum = (smearing)? randomGen->Gaus(_beamMomentum, _sigmaScalingFactor*_beamMomentumSigma) : _beamMomentum;
		_vertex.SetXYZ(projectedVertexX, projectedVertexY, z);
		_beam.SetXYZM(  beamMomentum * dxdz / norm,
						beamMomentum * dydz / norm,
						beamMomentum * 1  / norm,
						beam.particle.mass());
	return true;
}


bool beamAndVertexGenerator::loadNextEventFromBeamfile(){
	const long nEntries = _beamTree->GetEntries();

	if(_readBeamfileSequentially) {
		if(_currentBeamfileEntry >= nEntries) {
			printInfo << "reached end of beamfile, looping back to first event." << endl;
			_currentBeamfileEntry = 0;
		}
		_beamTree->GetEntry(_currentBeamfileEntry++);
	} else {
		TRandom3* randomGen = randomNumberGenerator::instance()->getGenerator();
		_currentBeamfileEntry = (long)randomGen->Uniform(0,nEntries); // because Uniform(a, b) is in ]a, b]
		randomGen->Integer(15);
		_beamTree->GetEntry(_currentBeamfileEntry);
	}
	return true;
}

double beamAndVertexGenerator::getVertexZ(const Target& target) const {
	if (_takeZpositionFromData){
		return _vertexZ; // assumes that the next event is loaded from the tree before this function is called
	} else {
		TRandom3* randomGen = randomNumberGenerator::instance()->getGenerator();
		double z;
		do {
			z = randomGen->Exp(target.interactionLength);
		} while(z > target.length);
		z = (target.position.Z() - target.length * 0.5) + z;
		return z;
	}
}


void
beamAndVertexGenerator::setTakeZpositionFromData(const bool takeZpositionFromData)
{
	if(takeZpositionFromData and _beamFileName.size()==0){
		printErr << "Cannot take Z-positions for vertices from data if no beam file is given!" << endl;
		throw;
	}
	_takeZpositionFromData = takeZpositionFromData;
}

void
beamAndVertexGenerator::setMomentumResolutionCorrection(const double resolution, const double momentumPDFsigma, const double momentumPDFmean){
	_momentumResolution = resolution;
	_momentumPDFSigma = momentumPDFsigma;
	_momentumPDFMean = momentumPDFmean;
}


ostream& beamAndVertexGenerator::print(ostream& out) const
{

	out << "Beam package information:" << endl;
	out << "    Using beam file ................ ";
	out << yesNo(_simulationMode == simpleSimulation) << endl;
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
	out << "    Simulation mode ................ ";
	if(_simulationMode == simulationFromMomenta) {
		out << "Using momenta" << endl;
	} else if (_simulationMode == simulationFromMomentaCorrectMomentumResolution) {
		out << "Using momenta" << endl;
		out << "    Correction for beam resolution: res= ";
		out << _momentumResolution << ", PDFsigma= " << _momentumPDFSigma;
		out << ", PDFmean= " << _momentumPDFMean << endl;
	}else {
		out << "Using beam inclination" << endl;
	}
	out << "    Sigmas found ................... ";
	if(_sigmasPresent) {
		out << "Yes" << endl;
	} else {
		out << "No" << endl;
	}
	out << "    Sigma scaling factor ........... " << _sigmaScalingFactor << endl;
	out << "    Z-positions from beam file ..... ";
	if(_takeZpositionFromData) {
		out << "Yes" << endl;
	} else {
		out << "No" << endl;
	}
	out << "    Read beam file sequentially .... ";
	if(_readBeamfileSequentially) {
		out << "Yes" << endl;
		out << "    Beamfile starting position ..... " << _currentBeamfileEntry << endl;
	} else {
		out << "No" << endl;
	}
	return out;

}

namespace{
TVectorD multivariateGauss(const TVectorD& mean, const TMatrixDSym& cov, const double sigmaScalingFactor)
{
	TRandom3* randomGen = randomNumberGenerator::instance()->getGenerator();
	const int ndim = mean.GetNrows();
	if ( ndim != cov.GetNrows() or ndim != cov.GetNcols()){
		printErr << "Covariance matrix or mean vector have wrong format" << endl;
		throw "Covariance matrix or mean vector have wrong format";
	}

	TMatrixDSymEigen eigenDecomposition(cov);
	const TVectorD& eigenvalues = eigenDecomposition.GetEigenValues();
	const TMatrixD& eigenvectorMatrix = eigenDecomposition.GetEigenVectors();

	TVectorD sample(ndim);
	for(int i=0; i<ndim; ++i){
		sample[i] = randomGen->Gaus(0,sigmaScalingFactor*sqrt(eigenvalues[i]));
	}
	sample = eigenvectorMatrix * sample;
	sample = sample + mean;
	return sample;
}
}
