
#include "phaseSpaceIntegral.h"

#include <boost/progress.hpp>
#include <boost/algorithm/string.hpp>

#include<TClonesArray.h>
#include<TDirectory.h>
#include<TF1.h>
#include<TFile.h>
#include<TGraphErrors.h>
#include<TObjString.h>
#include<TTree.h>

#include"factorial.hpp"
#include"isobarHelicityAmplitude.h"
#include"mathUtils.hpp"
#include"nBodyPhaseSpaceGen.h"
#include"physUtils.hpp"
#include"waveDescription.h"

using namespace std;
using namespace rpwa;

phaseSpaceIntegral* phaseSpaceIntegral::_instance = 0;

phaseSpaceIntegral* phaseSpaceIntegral::instance()
{
	if(not _instance) {
		printInfo << "singleton instantiated." << endl;
		_instance = new phaseSpaceIntegral();
	}
	return _instance;
}


complex<double> phaseSpaceIntegral::operator()(const isobarDecayVertex& vertex) {

	map<const isobarDecayVertex*, string>::iterator name_it = _vertexToSubwaveName.find(&vertex);
	if(name_it == _vertexToSubwaveName.end()) {
		string waveName = integralTableContainer::getSubWaveNameFromVertex(vertex);
		_vertexToSubwaveName[&vertex] = waveName;
	}
	const string& waveName = _vertexToSubwaveName[&vertex];
	map<string, integralTableContainer>::iterator cont_it = _subwaveNameToIntegral.find(waveName);
	if(cont_it == _subwaveNameToIntegral.end()) {
		printInfo << "adding new integralTableContainer for waveName=\"" << waveName << "\"." << endl;
		_subwaveNameToIntegral[waveName] = integralTableContainer(vertex);
	}

	// get Breit-Wigner parameters
	const particlePtr& parent = vertex.parent();
	const double       M      = parent->lzVec().M();             // parent mass
	const double       M0     = parent->mass();                  // resonance peak position
	const double       Gamma0 = parent->width();                 // resonance peak width

	return _subwaveNameToIntegral[waveName](M, M0, Gamma0);

}


void phaseSpaceIntegral::removeVertex(const isobarDecayVertex* vertex) {
	map<const isobarDecayVertex*, string>::iterator name_it = _vertexToSubwaveName.find(vertex);
	if(name_it != _vertexToSubwaveName.end()) {
		_vertexToSubwaveName.erase(name_it);
		printInfo << "removed vertex from map" << endl;
	}
}


namespace {

	// Re-Implementation of waveDescription::waveNameFromTopology which omits the
	// spin projection quantum number M for the top vertex, as it does not change
	// the integral.
	string
	__waveNameFromTopology(const isobarDecayTopology& topo)
	{
		ostringstream waveName;
		if (not topo.checkTopology() or not topo.checkConsistency()) {
			printWarn << "decay topology has issues. cannot construct wave name." << endl;
			return "";
		}
		// X quantum numbers
		const particle& X = *(topo.XParticle());
		waveName << "[" << spinQn(X.isospin()) << parityQn(X.G()) << ","
		         << spinQn(X.J()) << parityQn(X.P()) << parityQn(X.C()) << ","
		         << "M" << parityQn(X.reflectivity()) << "]"
		         << waveDescription::waveNameFromTopology(topo, topo.XIsobarDecayVertex());
		{
			string name = waveName.str();
			boost::replace_all(name, "(", "_");
			boost::replace_all(name, ")", "_");
			return name;
		}
	}
}


const string integralTableContainer::TREE_NAME = "psint";
const string integralTableContainer::DIRECTORY = "/home/kbicker/analysis/integralAmplitudesPwd";
const int integralTableContainer::N_POINTS = 50;
const int integralTableContainer::N_MC_EVENTS = 1000000;
const int integralTableContainer::N_MC_EVENTS_FOR_M0 = 10000000;
const int integralTableContainer::MC_SEED = 987654321;
const double integralTableContainer::UPPER_BOUND = 4.5;
const bool integralTableContainer::CALCULATE_ERRORS = true;


string integralTableContainer::getSubWaveNameFromVertex(const isobarDecayVertex& vertex)
{
	isobarDecayVertexPtr vertexPtr = isobarDecayVertexPtr();
	isobarDecayTopologyPtr subDecay = isobarDecayTopologyPtr();
	return integralTableContainer::getSubWaveNameFromVertex(vertex, vertexPtr, subDecay);
}


string integralTableContainer::getSubWaveNameFromVertex(const isobarDecayVertex& vertex,
                                                        isobarDecayVertexPtr& vertexPtr,
                                                        isobarDecayTopologyPtr& subDecay)
{

	vertexPtr = isobarDecayVertexPtr();
	const isobarDecayTopologyPtr mainTopo = boost::dynamic_pointer_cast<isobarDecayTopology>(vertex.decay());
	if(not mainTopo) {
		printErr << "got NULL pointer from vertex. Aborting..." << endl;
		throw;
	}
	const vector<interactionVertexPtr>& decayVertices = mainTopo->decayVertices();
	bool success = false;
	for(unsigned int i = 0; i < decayVertices.size(); ++i) {
		if(*(decayVertices[i]) == vertex) {
			vertexPtr = boost::dynamic_pointer_cast<isobarDecayVertex>(decayVertices[i]);
			success = true;
			break;
		}
	}
	if(not success) {
		printErr << "decay vertex of particle " << vertex.parent()->label() << " not found. Aborting..." << endl;
	}
	if(not vertexPtr) {
		printErr << "could not find decay vertex in decay topology. Aborting..." << endl;
		throw;
	}

	subDecay = createIsobarDecayTopology(mainTopo->subDecayConsistent(vertexPtr));
	return __waveNameFromTopology(*subDecay);

}


integralTableContainer::integralTableContainer(const isobarDecayVertex& vertex)
: _init(true)
{

	_subWaveName = getSubWaveNameFromVertex(vertex, _vertex, _subDecay);

	stringstream sstr;
	sstr<<DIRECTORY<<"/"<<__waveNameFromTopology(*_subDecay)<<".root";
	_fullPathToFile = sstr.str();

	readIntegralFile();

}


complex<double> integralTableContainer::operator()(double M, double M0, double Gamma0) {

	if(not _init) {
		printErr << "trying to use uninitialized integralTableContainer. Aborting..." << endl;
		throw;
	}

	const double Gamma = Gamma0 * dyn(M, M0);
	// A / (B - iC) = (A / (B^2 + C^2)) * (B + iC)
	const double A = M0 * Gamma0;
	const double B = M0 * M0 - M * M;
	const double C = M0 * Gamma;
	const complex<double> bw = (A / (B * B + C * C)) * complex<double>(B, C);

	return bw;
	// return (M0 * Gamma0) / (M0 * M0 - M * M - imag * M0 * Gamma);
}

double integralTableContainer::dyn(double M, double M0) {

	double psInt = interpolate(M);
	double psInt0 = getInt0(M0);
	return (psInt / psInt0);

}


double integralTableContainer::interpolate(const double& M) const
{

	double currentM = _integralTable[0].M;
	double psInt = _integralTable[0].integralValue;
	double lastM = currentM;
	double lastPsInt = psInt;
	for(unsigned int i = 1; i < _integralTable.size(); ++i) {
		currentM = _integralTable[i].M;
		psInt = _integralTable[i].integralValue;
		if((lastM <= M) and (M < currentM)) {
			break;
		}
		lastM = currentM;
		lastPsInt = psInt;
	}
	return (lastPsInt + (((psInt - lastPsInt) / (currentM - lastM)) * (M - lastM)));

}


double integralTableContainer::getInt0(const double& M0) {
	bool found = false;
	for(unsigned int i = 0; i < _M0s.size(); ++i) {
		if(fabs(_M0s[i] - M0) < 1e-10) {
			found = true;
			break;
		}
	}
	if(not found) {
		_M0s.push_back(M0);
	}
	for(unsigned int i = 0; i < _integralTable.size(); ++i) {
		if(fabs(_integralTable[i].M - M0) < 1e-10) {
			return interpolate(M0);
		}
	}
	printInfo << "adding new value for M0=" << M0 << " to integral table." << endl;
	addToIntegralTable(evalInt(M0, N_MC_EVENTS_FOR_M0));
	writeIntegralTableToDisk(true);
	return getInt0(M0);
}


integralTablePoint integralTableContainer::evalInt(const double& M, const unsigned int& nEvents) const {

	printInfo << "calculating integral for " << _subWaveName
	          << " at mother mass = " << M << "GeV with "
	          << nEvents << " events." << endl;

	const decayTopologyPtr origTopology = _vertex->decay();

	const massDependencePtr originalMassDep = _vertex->massDependence();
	_vertex->setMassDependence(createFlatMassDependence());

	const unsigned int nmbFsParticles = _subDecay->nmbFsParticles();

	isobarHelicityAmplitude subAmp = isobarHelicityAmplitude(_subDecay);
	TClonesArray prodNames("TObjString", 1);
	TClonesArray decayNames("TObjString", nmbFsParticles);
	new (prodNames[0]) TObjString(_vertex->parent()->name().c_str());
	vector<double> daughterMasses(nmbFsParticles, 0.);
	double fsParticlesMassSum = 0.;
	for(unsigned int i = 0; i < nmbFsParticles; ++i) {
		const particlePtr particle = _subDecay->fsParticles()[i];
		new (decayNames[i]) TObjString(particle->name().c_str());
		daughterMasses[i] = particle->mass();
		fsParticlesMassSum += daughterMasses[i];
	}
	subAmp.decayTopology()->initKinematicsData(prodNames, decayNames);
	subAmp.enableReflectivityBasis(false);
	subAmp.init();

	nBodyPhaseSpaceGen psGen;
	psGen.setDecay(daughterMasses);
	psGen.setMaxWeight(1.01 * psGen.estimateMaxWeight(M, 100000));
	TLorentzVector parent(0., 0., 0., M);

	double integral = 0.;
	vector<double> samples;
	if(CALCULATE_ERRORS) {
		samples.resize(nEvents, 0.);
	}

	boost::progress_display progressIndicator(nEvents, cout, "");
	for(unsigned int i = 0; i < nEvents; ++i) {
		psGen.generateDecay(parent);
		TClonesArray prodKinMom("TVector3", 1);
		TClonesArray decayKinMom("TVector3", nmbFsParticles);
		TLorentzVector particleSum(0., 0., 0., 0.);
		for(unsigned int j = 0; j < nmbFsParticles; ++j) {
			const TLorentzVector& daughter = psGen.daughter(j);
			new (decayKinMom[j]) TVector3(daughter.Vect());
			particleSum += daughter;
		}
		new (prodKinMom[0]) TVector3(particleSum.Vect());
		subAmp.decayTopology()->readKinematicsData(prodKinMom, decayKinMom);
		const double ampVal = norm(subAmp());
		const double sample = ampVal * psGen.eventWeight();
		integral += sample;
		if(CALCULATE_ERRORS) {
			samples[i] = sample;
		}
		++progressIndicator;
	}

	_vertex->setMassDependence(originalMassDep);
	const double V = fourPi * (1 / rpwa::factorial<double>(nmbFsParticles)) * pow((fourPi*(M - fsParticlesMassSum)), nmbFsParticles-2);
	integral /= (double)nEvents;

	double error = 0.;
	if(CALCULATE_ERRORS) {
		double sig2 = 0.;
		for(unsigned int i = 0; i < nEvents; ++i) {
			double diff = samples[i] - integral;
			sig2 += diff * diff;
		}
		sig2 /= (nEvents - 1);
		error = V * sqrt(sig2 / nEvents);
	}

	integral *= V;

	_subDecay->saveDecayToVertices(origTopology);

	printSucc << "calculated integral: " << integral << endl;
	return integralTablePoint(M, integral, error);
}


void integralTableContainer::readIntegralFile() {

	TDirectory* pwd = gDirectory;
	TFile* integralTableFile = TFile::Open(_fullPathToFile.c_str(), "READ");
	if(not integralTableFile) {
		printInfo << "no file '" << _fullPathToFile << "' with integral table found. Creating it..." << endl;
		fillIntegralTable();
		writeIntegralTableToDisk();
		printInfo << "created integral file '" << _fullPathToFile << "." << endl;
		integralTableFile = TFile::Open(_fullPathToFile.c_str(), "READ");
	} else {
		printInfo << "reading integral table from file '" << _fullPathToFile << "'." << endl;
	}
	TTree* tree = (TTree*)integralTableFile->Get(TREE_NAME.c_str());
	double M, psInt, psIntErr;
	tree->SetBranchAddress("M", &M);
	tree->SetBranchAddress("int", &psInt);
	tree->SetBranchAddress("err", &psIntErr);
	_integralTable.resize(tree->GetEntries());
	for(long i = 0; i < tree->GetEntries(); ++i) {
		tree->GetEntry(i);
		_integralTable[i].M = M;
		_integralTable[i].integralValue = psInt;
		_integralTable[i].integralError = psIntErr;
	}
	pwd->cd();

}


void integralTableContainer::writeIntegralTableToDisk(bool overwriteFile) const {

	TDirectory* pwd = gDirectory;
	TFile* integralFile = 0;
	if(overwriteFile) {
		integralFile = TFile::Open(_fullPathToFile.c_str(), "RECREATE");
	} else {
		integralFile = TFile::Open(_fullPathToFile.c_str(), "NEW");
	}
	integralFile->cd();
	TObjString filenameForRoot(_subWaveName.c_str());
	filenameForRoot.Write();
	TTree* outTree = new TTree(TREE_NAME.c_str(), TREE_NAME.c_str());
	double M, psInt, psIntErr;
	outTree->Branch("M", &M);
	outTree->Branch("int", &psInt);
	outTree->Branch("err", &psIntErr);

	const unsigned int nPoints = _integralTable.size();

	TGraphErrors* rawIntGraph = new TGraphErrors(nPoints);
	stringstream sstr;
	sstr<<_subWaveName<<"_integralValues";
	rawIntGraph->SetTitle(sstr.str().c_str());
	rawIntGraph->SetName(sstr.str().c_str());
	vector<TGraphErrors*> graphs(_M0s.size(), 0);
	for(unsigned int i = 0; i < graphs.size(); ++i) {
		graphs[i] = new TGraphErrors(nPoints);
		sstr.str("");
		sstr<<_subWaveName<<"_dyn_M0-"<<_M0s[i];
		graphs[i]->SetTitle(sstr.str().c_str());
		graphs[i]->SetName(sstr.str().c_str());
	}
	vector<integralTablePoint> psInt0s(_M0s.size());
	for(unsigned int j = 0; j < psInt0s.size(); ++j) {
		const double& M0 = _M0s[j];
		for(unsigned int i = 0; i < nPoints; ++i) {
			if(fabs(M0 - _integralTable[i].M) < 1e-10) {
				psInt0s[j] = _integralTable[i];
			}
		}
	}

	for(unsigned int i = 0; i < nPoints; ++i) {
		const integralTablePoint& point = _integralTable[i];
		M = point.M;
		psInt = point.integralValue;
		psIntErr = point.integralError;
		outTree->Fill();
		for(unsigned int j = 0; j < _M0s.size(); ++j) {
			const double dyn = psInt / psInt0s[j].integralValue;
			const double dynErr = ((psIntErr / psInt) + (psInt0s[j].integralError / psInt0s[j].integralValue)) * dyn;
			graphs[j]->SetPoint(i, M, dyn);
			graphs[j]->SetPointError(i, 0., dynErr);
		}
		rawIntGraph->SetPoint(i, M, psInt);
		rawIntGraph->SetPointError(i, 0., psIntErr);
	}

	integralFile->cd();
	rawIntGraph->Write();
	for(unsigned int i = 0; i < graphs.size(); ++i) {
		graphs[i]->Write();
	}
	outTree->Write();
	integralFile->Close();
	pwd->cd();

}


void integralTableContainer::fillIntegralTable() {

	if(_integralTable.size() > 0) {
		printErr << "trying to fill an already filled integral Table. Aborting..." << endl;
		throw;
	}

	double M = 0.;
	for(unsigned int i = 0; i < _subDecay->nmbFsParticles(); ++i) {
		M += _subDecay->fsParticles()[i]->mass();
	}
	const double step = (UPPER_BOUND - M) / (N_POINTS);
	_integralTable.push_back(integralTablePoint(M));
	M += step;

	const double& M0 = _vertex->parent()->mass();
	printInfo << "filling the integration table assuming the decaying isobar to be "
	          << _vertex->parent()->name() << " with M0=" << M0 << endl;
	_M0s.push_back(M0);

	for(; M <= UPPER_BOUND; M += step) {
		if((M0 > (M-step)) and (M0 < M)) {
			_integralTable.push_back(evalInt(M0, N_MC_EVENTS_FOR_M0));
        }
		_integralTable.push_back(evalInt(M, N_MC_EVENTS));
	}

}


void integralTableContainer::addToIntegralTable(const integralTablePoint& newPoint) {

	vector<integralTablePoint> newIntegralTable(1, _integralTable[0]);
	const double& newM = newPoint.M;
	bool added = false;
	for(unsigned int i = 1; i < _integralTable.size(); ++i) {
		const double& lastM = _integralTable[i-1].M;
		const double& currentM = _integralTable[i].M;
		if((newM >= lastM) and (newM < currentM)) {
			newIntegralTable.push_back(newPoint);
			added = true;
		}
		newIntegralTable.push_back(_integralTable[i]);
	}
	if(not added) {
		printErr << "could not add new point to integral table (out of range?). Aborting..." << endl;
		throw;
	}
	_integralTable = newIntegralTable;

}
