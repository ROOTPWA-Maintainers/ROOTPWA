
#include "phaseSpaceIntegral.h"

#include <boost/progress.hpp>

#include<TClonesArray.h>
#include<TDirectory.h>
#include<TF1.h>
#include<TFile.h>
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
const std::string phaseSpaceIntegral::TREE_NAME = "psint";
const std::string phaseSpaceIntegral::DIRECTORY = "/home/kbicker/analysis/integralAmplitudesPwd";

phaseSpaceIntegral* phaseSpaceIntegral::instance() {

	if(not _instance) {
		_instance = new phaseSpaceIntegral();
		randomNumberGenerator::instance()->setSeed(MC_SEED);
	}

	return _instance;

}

complex<double> phaseSpaceIntegral::operator()(const isobarDecayVertex& vertex) {

	_vertex = isobarDecayVertexPtr();
	const isobarDecayTopologyPtr mainTopo = boost::dynamic_pointer_cast<isobarDecayTopology>(vertex.decay());
	if(not mainTopo) {
		printErr << "got NULL pointer from vertex. Aborting..." << std::endl;
		throw;
	}
	const std::vector<interactionVertexPtr>& decayVertices = mainTopo->decayVertices();
	bool success = false;
	for(unsigned int i = 0; i < decayVertices.size(); ++i) {
		if(*(decayVertices[i]) == vertex) {
			_vertex = boost::dynamic_pointer_cast<isobarDecayVertex>(decayVertices[i]);
			success = true;
			break;
		}
	}
	if(not success) {
		printErr << "decay vertex of particle " << vertex.parent()->label() << " not found. Aborting..." << endl;
	}
	if(not _vertex) {
		printErr << "coud not find decay vertex in decay topology. Aborting..." << endl;
		throw;
	}

	_subDecay = createIsobarDecayTopology(mainTopo->subDecayConsistent(_vertex));
	std::stringstream sstr;
	sstr<<DIRECTORY<<"/"<<waveDescription::waveNameFromTopology(*_subDecay, NEW_FILENAME_CONVENTION)<<".root";
	_filename = sstr.str();

	const particlePtr& parent = _vertex->parent();

	// get Breit-Wigner parameters
	const double       M      = parent->lzVec().M();         // parent mass
	const double       M0     = parent->mass();              // resonance peak position
	const double       Gamma0 = parent->width();             // resonance peak width

	const double Gamma = Gamma0 * dyn();
	// A / (B - iC) = (A / (B^2 + C^2)) * (B + iC)
	const double A = M0 * Gamma0;
	const double B = M0 * M0 - M * M;
	const double C = M0 * Gamma;
	return (A / (B * B + C * C)) * std::complex<double>(B, C);
	// return (M0 * Gamma0) / (M0 * M0 - M * M - imag * M0 * Gamma);
}

double phaseSpaceIntegral::dyn() {

	TFile* integralFile = 0;

	integralFile = TFile::Open(_filename.c_str(), "READ");

	if(not integralFile) {
		printInfo << "no file '" << _filename << "' with integral table found. Creating it..." << endl;
		createIntegralFile();
		printInfo << "created integral file '" << _filename << "." << endl;
		integralFile = TFile::Open(_filename.c_str(), "READ");
	}
	TTree* tree = (TTree*)integralFile->Get(TREE_NAME.c_str());
	double psInt = readIntegralValueFromTree(_vertex->parent()->lzVec().M(), tree);
	double psInt0 = readIntegralValueFromTree(_vertex->parent()->mass(), tree);
	integralFile->Close();
	return (psInt / psInt0);
}


double phaseSpaceIntegral::readIntegralValueFromTree(const double& M, TTree* tree) const
{

	double currentM, psInt;
	tree->SetBranchAddress("M", &currentM);
	tree->SetBranchAddress("int", &psInt);
	tree->GetEntry(0);
	double lastM = currentM;
	double lastPsInt = psInt;
	for(long i = 1; i < tree->GetEntries(); ++i) {
		tree->GetEntry(i);
		if((lastM <= M) and (M < currentM)) {
			break;
		}
	}
	return (lastPsInt + (((psInt - lastPsInt) / (currentM - lastM)) * (M - lastM)));

}


double phaseSpaceIntegral::evalInt(const double& M, const unsigned int& nEvents) const {

	printInfo << "calculating integral for " << _filename << " at mother mass = " << M << "GeV." << endl;

	const massDependencePtr originalMassDep = _vertex->massDependence();
	const decayTopologyPtr originalDecayTopology = _vertex->decay();
	_vertex->setMassDependence(createFlatMassDependence());

	const unsigned int nmbFsParticles = _subDecay->nmbFsParticles();

	isobarHelicityAmplitude subAmp = isobarHelicityAmplitude(_subDecay);
	TClonesArray prodNames("TObjString", 1);
	TClonesArray decayNames("TObjString", nmbFsParticles);
	new (prodNames[0]) TObjString(_vertex->parent()->name().c_str());
	std::vector<double> daughterMasses(nmbFsParticles, 0.);
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
		integral += ampVal * psGen.eventWeight();
		++progressIndicator;
	}

	_vertex->setMassDependence(originalMassDep);
	_vertex->setDecay(originalDecayTopology);
	const double V = fourPi * (1 / rpwa::factorial<double>(nmbFsParticles)) * pow((fourPi*(M - fsParticlesMassSum)), nmbFsParticles-2);
	integral *= (V / (double)nEvents);
	printSucc << "calculated integral: " << integral << std::endl;
	return integral;
}


void phaseSpaceIntegral::createIntegralFile() const {
	TDirectory* pwd = gDirectory;
	TFile* integralFile = TFile::Open(_filename.c_str(), "NEW");
	integralFile->cd();
	TObjString filenameForRoot(_filename.c_str());
	filenameForRoot.Write();
	TTree* outTree = new TTree(TREE_NAME.c_str(), TREE_NAME.c_str());
	double M, psInt;
	outTree->Branch("M", &M);
	outTree->Branch("int", &psInt);

	const double step = (UPPER_BOUND - LOWER_BOUND) / (N_POINTS - 1);

	const double M0 = _vertex->parent()->mass();

	M = LOWER_BOUND;
	for(unsigned int i = 0; i < N_POINTS; ++i) {
		psInt = evalInt(M, N_MC_EVENTS);
		outTree->Fill();
		if((M0 > M) and (M0 < (M+step))) {
			const double oldM = M;
			M = M0;
			psInt = evalInt(M, N_MC_EVENTS_FOR_M0);
			outTree->Fill();
			M = oldM;
		}
		M += step;
	}
/*
	for(double M = LOWER_BOUND; M <= UPPER_BOUND; M += step) {
		psInt = evalInt(M);
		outTree->Fill();
//		std::cout<<M<<std::endl;
	}
*/
	outTree->Write();
	integralFile->Close();
	pwd->cd();

}
