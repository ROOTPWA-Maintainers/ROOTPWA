
#include "phaseSpaceIntegral.h"

#include<TDirectory.h>
#include<TF1.h>
#include<TFile.h>
#include<TObjString.h>
#include<TTree.h>

#include"physUtils.hpp"

using namespace std;
using namespace rpwa;

phaseSpaceIntegral* phaseSpaceIntegral::_instance = 0;
const std::string phaseSpaceIntegral::TREE_NAME = "psint";
const std::string phaseSpaceIntegral::DIRECTORY = "/home/kbicker/analysis/integralAmplitudesPwd";

phaseSpaceIntegral* phaseSpaceIntegral::instance() {

	if(not _instance) {
		_instance = new phaseSpaceIntegral();
	}

	return _instance;

}

complex<double> phaseSpaceIntegral::operator()(const isobarDecayVertex& vertex) {

	_vertex = &vertex;

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

	std::vector<std::string> filenames = getFilenames();
	for(unsigned int i = 0; i < filenames.size(); ++i) {
		integralFile = TFile::Open(filenames[i].c_str(), "READ");
		if(integralFile) {
			break;
		}
	}

	if(not integralFile) {
		createIntegralFile(filenames[0]);
		integralFile = TFile::Open(filenames[0].c_str(), "READ");
	}
	TTree* tree = (TTree*)integralFile->Get(TREE_NAME.c_str());
	double psInt = readIntegralValueFromTree(_vertex->parent()->lzVec().M(), tree);
	double psInt0 = readIntegralValueFromTree(_vertex->parent()->mass(), tree);
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
	return (lastPsInt + (((psInt - lastPsInt) / (currentM - lastM)) * (M - currentM)));

}


std::vector<std::string> phaseSpaceIntegral::getFilenames() const {

	std::vector<std::string> retval;

	std::string parentName = _vertex->parent()->name();
	std::string daughter1Name = _vertex->daughter1()->name();
	std::string daughter2Name = _vertex->daughter2()->name();

	std::stringstream sstr;

	sstr<<DIRECTORY<<"/"<<parentName<<"---"<<daughter1Name<<"-"<<daughter2Name<<".root";
	retval.push_back(sstr.str());
	sstr.str("");

	sstr<<DIRECTORY<<"/"<<parentName<<"---"<<daughter2Name<<"-"<<daughter1Name<<".root";
	retval.push_back(sstr.str());
	sstr.str("");

	return retval;

}

double phaseSpaceIntegral::phaseSpace1D(double* x,
                                        double* p) const
{
	const double twoBodyMass   = x[0];
	const double threeBodyMass = p[0];

	// breakup momentum between daughter particles of isobar

	boost::shared_ptr<particle> isobarParticle = _vertex->daughter1();
	boost::shared_ptr<const particle> stableParticle = _vertex->daughter2();
	if(isobarParticle->isStable()) {
		isobarParticle = _vertex->daughter2();
		stableParticle = _vertex->daughter1();
	}
/*
	std::cout<<"threeBodyMass="<<threeBodyMass<<std::endl;
	std::cout<<"twoBodyMass="<<twoBodyMass<<std::endl;
	std::cout<<"stableParticle->mass()="<<stableParticle->mass()<<std::endl;
*/
	// breakup momentum between isobar and bachelor particle
	const double q   = breakupMomentum(threeBodyMass, stableParticle->mass(), twoBodyMass);

	isobarDecayVertexPtr isobarVertex = (boost::dynamic_pointer_cast<isobarDecayVertex>(_decay->toVertex(isobarParticle)))->clone(true, true);
	TLorentzVector lzv;
	TVector3 vec(isobarVertex->parent()->momentum());
	lzv.SetXYZM(vec.X(), vec.Y(), vec.Z(), twoBodyMass);
	isobarVertex->parent()->setLzVec(lzv);

/*	std::cout<<"isobarVertex->parent()->lzVec().M()="<<isobarVertex->parent()->lzVec().M()<<std::endl;

	std::cout<<std::endl;
	std::cout<<"isobarVertex->daughter1()->lzVec().M()="<<isobarVertex->daughter1()->lzVec().M()<<std::endl;
	std::cout<<"isobarVertex->daughter2()->lzVec().M()="<<isobarVertex->daughter2()->lzVec().M()<<std::endl;
	std::cout<<std::endl;
*/
	const double q12 = rpwa::breakupMomentum(twoBodyMass, isobarVertex->daughter1()->lzVec().M(), isobarVertex->daughter2()->lzVec().M());

//	std::cout<<"before massDepAmplitude()"<<std::endl;

	const std::complex<double> bwAmp = isobarVertex->massDepAmplitude();

//	std::cout<<"after massDepAmplitude()"<<std::endl;

	return q * rpwa::barrierFactorSquared(_vertex->L(), q) * q12 * std::norm(bwAmp);
}

double phaseSpaceIntegral::evalInt(const double& M) const {
	TF1 f("f", this, &phaseSpaceIntegral::phaseSpace1D, 0, 1, 1, "threeBodyDynAmpInt", "phaseSpace");
	f.SetParameter(0, M);
	const double pionMass = 0.13957;
	return f.Integral(pionMass * 2., M - pionMass);
}


void phaseSpaceIntegral::createIntegralFile(std::string filename) const {
	TDirectory* pwd = gDirectory;
	TFile* integralFile = TFile::Open(filename.c_str(), "NEW");
	integralFile->cd();
	TObjString filenameForRoot(filename.c_str());
	filenameForRoot.Write();
	TTree* outTree = new TTree(TREE_NAME.c_str(), TREE_NAME.c_str());
	double M, psInt;
	outTree->Branch("M", &M);
	outTree->Branch("int", &psInt);

	const double step = (UPPER_BOUND - LOWER_BOUND) / (N_POINTS - 1);

	M = LOWER_BOUND;
	for(unsigned int i = 0; i < N_POINTS; ++i) {
		psInt = evalInt(M);
		outTree->Fill();
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
