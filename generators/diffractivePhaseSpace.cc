///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2009 Sebastian Neubert
//
//    This file is part of rootpwa
//
//    rootpwa is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    rootpwa is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////


#include <iomanip>
#include <fstream>
#include <limits>

#include "TVector3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TFile.h"

#include "randomNumberGenerator.h"
#include "reportingUtils.hpp"
#include "physUtils.hpp"
#include "diffractivePhaseSpace.h"

using namespace std;
using namespace rpwa;


// /////////////////////////////////////////////////////////////////////////////////
// // global constants


// // target position
// const double _targetZPos = -300.0; // [cm]

// // beam parameters:
// const double _beamMomSigma = 1.2;  // [GeV/c]
// const double _beamMom      = 189;  // [GeV/c]
// // 2004 beam:
// // const double _beamDxDz      = 0.00026; // tilt from Quirin was in mrad
// // const double _beamDxDzSigma = 0.00010;
// // const double _beamDyDz      = 0.00001; // tilt from Quirin was in mrad
// // const double _beamDyDzSigma = 0.00018;
// // ideal beam:
// const double _beamDxDz      = 0.0;
// const double _beamDxDzSigma = 0.0;
// const double _beamDyDz      = 0.0;
// const double _beamDyDzSigma = 0.0;

// // cut on t-distribution
// const double _tMin = 0.001;  // [(GeV/c)^2]


diffractivePhaseSpace::diffractivePhaseSpace()
	: generator(),
//	  _primaryVertexGen(NULL),
	  _tPrime(0.),
//	  _invSlopePar(NULL),
//	  _invM(NULL),
	  _tMin(0.)//,
//	  _tprimeMin(0.),
//	  _tprimeMax(numeric_limits<double>::max()),
//	  _xMassMin(0),
//	  _xMassMax(0),
//	  _protonMass(0.938272013),
//	  _pionMass(0.13957018),
//	  _pionMass2(_pionMass * _pionMass)
{
	_phaseSpace.setWeightType    (nBodyPhaseSpaceGen::S_U_CHUNG);
	_phaseSpace.setKinematicsType(nBodyPhaseSpaceGen::BLOCK);
}

diffractivePhaseSpace::~diffractivePhaseSpace() {
	delete _primaryVertexGen;
}


// /////////////////////////////////////////////////////////////////////////////////
// helper functions

// constructs beam Lorentz vector
TLorentzVector
diffractivePhaseSpace::makeBeam()
{
	TRandom3* random = randomNumberGenerator::instance()->getGenerator();
	// throw magnituide of beam momentum
	const double pBeam = random->Gaus(_beam.momentum, _beam.momentumSigma);
	// throw beam inclination
	const double dxdz = random->Gaus(_beam.DxDz, _beam.DxDzSigma);
	const double dydz = random->Gaus(_beam.DyDz, _beam.DyDzSigma);
	// construct tilted beam momentum Lorentz vector
	const double pz        = pBeam / sqrt(1 + dxdz * dxdz + dydz * dydz);
	const double px        = dxdz * pz;
	const double py        = dydz * pz;
	const double EBeam     = sqrt(pBeam * pBeam + _beam.particle.mass2());
	return TLorentzVector(px, py, pz, EBeam);
}

/*
// writes event to ascii file read by gamp
bool
diffractivePhaseSpace::writePwa2000Ascii(ostream&  out,
                                          const int beamGeantId,
                                          const int beamCharge)
{
	if(!out) {
		printErr << "output stream is not writable." << endl;
		return false;
	}
	unsigned int nmbDaughters = _decayProducts.size();
	out << nmbDaughters + 1 << endl;
	// beam particle: geant ID, charge, p_x, p_y, p_z, E
	out << beamGeantId << " " << beamCharge
	    << setprecision(numeric_limits<double>::digits10 + 1)
	    << " " << _beamLab.Px() << " " << _beamLab.Py() << " " << _beamLab.Pz()
	    << " " << _beamLab.E() << endl;
	for(unsigned int i = 0; i < nmbDaughters; ++i) {
		const TLorentzVector& hadron = _phaseSpace.daughter(i);
		// hadron: geant ID, charge, p_x, p_y, p_z, E
		out << _decayProducts[i].geantId() << " " << _decayProducts[i].charge()
		    << setprecision(numeric_limits<double>::digits10 + 1)
		    << " " << hadron.Px() << " " << hadron.Py() << " " << hadron.Pz()
		    << " " << hadron.E() << endl;
	}
	return true;
}


bool
diffractivePhaseSpace::writeComgeantAscii(ostream& out, bool  binary) {

	if(!out) {
		cerr << "Output stream is not writable." << endl;
		return false;
	}

	if(not binary) { // Write text file.
		// total number of particles including recoil proton and beam particle
		unsigned int nmbDaughters = _decayProducts.size();
		out << nmbDaughters+1+1 << endl;
		// vertex position in cm
		// note that Comgeant's coordinate system is different
		out << _vertex.Z() << " " << _vertex.X() << " " << _vertex.Y() << endl;
		// beam particle: geant ID , -p_z, -p_x, -p_y must go the opposite direction upstream and should be defined as mulike with PID 44 in Comgeant
		out << setprecision(numeric_limits<double>::digits10 + 1)
		    << "44 " << -_beamLab.Pz() << " " << -_beamLab.Px() << " " << -_beamLab.Py() << endl;// << " " << beam.E() << endl;
		// the recoil proton
		out << setprecision(numeric_limits<double>::digits10 + 1)
		    << "14 " << _recoilprotonLab.Pz() << " " << _recoilprotonLab.Px() << " " << _recoilprotonLab.Py() << endl;// << " " << beam.E() << endl;
		for (unsigned int i = 0; i < nmbDaughters; ++i) {
			const TLorentzVector& hadron = _phaseSpace.daughter(i);
			// hadron: geant ID, p_z, p_x, p_y
			out << setprecision(numeric_limits<double>::digits10 + 1)
			    << _decayProducts[i].geantId() << " "
			    << hadron.Pz() << " "
			    << hadron.Px() << " "
			    << hadron.Py() << endl;// << " " << hadron->E() << endl;
		}
	} else {
		int intval;
		float floatval;
		// total number of particles including recoil proton and beam particle
		unsigned int nmbDaughters = _decayProducts.size();
		//out << nmbDaughters+1+1 << endl;
		intval = (int)nmbDaughters+1+1; out.write((char*)&intval,4);
		// vertex position in cm
		// note that Comgeant's coordinate system is different
		floatval = (float)_vertex.Z(); out.write((char*)&floatval,4);
		floatval = (float)_vertex.X(); out.write((char*)&floatval,4);
		floatval = (float)_vertex.Y(); out.write((char*)&floatval,4);
		//out << _vertex.Z() << " " << _vertex.X() << " " << _vertex.Y() << endl;
		// beam particle: geant ID , -p_z, -p_x, -p_y must go the opposite direction upstream and should be defined as mulike with PID 44 in Comgeant
		intval = 44; out.write((char*)&intval,4);
		floatval = (float)-_beamLab.Pz(); out.write((char*)&floatval,4);
		floatval = (float)-_beamLab.Px(); out.write((char*)&floatval,4);
		floatval = (float)-_beamLab.Py(); out.write((char*)&floatval,4);
		//out << setprecision(numeric_limits<double>::digits10 + 1)
		//  << "44 " << -_beamLab.Pz() << " " << -_beamLab.Px() << " " << -_beamLab.Py() << endl;// << " " << beam.E() << endl;
		// the recoil proton
		intval = 14; out.write((char*)&intval,4);
		floatval = (float)_recoilprotonLab.Pz(); out.write((char*)&floatval,4);
		floatval = (float)_recoilprotonLab.Px(); out.write((char*)&floatval,4);
		floatval = (float)_recoilprotonLab.Py(); out.write((char*)&floatval,4);
		//out << setprecision(numeric_limits<double>::digits10 + 1)
		//  << "14 " << _recoilprotonLab.Pz() << " " << _recoilprotonLab.Px() << " " << _recoilprotonLab.Py() << endl;// << " " << beam.E() << endl;
		for (unsigned int i = 0; i < nmbDaughters; ++i) {
			const TLorentzVector& hadron = _phaseSpace.daughter(i);
			// hadron: geant ID, p_z, p_x, p_y
			intval = (int)_decayProducts[i].geantId(); out.write((char*)&intval,4);
			floatval = (float)hadron.Pz(); out.write((char*)&floatval,4);
			floatval = (float)hadron.Px(); out.write((char*)&floatval,4);
			floatval = (float)hadron.Py(); out.write((char*)&floatval,4);
			//out << setprecision(numeric_limits<double>::digits10 + 1)
			//<< _decayProducts[i]._gId << " "
			//<< hadron.Pz() << " "
			//<< hadron.Px() << " "
			//<< hadron.Py() << endl;// << " " << hadron->E() << endl;
		}
	}
	return true;
}
*/
/*
void
diffractivePhaseSpace::setSeed(int seed)
{
	gRandom->SetSeed(seed);
	_phaseSpace.setSeed(seed);
}
*/

void
diffractivePhaseSpace::setDecayProducts(const vector<particle>& particles)
{
	_decayProducts.clear();
	_decayProducts = particles;
	buildDaughterList();
}


void
diffractivePhaseSpace::addDecayProduct(const particle& particle)
{
	_decayProducts.push_back(particle);
	buildDaughterList();
}


void
diffractivePhaseSpace::buildDaughterList()
{
	const unsigned int nmbDaughters = _decayProducts.size();
	vector<double> daughterMasses(nmbDaughters, 0);
	for(unsigned int i = 0; i < nmbDaughters; ++i) {
		daughterMasses[i] = _decayProducts[i].mass();
	}
	if(nmbDaughters > 1) {
		_phaseSpace.setDecay(daughterMasses);
		if(not _pickerFunction) {
			printErr << "mass- and t'-picker function has not been set." <<endl;
			throw;
		}
		const double xMassMin = _pickerFunction->massRange().first;
		const double xMassMax = _pickerFunction->massRange().second;
		if((xMassMax < xMassMin) || (xMassMax <= 0.)) {
			printErr << "mass range [" << xMassMin << ", " << xMassMax << "] GeV/c^2 "
			         << "does mot make sense. exiting." << endl;
			throw;
		} else {
			printInfo << "calculating max weight (" << nmbDaughters << " FS particles) "
			          << "for m = " << xMassMax << " GeV/c^2:";
			_phaseSpace.setMaxWeight(1.01 * _phaseSpace.estimateMaxWeight(xMassMax, 1000000));
			cout << " max weight = " << _phaseSpace.maxWeight() << endl;
		}
	}
}


/////////////////////////////////////////////////////////////////////////////////
// main routine
// void
// diffractivePhaseSpace::genPhaseSpaceData(const double   xMassMin          = 2.100,  // lower bound of mass bin [GeV/c^2]
// 		  const double   xMassMax          = 2.140,  // upper bound of mass bin [GeV/c^2]
// 		  const TString& outFileName       = "2100.2140.genbod.evt",
// 		  const TString& thetaHistFileName = "./hTheta.root",  // histogram with experimental distribution of scattering angle
// 		  const int      nmbEvent          = 2000,
// 		  const bool     plot              = false)
// {
//   Double_t daughterMasses[3] = {_pionMass, _pionMass, _pionMass};
//
//   gRandom->SetSeed(12345);
//
//
//   // setup histograms
//   TH1D* ht;
//   TH1D* hm;
//   TH1D* hTheta;
//   TH3D* hVertex3D;
//   TH2D* hVertex;
//   TH1D* hVz;
//   TH1D* hE;
//   if (plot) {
//     ht        = new TH1D("ht", "t", 10000, -0.1, 1);
//     hm        = new TH1D("hm", "3pi mass", 1000, 0.5, 2.5);
//     hTheta    = new TH1D("hThetaGen", "cos theta", 100, 0.99985, 1);
//     hVertex3D = new TH3D("hVertex3D", "Vertex xyz", 100, -2, 2, 100, -2, 2, 200, _targetZPos - 5, _targetZPos + 40);
//     hVertex   = new TH2D("hVertex", "Vertex xy", 100, -2, 2, 100, -2, 2);
//     hVz       = new TH1D("hVz","Vertex z", 1000, _targetZPos - 40, _targetZPos + 40);
//     hE        = new TH1D("hE", "E", 100, 180, 200);
//   }
//
//   // open output file
//   ofstream outFile(outFileName);
//   cout << "Writing " << nmbEvent << " events to file '" << outFileName << "'." << endl;
//
//   // get theta histogram
//   TH1* thetaDist = NULL;
//   {
//     TFile* thetaHistFile = TFile::Open(thetaHistFileName, "READ");
//     if (!thetaHistFile || thetaHistFile->IsZombie()) {
//       cerr << "Cannot open histogram file '" << thetaHistFileName << "'. exiting." << endl;
//       return;
//     }
//     thetaHistFile->GetObject("h1", thetaDist);
//     if (!thetaDist) {
//       cerr << "Cannot find theta histogram in file '" << thetaHistFileName << "'. exiting." << endl;
//       return;
//     }
//   }
//
//   int countEvent = 0;
//   int attempts   = 0;
//   int tenpercent = (int)(nmbEvent * 0.1);
//   while (countEvent < nmbEvent) { // loop over events
//     ++attempts;
//
//
//   }
// }

/*
double
diffractivePhaseSpace::getInvSlopePar(double invariant_M) {
	double result = 1.;
	if(!_invSlopePar) {
		return result;
	}
	if(_ninvSlopePar == 1) {
		return _invSlopePar[0];
	}
	if(invariant_M < 0.) {
		return _invSlopePar[0];
	}
	// assuming to have sorted entries
	// case of linear extrapolation
	int i_1 = -1;
	int i_2 = -1;
	if(invariant_M < _invM[0]) {
		i_1 = 0;
		i_2 = 1;
	}
	if(invariant_M >= _invM[_ninvSlopePar - 2]) {
		i_1 = _ninvSlopePar - 2;
		i_2 = _ninvSlopePar - 1;
	}
	// case of linear interpolation
	if(i_1 < 0 || i_2 < 0) {
		// search for the matching two points
		for(int i = 0; i < _ninvSlopePar-2; i++) {
			if(_invM[i] <= invariant_M && invariant_M < _invM[i+1]) {
				i_1=i;
				i_2=i+1;
				break;
			}
		}
	}
	// extra/interpolate lineary
	double m_1 = _invM[i_1];
	double m_2 = _invM[i_2];
	double invt_1 = _invSlopePar[i_1];
	double invt_2 = _invSlopePar[i_2];
	result = invt_1 + (((invt_2-invt_1) / (m_2-m_1)) * (invariant_M-m_1));
	return result;
}
*/


// based on Dima's prod_decay_split.f
unsigned int
diffractivePhaseSpace::event()
{

	TRandom3* random = randomNumberGenerator::instance()->getGenerator();

	unsigned long int attempts = 0;
  // construct primary vertex and beam
  // use the primary Vertex Generator if available
	if(_primaryVertexGen) {
		// as documented one may need several attempts to get a vertex position
		// which is valid also for the beam direction measurement
		while(attempts < 1000) {
			_vertex = _primaryVertexGen->getVertex();
			TVector3 beam_dir = _primaryVertexGen->getBeamDir(_vertex);
			if(beam_dir.Mag() == 0) {
				++attempts;
				//cout << " skipping " << endl;
				continue;
			}
			_beam.particle.setLzVec(_primaryVertexGen->getBeamPart(beam_dir));
			break;
		}
		// Just a few events should contain a vertex position with no beam direction information
		if(attempts == 999) {
			cerr << " Error in beam construction. Please check the beam properties loaded correctly!" << endl;
			_beam.particle.setLzVec(makeBeam());
		}
	} else {
		double x;
		double y;
		random->Circle(x, y, _target.radius);
		_vertex.SetXYZ(_target.position.X() + x,
		               _target.position.Y() + y,
		               _target.position.Z() + random->Uniform(-_target.length * 0.5, _target.length * 0.5));
		_beam.particle.setLzVec(makeBeam());
	}

	if(not _pickerFunction) {
		printErr << "mass- and t'-picker function has not been set." <<endl;
		throw;
	}

	const TLorentzVector& beamLorentzVector = _beam.particle.lzVec();

	const double xMassMax = _pickerFunction->massRange().second;
	const TLorentzVector targetLab(0, 0, 0, _target.targetParticle.mass());
	const TLorentzVector overallCm = beamLorentzVector + targetLab;  // beam-target center-of-mass system
	// check
	if(xMassMax + _target.targetParticle.mass()  > overallCm.M()) {
		printErr << "Max Mass out of kinematic range." <<  endl
		         << "Limit = " << overallCm.M()  - _target.targetParticle.mass() << "GeV/c2" << endl
		         << " ABORTING " << flush<< endl;
		throw;
	}


	bool done = false;
	while(!done) {

/*		// _xMassMin == _xMassMax
		double xMass = _xMassMin;
		if(_xMassMin < _xMassMax) {
			xMass = gRandom->Uniform(_xMassMin, _xMassMax); // pick random X mass
		} else if(_xMassMin > _xMassMax) {
			xMass = gRandom->Gaus( _xMassMin, _xMassMax); // pick random X mass arround _xMassMin
		}

		double tPrime = _tprimeMin;
		if(_tprimeMax < _tprimeMin) {
			// calculate the slope parameter depending on the invariant mass
			const double calc_invSlopePar = getInvSlopePar(xMass);
			tPrime = -gRandom->Exp(calc_invSlopePar);  // pick random t'
			//cout << " inv slope par " << _invSlopePar << " gradient " << _invSlopeParGradient << " t' is " << tPrime << endl;
		}
*/
		double tPrime;
		double xMass;
		assert((*_pickerFunction)(xMass, tPrime));

		// make sure that X mass is not larger than maximum allowed mass
		if(xMass + _target.recoilParticle.mass() > overallCm.M()) {
			continue;
		}

		// calculate t from t' in center-of-mass system of collision
		const double s            = overallCm.Mag2();
		const double sqrtS        = sqrt(s);
		const double recoilMass2  = _target.recoilParticle.mass2();
		const double xMass2       = xMass * xMass;
		const double xEnergyCM    = (s - recoilMass2 + xMass2) / (2 * sqrtS);  // breakup energy
		const double xMomCM       = sqrt(xEnergyCM * xEnergyCM - xMass2);      // breakup momentum
		const double beamMass2    = _beam.particle.mass2();
		const double targetMass2  = _target.targetParticle.mass2();
		const double beamEnergyCM = (s - targetMass2 + beamMass2) / (2 * sqrtS);    // breakup energy
		const double beamMomCM    = sqrt(beamEnergyCM * beamEnergyCM - beamMass2);  // breakup momentum
		const double t0           = (xEnergyCM - beamEnergyCM) * (xEnergyCM - beamEnergyCM) -
		                            (xMomCM - beamMomCM) * (xMomCM - beamMomCM);
		const double t            = t0 + tPrime;
		// reject events outside of allowed kinematic region
		if(t > t0) {
			continue;
		}

		// construct X Lorentz-vector in lab frame (= target RF)
		// convention used here: Reggeon = X - beam = target - recoil (momentum transfer from target to beam vertex)
		const double beamEnergy       = beamLorentzVector.E();
		const double beamMomentum     = beamLorentzVector.P();
		const double reggeonEnergyLab = (targetMass2 - recoilMass2 + t) / (2 * _target.targetParticle.mass());  // breakup energy
		const double xEnergyLab       = beamEnergy + reggeonEnergyLab;
		const double xMomLab          = sqrt(xEnergyLab * xEnergyLab - xMass2);
		const double xCosThetaLab     = (t - xMass2 - beamMass2 + 2 * beamEnergy * xEnergyLab) / (2 * beamMomentum * xMomLab);
		const double xSinThetaLab     = sqrt(1 - xCosThetaLab * xCosThetaLab);
		const double xPtLab           = xMomLab * xSinThetaLab;
		const double xPhiLab          = random->Uniform(0., TMath::TwoPi());

		// xSystemLab is defined w.r.t. beam direction
		TLorentzVector xSystemLab = TLorentzVector(xPtLab  * cos(xPhiLab),
		                                           xPtLab  * sin(xPhiLab),
		                                           xMomLab * xCosThetaLab,
		                                           xEnergyLab);
		// rotate according to beam tilt
		TVector3 beamDir = _beam.particle.momentum().Unit();
		xSystemLab.RotateUz(beamDir);
		// calculate the recoil proton properties
		_target.recoilParticle.setLzVec((beamLorentzVector + targetLab) - xSystemLab); // targetLab

		/* check for coplanarity
		   cout << " Energy balance is " << (_beamLab + targetLab).E() << " vs. " << (_recoilprotonLab+xSystemLab).E() << endl;
		   cout << " Momentum balance is " << (_beamLab + targetLab).Mag() << " vs. " << (_recoilprotonLab+xSystemLab).Mag() << endl;
		   cout << " Direction X balance is " << (_beamLab + targetLab).Px() << " vs. " << (_recoilprotonLab+xSystemLab).Px() << endl;
		   cout << " Direction Y balance is " << (_beamLab + targetLab).Py() << " vs. " << (_recoilprotonLab+xSystemLab).Py() << endl;
		   cout << " Direction Z balance is " << (_beamLab + targetLab).Pz() << " vs. " << (_recoilprotonLab+xSystemLab).Pz() << endl;

		// own rotation
		TVector3 vec_direction = _beamLab.Vect().Unit();
		TVector3 vec_origin(0.,0.,1.);
		// get the angle of the vector
		double angle = vec_origin.Angle(vec_direction);
		// get the rotation axis perpendicular to the plane between these both
		TVector3 vec_rotation = vec_origin.Cross(vec_direction);
		vec_rotation = vec_rotation.Unit();
		// rotate around this axis by the given angle
		//particle_null.Rotate  (-angle, vec_rotation);
		_recoilprotonLab.Rotate(-angle, vec_rotation);
		xSystemLab.Rotate(-angle, vec_rotation);

		cout << " delta phi is " << (_recoilprotonLab.Phi()-xSystemLab.Phi())/3.141592654 << endl;
		 */

		// recalculate t' for xcheck or save the generated
		// number directly if you change if (1) to (0) to
		// speed up the process a bit
		if(0) {
			_tPrime = calcTPrime(beamLorentzVector, xSystemLab);
		} else {
			_tPrime = tPrime;
		}

		// apply t cut
		if(t > _tMin /*|| _tPrime > _tprimeMin || _tPrime < _tprimeMax*/) {
			continue;
		}

		// generate n-body phase space for X system
		++attempts;
		{
			_phaseSpace.pickMasses(xMass);

			// correct weight for phase space splitting
			// and for 2-body phase space beam-target
			// (1 / 4pi) * q(sqrt(s), m_x, m_recoil) / sqrt(s)
			const double ps2bodyWMax = breakupMomentum(sqrtS, xMassMax, _target.targetParticle.mass())/sqrtS;
			const double ps2bodyW = breakupMomentum(sqrtS, xMass, _target.targetParticle.mass())/sqrtS;

			const double maxPsWeight = _phaseSpace.maxWeight()  * xMassMax * ps2bodyWMax;
			const double psWeight    = _phaseSpace.calcWeight() * xMass* ps2bodyW;

			if((psWeight / maxPsWeight) < random->Rndm()) {
				continue;
			}
			_phaseSpace.pickAngles();
			_phaseSpace.calcEventKinematics(xSystemLab);
		}

		done = true;
	}
	// event was accepted

	const std::vector<TLorentzVector>& daughters = _phaseSpace.daughters();
	assert(daughters.size() == _decayProducts.size());
	for(unsigned int i = 0; i < daughters.size(); ++i) {
		_decayProducts[i].setLzVec(daughters[i]);
	}

	return attempts;
}

/*
unsigned int
diffractivePhaseSpace::event(ostream& stream)
{
	unsigned int attempts = event();
	//writePwa2000Ascii(stream, 9, -1);  // use pi^- beam
	// use the first particle as the beam particle
	writePwa2000Ascii(stream, _decayProducts[0].geantId(), _decayProducts[0].charge());
	return attempts;
}


unsigned int
diffractivePhaseSpace::event(ostream& stream, ostream& streamComGeant)
{
	unsigned int attempts = event();
	//writePwa2000Ascii(stream, 9, -1);  // use pi^- beam
	// use the first particle as the beam particle
	writePwa2000Ascii(stream, _decayProducts[0].geantId(), _decayProducts[0].charge());
	writeComgeantAscii(streamComGeant, false);
	return attempts;
}
*/
/*
void
diffractivePhaseSpace::setBeam(const Beam& beam)
{
	_beam = beam;
}
*/

float
diffractivePhaseSpace::calcTPrime(const TLorentzVector& particle_In, const TLorentzVector& particle_Out){
	float result = 0.;
	result = (particle_Out.M2()-particle_In.M2());
	result = pow(result,2);
	result /= 4*pow(particle_In.P(),2);
	result = fabs((particle_In-particle_Out).M2())-fabs(result);
	return result;
}

