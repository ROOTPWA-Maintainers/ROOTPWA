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


diffractivePhaseSpace::diffractivePhaseSpace()
	: generator()
{
	_phaseSpace.setWeightType    (nBodyPhaseSpaceGen::S_U_CHUNG);
	_phaseSpace.setKinematicsType(nBodyPhaseSpaceGen::BLOCK);
}


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
			printErr << "mass- and t'-picker function has not been set. Aborting..." << endl;
			throw;
		}
		const double xMassMin = _pickerFunction->massRange().first;
		const double xMassMax = _pickerFunction->massRange().second;
		if((xMassMax < xMassMin) || (xMassMax <= 0.)) {
			printErr << "mass range [" << xMassMin << ", " << xMassMax << "] GeV/c^2 "
			         << "does mot make sense. Aborting..." << endl;
			throw;
		} else {
			printInfo << "calculating max weight (" << nmbDaughters << " FS particles) "
			          << "for m = " << xMassMax << " GeV/c^2:";
			_phaseSpace.setMaxWeight(1.01 * _phaseSpace.estimateMaxWeight(xMassMax, 1000000));
			cout << " max weight = " << _phaseSpace.maxWeight() << endl;
		}
	}
}


// based on Dima's prod_decay_split.f
unsigned int
diffractivePhaseSpace::event()
{

	TRandom3* random = randomNumberGenerator::instance()->getGenerator();

	unsigned long int attempts = 0;
  // construct primary vertex and beam
  // use the primary Vertex Generator if available
	if(_beamAndVertexGenerator) {
		_beamAndVertexGenerator->event();
		_vertex = _beamAndVertexGenerator->getVertex();
		_beam.particle.setLzVec(_beamAndVertexGenerator->getBeam());
	} else {
		double x;
		double y;
		double radius = std::sqrt(random->Uniform(0, _target.radius * _target.radius));
		random->Circle(x, y, radius);
		double z;
		do {
			z = random->Uniform();
		} while (random->Uniform() < z*_target.interactionLength);
		z = (_target.position.Z() - _target.length * 0.5) + z * _target.length;
		_vertex.SetXYZ(_target.position.X() + x,
		               _target.position.Y() + y,
		               z);
		_beam.particle.setLzVec(makeBeam());
	}

	if(not _pickerFunction) {
		printErr << "mass- and t'-picker function has not been set. Aborting..." << endl;
		throw;
	}

	const TLorentzVector& beamLorentzVector = _beam.particle.lzVec();

	const double xMassMax = _pickerFunction->massRange().second;
	const TLorentzVector targetLab(0, 0, 0, _target.targetParticle.mass());
	const TLorentzVector overallCm = beamLorentzVector + targetLab;  // beam-target center-of-mass system
	// check
	if(xMassMax + _target.targetParticle.mass()  > overallCm.M()) {
		printErr << "Maximum mass out of kinematic range. "
		         << "Limit = " << overallCm.M() - _target.targetParticle.mass() << " GeV/c2. "
		         << "Aborting..." << endl;
		throw;
	}

	double tPrime;
	double xMass;
	do {
		assert((*_pickerFunction)(xMass, tPrime));
	} while(xMass + _target.recoilParticle.mass() > overallCm.M());
	// t' should be negative (why?)
	tPrime *= -1;

	bool done = false;
	while(!done) {

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


double
diffractivePhaseSpace::calcTPrime(const TLorentzVector& particle_In, const TLorentzVector& particle_Out) {
	double result = 0.;
	result = (particle_Out.M2()-particle_In.M2());
	result = pow(result,2);
	result /= 4*pow(particle_In.P(),2);
	result = fabs((particle_In-particle_Out).M2())-fabs(result);
	return result;
}
