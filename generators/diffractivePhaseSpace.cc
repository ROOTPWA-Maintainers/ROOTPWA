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
	: generator(),
	  _phaseSpace(),
	  _maxXMassSlices(),
	  _maxWeightsForXMasses()
{
	_phaseSpace.setWeightType    (nBodyPhaseSpaceGen::S_U_CHUNG);
	_phaseSpace.setKinematicsType(nBodyPhaseSpaceGen::BLOCK);
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
			unsigned int numberOfMassSlices = _numberOfMassSlices;
			if((xMassMax - xMassMin) < 0.5) {
				numberOfMassSlices = 1;
			}
			double massRangeSliceWidth = (xMassMax - xMassMin) / numberOfMassSlices;
			for(unsigned int i = 1; i <= numberOfMassSlices; ++i) {
				_maxXMassSlices.push_back(xMassMin + (massRangeSliceWidth * i));
			}
			for(unsigned int i = 0; i < _maxXMassSlices.size(); ++i) {
				printInfo << "calculating max weight (" << nmbDaughters << " FS particles) "
				          << "for m = " << _maxXMassSlices[i] << " GeV/c^2:";
				_phaseSpace.setMaxWeight(1.01 * _phaseSpace.estimateMaxWeight(_maxXMassSlices[i], 1000000));
				_maxWeightsForXMasses.push_back(_phaseSpace.maxWeight());
				cout << " max weight = " << _phaseSpace.maxWeight() << endl;
			}
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
	if(not _beamAndVertexGenerator) {
		printErr << "beam and vertex generator has not been set. Aborting..." << endl;
		throw;
	}
	if(not _beamAndVertexGenerator->event(_target, _beam)) {
		printErr << "could not generate vertex/beam. Aborting..." << endl;
		throw;
	}

	_vertex = _beamAndVertexGenerator->getVertex();
	_beam.particle.setLzVec(_beamAndVertexGenerator->getBeam());

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

	bool done = false;
	do {
		do {
			if(not (*_pickerFunction)(_xMass, _tPrime)) {
				printErr << "could not generate X mass and t'. Aborting..." << endl;
				throw;
			}
		} while(_xMass + _target.recoilParticle.mass() > overallCm.M());

		{
			unsigned int i = 0;
			for(; _xMass > _maxXMassSlices[i]; ++i);
			_phaseSpace.setMaxWeight(_maxWeightsForXMasses[i]);
		}

		// calculate t from t' in center-of-mass system of collision
		const double s            = overallCm.Mag2();
		const double sqrtS        = sqrt(s);
		const double recoilMass2  = _target.recoilParticle.mass2();
		const double xMass2       = _xMass * _xMass;
		const double xEnergyCM    = (s - recoilMass2 + xMass2) / (2 * sqrtS);  // breakup energy
		const double xMomCM       = sqrt(xEnergyCM * xEnergyCM - xMass2);      // breakup momentum
		const double beamMass2    = _beam.particle.mass2();
		const double targetMass2  = _target.targetParticle.mass2();
		const double beamEnergyCM = (s - targetMass2 + beamMass2) / (2 * sqrtS);    // breakup energy
		const double beamMomCM    = sqrt(beamEnergyCM * beamEnergyCM - beamMass2);  // breakup momentum
		const double t0           = (xEnergyCM - beamEnergyCM) * (xEnergyCM - beamEnergyCM) -
		                            (xMomCM - beamMomCM) * (xMomCM - beamMomCM); // t0 <= 0
		const double t            = t0 - _tPrime;
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

		do {
			// generate n-body phase space for X system
			++attempts;

			_phaseSpace.pickMasses(_xMass);

			// correct weight for phase space splitting
			// and for 2-body phase space beam-target
			// (1 / 4pi) * q(sqrt(s), m_x, m_recoil) / sqrt(s)
			const double ps2bodyWMax = breakupMomentum(sqrtS, xMassMax, _target.targetParticle.mass())/sqrtS;
			const double ps2bodyW = breakupMomentum(sqrtS, _xMass, _target.targetParticle.mass())/sqrtS;

			const double maxPsWeight = _phaseSpace.maxWeight()  * xMassMax * ps2bodyWMax;
			const double psWeight    = _phaseSpace.calcWeight() * _xMass   * ps2bodyW;

			if((psWeight / maxPsWeight) < random->Rndm()) {
				continue;
			}
			_phaseSpace.pickAngles();
			_phaseSpace.calcEventKinematics(xSystemLab);

			done = true;
		} while(!done);
	} while(!done);
	// event was accepted

	const std::vector<TLorentzVector>& daughters = _phaseSpace.daughters();
	if(daughters.size() != _decayProducts.size()) {
		printErr << "size mismatch between daughters and _decayProducts. Aborting..." << endl;
		throw;
	}
	for(unsigned int i = 0; i < daughters.size(); ++i) {
		_decayProducts[i].setLzVec(daughters[i]);
	}

	return attempts;
}
