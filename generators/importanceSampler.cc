#include<cmath>
#include<iostream>
#include<sstream>

#include<TFile.h>

#include<BAT/BCMath.h>

#include"importanceSampler.h"
#include"nBodyPhaseSpaceKinematics.h"
#include"physUtils.hpp"
#include"randomNumberGenerator.h"


size_t rpwa::importanceSampler::_nCalls = 0;


rpwa::importanceSampler::importanceSampler(const double            mMin,
                                           const double            mMax,
                                           rpwa::modelIntensityPtr model)
	: BCModel("phaseSpaceImportanceSampling"),
	  _phaseSpaceOnly(false),
	  _productionGeneratorInitialized(false),
	  _nPart(0),
	  _mMin(mMin),
	  _mMax(mMax),
	  _masses(model->finalStateMasses()),
	  _mSum(0.),
	  _model(model),
	  _fileWriter(),
	  _beamAndVertexGenerator(beamAndVertexGeneratorPtr(new beamAndVertexGenerator())),
	  _pickerFunction(massAndTPrimePickerPtr())
{
	_nPart = _model->nFinalState();
	if (_nPart < 2) {
		printErr << "less than two particles to sample. Four-momentum conservation leaves no d.o.f.. No sampling necessary. Aborting..." << std::endl;
		throw;
	}
	if (_mMin > _mMax) {
		printErr << "_mMin = " << _mMin << " > _mMax " << _mMax << ". Aborting..." << std::endl;
		throw;
	}

	// follow the convention of the nBodyPhaseSpace classes
	double mMinPar = _masses[0];
	_mSum          = _masses[0];
	for (size_t part = 1; part < _nPart; ++part) {
		mMinPar += _masses[part];
		_mSum   += _masses[part];

		std::stringstream nameM;
		std::stringstream texM;
		nameM << "m";
		texM << "m_{";
		for (size_t p = 1; p <= part+1; ++p) {
			nameM << p;
			texM << p;
		}
		texM << "}";
		BCParameter param(nameM.str(), mMinPar, _mMax, texM.str());
		if (part == _nPart-1) {
			if (_mMin == _mMax) {
				// X mass is fixed to one value
				param.Fix(_mMax);
			} else {
				param.SetLimits(_mMin, _mMax);
			}
		}
		AddParameter(param);

		// sample cos(theta), because dOmega = sin(theta) dtheta dphi = dcos(theta) dphi
		std::stringstream nameTheta;
		std::stringstream texTheta;
		nameTheta << "cosTheta" << part+1;
		texTheta << "cos(#theta_{" << part+1 << "})";
		AddParameter(nameTheta.str(), -1., 1., texTheta.str());

		std::stringstream namePhi;
		std::stringstream texPhi;
		namePhi << "phi" << part+1;
		texPhi << "#phi_{" << part+1 << "}";
		AddParameter(namePhi.str(), 0., 2*M_PI, texPhi.str());
	}
	for (size_t part = 0; part < _nPart; ++part) {
		for (size_t qart = 0; qart < part; ++qart) {
			std::stringstream name;
			std::stringstream tex;
			name << "M" << qart+1 << part+1 << 2;
			tex << "m^{2}_{" << qart+1 << part+1 << "}";
			AddObservable(name.str(), 0., _mMax*_mMax, tex.str());
		}
	}
	BCAux::SetStyle();
}


double
rpwa::importanceSampler::LogAPrioriProbability(const std::vector<double>& parameters)
{
	rpwa::nBodyPhaseSpaceKinematics nBodyPhaseSpace;
	if (not initializeNBodyPhaseSpace(nBodyPhaseSpace, parameters, false)) {
		return -std::numeric_limits<double>::infinity();
	}

	const double phaseSpace = nBodyPhaseSpace.calcWeight() / std::pow(parameters[3*_nPart-6] - _mSum, _nPart-2);
	return std::log(phaseSpace);
}


double
rpwa::importanceSampler::LogLikelihood(const std::vector<double>& parameters)
{
	++_nCalls;
	if (_phaseSpaceOnly) {
		return -1;
	}

	rpwa::nBodyPhaseSpaceKinematics nBodyPhaseSpace;
	if (not initializeNBodyPhaseSpace(nBodyPhaseSpace, parameters, true)) {
		printErr << "error while initializing n-body phase space. Aborting..." << std::endl;
		throw;
	}

	nBodyPhaseSpace.calcBreakupMomenta();
	nBodyPhaseSpace.calcEventKinematics(TLorentzVector(0., 0., 0., parameters[3*_nPart-6]));
	std::vector<TVector3> decayKinMomenta(_nPart);
	for (size_t part = 0; part < _nPart; ++part) {
		decayKinMomenta[part] = nBodyPhaseSpace.daughter(part).Vect();
	}

	const double intensity = _model->getIntensity(decayKinMomenta);
	return std::log(intensity);
}


void
rpwa::importanceSampler::CalculateObservables(const std::vector<double>& parameters)
{
	if (!_productionGeneratorInitialized) {
		printErr << "Production kinematics not initilized. Cannot generate beam." << std::endl;
		throw;
	}


	rpwa::nBodyPhaseSpaceKinematics nBodyPhaseSpace;
	if (not initializeNBodyPhaseSpace(nBodyPhaseSpace, parameters, true)) {
		printErr << "error while initializing n-body phase space. Aborting..." << std::endl;
		throw;
	}

	nBodyPhaseSpace.calcBreakupMomenta();
	nBodyPhaseSpace.calcEventKinematics(TLorentzVector(0., 0., 0., parameters[3*_nPart-6]));

	size_t count = 0;
	for (size_t part = 0; part < _nPart; ++part) {
		for (size_t qart = 0; qart < part; ++qart) {
			GetObservable(count) = (nBodyPhaseSpace.daughter(part) + nBodyPhaseSpace.daughter(qart)).M2();
			++count;
		}
	}

	if (_fileWriter.initialized() && GetPhase() > 0) { // Put file writer stuff here
		if (fMCMCCurrentIteration % fMCMCNLag != 0) {
			return; // apply lag
		}

		std::pair<std::pair<bool, double>, std::vector<TLorentzVector> > production = getProductionKinematics(parameters[3*_nPart-6]);
		if (!production.first.first) {
			printErr << "generation of production failed" << std::endl;
			throw;
		}

		const TLorentzRotation toLab = rpwa::isobarAmplitude::gjTransform(production.second[0], production.second[1]).Inverse();
		std::vector<TVector3> decayKinMomenta(_nPart);
		TLorentzVector pX(0.,0.,0.,0.);
		for (size_t part = 0; part < _nPart; ++part) {
			TLorentzVector p = nBodyPhaseSpace.daughter(part);
			p.Transform(toLab);
			decayKinMomenta[part] = p.Vect();
			pX += p;
		}

		std::vector<double>var(3+_nPart);
		var[0] = production.first.second;
		var[1] = pX.M();
		var[2] = pX.E();
		for (size_t part = 0; part < _nPart; ++part) {
			var[part + 3] = nBodyPhaseSpace.daughter(part).E();
		}
		std::vector<TVector3> prodKinMomenta(1, production.second[0].Vect());
		_fileWriter.addEvent(prodKinMomenta, decayKinMomenta, var);
	}
}


bool
rpwa::importanceSampler::initializeFileWriter(TFile* outFile)
{
	std::string userString("importanceSampledEvents");
	std::map<std::string, std::pair<double, double> > binning;
	binning["mass"] = std::pair<double, double> (_mMin, _mMax); // Use own mass range, since the range of the picker will be overridden
	if (_productionGeneratorInitialized) {
		binning["tPrime"] = _pickerFunction->tPrimeRange(); // Use tPicker t' range, since the actual sampling is decoupled from t'.
	}
	std::vector<std::string> var(1, "tPrime");
	var.push_back("mass");
	var.push_back("EXlab");
	for (size_t part = 0; part < _nPart; ++part) {
		std::stringstream name;
		name << "E" << part+1;
		var.push_back(name.str());
	}
	return _fileWriter.initialize(*outFile,
	                              userString,
	                              rpwa::eventMetadata::GENERATED,
	                              _model->initialStateParticles(),
	                              _model->finalStateParticles(),
	                              binning,
	                              var);
}


bool
rpwa::importanceSampler::finalizeFileWriter()
{
	return _fileWriter.finalize();
}


bool
rpwa::importanceSampler::initializeProductionGenerator(rpwa::Beam&                     beam,
                                                       rpwa::Target&                   target,
                                                       rpwa::beamAndVertexGeneratorPtr generator,
                                                       rpwa::massAndTPrimePickerPtr    picker)
{
	_productionGeneratorInitialized = false;
	_target = rpwa::Target(target);
	_beam   = rpwa::Beam(beam);
	_beamAndVertexGenerator = rpwa::beamAndVertexGeneratorPtr(generator);
	_pickerFunction = rpwa::massAndTPrimePickerPtr(picker);
	_productionGeneratorInitialized = true;
	return true;
}


std::pair<std::pair<bool, double>, std::vector<TLorentzVector> >
rpwa::importanceSampler::getProductionKinematics(const double xMass)
{
	if (!_productionGeneratorInitialized) {
		printErr << "poduction generators not initalized" << std::endl;
		return std::pair<std::pair<bool, double>, std::vector<TLorentzVector> >(std::pair<bool, double>(false, 0.), std::vector<TLorentzVector>(2));
	}
	if(not _beamAndVertexGenerator->event(_target, _beam)) {
		printErr << "could not generate vertex/beam. Aborting..." << std::endl;
		return std::pair<std::pair<bool, double>, std::vector<TLorentzVector> >(std::pair<bool, double>(false, 0.), std::vector<TLorentzVector>(2));;
	}
	const TLorentzVector targetLab(0., 0., 0., _target.targetParticle.mass());
	_beam.particle.setLzVec(_beamAndVertexGenerator->getBeam());
	const TLorentzVector& beamLorentzVector = _beam.particle.lzVec();
	const TLorentzVector overallCm = beamLorentzVector + targetLab;  // beam-target center-of-mass system


	double tPrime;
	if (not _pickerFunction->pickTPrimeForMass(xMass, tPrime)) {
		printErr << "t' pick failed" << std::endl;
		return std::pair<std::pair<bool, double>, std::vector<TLorentzVector> >(std::pair<bool, double>(false, 0.), std::vector<TLorentzVector>(2));
	}
	TRandom3* random = randomNumberGenerator::instance()->getGenerator();

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
	const double t            = t0 - tPrime;

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

	TVector3 beamDir = _beam.particle.momentum().Unit();
	xSystemLab.RotateUz(beamDir);
	std::vector<TLorentzVector> kinematics(2);
	kinematics[0] = beamLorentzVector;
	kinematics[1] = xSystemLab;
	return std::pair<std::pair<bool, double>, std::vector<TLorentzVector> >(std::pair<bool, double>(true, tPrime), kinematics);
}


bool
rpwa::importanceSampler::initializeNBodyPhaseSpace(rpwa::nBodyPhaseSpaceKinematics& nBodyPhaseSpace,
                                                   const std::vector<double>&       parameters,
                                                   const bool                       angles) const
{
	// currently these options are in the default values, but just to be safe
	nBodyPhaseSpace.setWeightType(rpwa::nBodyPhaseSpaceKinematics::S_U_CHUNG);
	nBodyPhaseSpace.setKinematicsType(rpwa::nBodyPhaseSpaceKinematics::BLOCK);

	if (not nBodyPhaseSpace.setDecay(_masses)) {
		return false;
	}

	std::vector<double> masses(_nPart);
	masses[0] = _masses[0];
	for (size_t part = 1; part < _nPart; ++part) {
		masses[part] = parameters[3*part-3];
		if ((masses[part-1] + _masses[part]) > masses[part])
			return false;
	}
	nBodyPhaseSpace.setMasses(masses);

	if (angles) {
		std::vector<double> cosTheta(_nPart);
		std::vector<double> phi     (_nPart);

		for (size_t part = 1; part < _nPart; ++part) {
			cosTheta[part] = parameters[3*part-2];
			phi     [part] = parameters[3*part-1];
		}

		nBodyPhaseSpace.setAngles(cosTheta, phi);
	}

	return true;
}
