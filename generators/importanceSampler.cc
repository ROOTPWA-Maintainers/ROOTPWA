#include<cmath>
#include<iostream>
#include<sstream>

#include<TFile.h>
#include<TStopwatch.h>

#include<BAT/BCMath.h>

#include"importanceSampler.h"
#include"nBodyPhaseSpaceKinematics.h"
#include"randomNumberGenerator.h"


rpwa::importanceSampler::importanceSampler(rpwa::modelIntensityPtr         model,
                                           rpwa::beamAndVertexGeneratorPtr beamAndVertexGenerator,
                                           rpwa::massAndTPrimePickerPtr    massAndTPrimePicker,
                                           const rpwa::Beam&               beam,
                                           const rpwa::Target&             target,
                                           const rpwa::FinalState&         finalState)
	: BCModel("phaseSpaceImportanceSampling"),
	  _model(model),
	  _beamAndVertexGenerator(beamAndVertexGenerator),
	  _massAndTPrimePicker(massAndTPrimePicker),
	  _beam(beam),
	  _target(target),
	  _finalState(finalState),
	  _nPart(_finalState.particles.size()),
	  _phaseSpaceOnly(false)
{
	std::vector<std::string> decayKinParticleNames(_nPart);
	_masses.resize(_nPart);
	_mSum = 0;
	for (size_t part = 0; part < _nPart; ++part) {
		decayKinParticleNames[part] = _finalState.particles[part].name();
		_masses              [part] = _finalState.particles[part].mass();
		_mSum                      += _finalState.particles[part].mass();
	}

	_model->initAmplitudes(decayKinParticleNames);

	if (_nPart < 2) {
		printErr << "less than two particles to sample. Four-momentum conservation leaves no d.o.f.. No sampling necessary. Aborting..." << std::endl;
		throw;
	}
	const double mMin = _massAndTPrimePicker->massRange().first;
	const double mMax = _massAndTPrimePicker->massRange().second;
	if (mMin > mMax) {
		printErr << "mMin = " << mMin << " > mMax = " << mMax << ". Aborting..." << std::endl;
		throw;
	}

	// follow the convention of the nBodyPhaseSpace classes
	double mMinPar = _masses[0];
	for (size_t part = 1; part < _nPart; ++part) {
		mMinPar += _masses[part];

		std::stringstream nameM;
		std::stringstream texM;
		nameM << "m";
		texM << "m_{";
		for (size_t p = 1; p <= part+1; ++p) {
			nameM << p;
			texM << p;
		}
		texM << "}";
		BCParameter param(nameM.str(), mMinPar, mMax, texM.str());
		if (part == _nPart-1) {
			if (mMin == mMax) {
				// X mass is fixed to one value
				param.Fix(mMax);
			} else {
				param.SetLimits(mMin, mMax);
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
			AddObservable(name.str(), 0., mMax*mMax, tex.str());
		}
	}

	resetFuncInfo();
	BCAux::SetStyle();
}


double
rpwa::importanceSampler::LogAPrioriProbability(const std::vector<double>& parameters)
{
	++(_funcCallInfo[LOGAPRIORIPROBABILITY].nmbCalls);
	TStopwatch timerTot;
	timerTot.Start();

	rpwa::nBodyPhaseSpaceKinematics nBodyPhaseSpace;
	if (not initializeNBodyPhaseSpace(nBodyPhaseSpace, parameters, false)) {
		return -std::numeric_limits<double>::infinity();
	}

	const double phaseSpace = nBodyPhaseSpace.calcWeight() / std::pow(parameters[3*_nPart-6] - _mSum, _nPart-2);

	timerTot.Stop();
	_funcCallInfo[LOGAPRIORIPROBABILITY].totalTime += timerTot.RealTime();

	return std::log(phaseSpace);
}


double
rpwa::importanceSampler::LogLikelihood(const std::vector<double>& parameters)
{
	++(_funcCallInfo[LOGLIKELIHOOD].nmbCalls);
	TStopwatch timerTot;
	timerTot.Start();

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

	timerTot.Stop();
	_funcCallInfo[LOGLIKELIHOOD].totalTime += timerTot.RealTime();

	return std::log(intensity);
}


void
rpwa::importanceSampler::CalculateObservables(const std::vector<double>& parameters)
{
	++(_funcCallInfo[CALCULATEOBSERVABLES].nmbCalls);
	TStopwatch timerTot;
	timerTot.Start();

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

		bool           valid;
		double         tPrime;
		TLorentzVector pBeam;
		TLorentzVector pX;
		boost::tuples::tie(valid, tPrime, pBeam, pX) = getProductionKinematics(parameters[3*_nPart-6]);
		if (not valid) {
			printErr << "generation of production failed" << std::endl;
			throw;
		}

		const TLorentzRotation toLab = rpwa::isobarAmplitude::gjTransform(pBeam, pX).Inverse();
		std::vector<TVector3> decayKinMomenta(_nPart);
		for (size_t part = 0; part < _nPart; ++part) {
			TLorentzVector p = nBodyPhaseSpace.daughter(part);
			p.Transform(toLab);
			decayKinMomenta[part] = p.Vect();
		}

		std::vector<double> var;
		if (_storeMassAndTPrime) {
			var.push_back(pX.M());
			var.push_back(tPrime);
		}
		std::vector<TVector3> prodKinMomenta(1, pBeam.Vect());
		_fileWriter.addEvent(prodKinMomenta, decayKinMomenta, var);
	}

	timerTot.Stop();
	_funcCallInfo[CALCULATEOBSERVABLES].totalTime += timerTot.RealTime();
}


bool
rpwa::importanceSampler::initializeFileWriter(TFile*             outFile,
                                              const std::string& userString,
                                              const bool         storeMassAndTPrime,
                                              const std::string& massVariableName,
                                              const std::string& tPrimeVariableName)
{
	_storeMassAndTPrime = storeMassAndTPrime;

	std::vector<std::string> prodKinParticleNames(1);
	prodKinParticleNames[0] = _beam.particle.name();
	std::vector<std::string> decayKinParticleNames(_finalState.particles.size());
	for (size_t part = 0; part < _nPart; ++part) {
		decayKinParticleNames[part] = _finalState.particles[part].name();
	}

	std::map<std::string, std::pair<double, double> > binning;
	binning[massVariableName]   = _massAndTPrimePicker->massRange();
	binning[tPrimeVariableName] = _massAndTPrimePicker->tPrimeRange();
	std::vector<std::string> var;
	if (_storeMassAndTPrime) {
		var.push_back(massVariableName);
		var.push_back(tPrimeVariableName);
	}

	const bool valid = _fileWriter.initialize(*outFile,
	                                          userString,
	                                          rpwa::eventMetadata::REAL,
	                                          prodKinParticleNames,
	                                          decayKinParticleNames,
	                                          binning,
	                                          var);

	return valid;
}


bool
rpwa::importanceSampler::finalizeFileWriter()
{
	return _fileWriter.finalize();
}


boost::tuples::tuple<bool, double, TLorentzVector, TLorentzVector>
rpwa::importanceSampler::getProductionKinematics(const double xMass) const
{
	if(not _beamAndVertexGenerator->event(_target, _beam)) {
		printErr << "could not generate vertex/beam. Aborting..." << std::endl;
		return boost::tuples::make_tuple(false, 0., TLorentzVector(), TLorentzVector());
	}
	const TLorentzVector targetLab(0., 0., 0., _target.targetParticle.mass());
	const TLorentzVector beamLorentzVector = _beamAndVertexGenerator->getBeam();
	const TLorentzVector overallCm = beamLorentzVector + targetLab;  // beam-target center-of-mass system


	double tPrime;
	if (not _massAndTPrimePicker->pickTPrimeForMass(xMass, tPrime)) {
		printErr << "t' pick failed" << std::endl;
		return boost::tuples::make_tuple(false, 0., TLorentzVector(), TLorentzVector());
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

	TVector3 beamDir = beamLorentzVector.Vect().Unit();
	xSystemLab.RotateUz(beamDir);
	return boost::tuples::make_tuple(true, tPrime, beamLorentzVector, xSystemLab);
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


void
rpwa::importanceSampler::resetFuncInfo()
{
	for (unsigned int i = 0; i < NMB_FUNCTIONCALLENUM; ++i) {
		_funcCallInfo[i].nmbCalls  = 0;
		_funcCallInfo[i].totalTime = 0;
	}
}


std::ostream&
rpwa::importanceSampler::printFuncInfo(std::ostream& out) const
{
	const std::string funcNames[NMB_FUNCTIONCALLENUM] = {"LogAPrioriProbability", "LogLikelihood", "CalculateObservables"};
	for (unsigned int i = 0; i < NMB_FUNCTIONCALLENUM; ++i)
		if (_funcCallInfo[i].nmbCalls > 0)
			out << "importanceSampler::" << funcNames[i] << "():" << std::endl
			    << "    number of calls ... " << _funcCallInfo[i].nmbCalls << std::endl
			    << "    total time ........ " << _funcCallInfo[i].totalTime << " sec" << std::endl;
	return out;
}
