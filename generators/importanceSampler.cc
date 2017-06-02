//-------------------------------------------------------------------------
//
// This class creates events according to a model by sampling the phase
// space using the Bayesian Analysis Toolkit (BAT). To this end this class
// inherits from the BAT model class and overrides the
// 'LogAPrioriProbability', 'LogLikelihood', and 'CalculateObservables'
// class functions:
// * 'LogAPrioriProbability' ensures that the events are inside the
//   kinematic boundaries
// * 'LogLikelihood' calculates the weight of an event given the model by
//   calculating the decay amplitudes
// * 'CalculateObservables' is used to write the events to an event file
//   after boosting it to the laboratory system
//
// For a n-particle final state 3*(n-1) parameters are defined in BAT.
// Each of the (n-1) triples consist of one mass and two angles
// (cos(theta) and phi). The ordering of the triples follows the ordering
// of the parameters in the 'nBodyPhaseSpace...' classes and is sketched
// below:
//
// n-body                              3-body                      2-body                      single daughter
//
// _masses[n - 1]                      _masses[2]                  _masses[1]
// ^                                   ^                           ^
// |                                   |                           |
// |                                   |                           |
// M         [n - 1] = p[3*n - 6]  --> M         [2] = p[3]    --> M         [1] = p[0]    --> M         [0] = _masses[0]
// cos(theta)[n - 1] = p[3*n - 5]      cos(theta)[2] = p[4]        cos(theta)[1] = p[1]        cos(theta)[0] (not used)
// phi       [n - 1] = p[3*n - 4]      phi       [2] = p[5]        phi       [1] = p[2]        phi       [0] (not used)
// _mSum[n - 1]                        _mSum[2]                    _mSum[1]                    _mSum[0] = _masses[0]
// = _mSum[n - 2] + _masses[n - 1]     = _mSum[1] + _masses[2]     = _mSum[0] + _masses[1]
//
//
//-------------------------------------------------------------------------


#include<cmath>
#include<iostream>
#include<sstream>

#include<TFile.h>
#include<TStopwatch.h>

#include<BAT/BCMath.h>

#include"diffractivePhaseSpace.h"
#include"importanceSampler.h"
#include"multibinTypes.h"
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
	  _phaseSpaceOnly(false),
	  _massPrior(0)
{
	std::vector<std::string> decayKinParticleNames(_nPart);
	_masses.resize(_nPart);
	_mSum = 0;
	for (size_t part = 0; part < _nPart; ++part) {
		decayKinParticleNames[part] = _finalState.particles[part].name();
		_masses              [part] = _finalState.particles[part].mass();
		_mSum                      += _finalState.particles[part].mass();
	}

	_model->initDecayAmplitudes(decayKinParticleNames);

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

		// define parameter corresponding to (part+1)-body mass
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

		// the parameter for the X mass is special, its range is taken
		// from the generator options
		if (part == _nPart-1) {
			if (mMin == mMax) {
				// X mass is fixed to one value
				param.Fix(mMax);
			} else {
				param.SetLimits(mMin, mMax);
			}
		}

		AddParameter(param);

		// define parameters corresponding to angles of the bachelor
		// particle in (part+1)-body rest frame
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

	// also show the squared masses of all possible two-body sub-systems
	// and their correlations in the plots created by BAT, in particular
	// for three-body final-states this creates Dalitz plots
	for (size_t partA = 0; partA < _nPart; ++partA) {
		for (size_t partB = 0; partB < partA; ++partB) {
			std::stringstream name;
			std::stringstream tex;
			name << "M" << partB+1 << partA+1 << 2;
			tex << "m^{2}_{" << partB+1 << partA+1 << "}";
			AddObservable(name.str(), 0., mMax*mMax, tex.str());
		}
	}

	resetFuncInfo();
	BCAux::SetStyle();
}


double
rpwa::importanceSampler::LogAPrioriProbability(const std::vector<double>& parameters)
{
#ifdef _OPENMP
#pragma omp atomic
#endif
	++(_funcCallInfo[LOGAPRIORIPROBABILITY].nmbCalls);
	TStopwatch timerTot;
	timerTot.Start();

	rpwa::nBodyPhaseSpaceKinematics nBodyPhaseSpace;
	if (not initializeNBodyPhaseSpace(nBodyPhaseSpace, parameters, false)) {
		return -std::numeric_limits<double>::infinity();
	}

	const double mX = parameters[3*_nPart-6];
	double phaseSpace = nBodyPhaseSpace.calcWeight() / std::pow(mX - _mSum, _nPart-2);
	if (_massPrior) {
#ifdef _OPENMP
#pragma omp critical(massPrior)
#endif
		{
			phaseSpace *= _massPrior->Eval(mX);
		}
	}

	timerTot.Stop();
#ifdef _OPENMP
#pragma omp atomic
#endif
	_funcCallInfo[LOGAPRIORIPROBABILITY].totalTime += timerTot.RealTime();

	return std::log(phaseSpace);
}


double
rpwa::importanceSampler::LogLikelihood(const std::vector<double>& parameters)
{
#ifdef _OPENMP
#pragma omp atomic
#endif
	++(_funcCallInfo[LOGLIKELIHOOD].nmbCalls);
	TStopwatch timerTot;
	timerTot.Start();

	if (_phaseSpaceOnly) {
		return 0;
	}

	rpwa::nBodyPhaseSpaceKinematics nBodyPhaseSpace;
	if (not initializeNBodyPhaseSpace(nBodyPhaseSpace, parameters, true)) {
		printErr << "error while initializing n-body phase space. Aborting..." << std::endl;
		throw;
	}

	const double mX = parameters[3*_nPart-6];
	nBodyPhaseSpace.calcBreakupMomenta();
	nBodyPhaseSpace.calcEventKinematics(TLorentzVector(0., 0., 0., mX));
	std::vector<TVector3> decayKinMomenta(_nPart);
	for (size_t part = 0; part < _nPart; ++part) {
		decayKinMomenta[part] = nBodyPhaseSpace.daughter(part).Vect();
	}

	double intensity;
#ifdef _OPENMP
#pragma omp critical(model)
#endif
	{
		intensity = _model->getIntensity(decayKinMomenta);
	}

	timerTot.Stop();
#ifdef _OPENMP
#pragma omp atomic
#endif
	_funcCallInfo[LOGLIKELIHOOD].totalTime += timerTot.RealTime();

	return std::log(intensity);
}


void
rpwa::importanceSampler::CalculateObservables(const std::vector<double>& parameters)
{
#ifdef _OPENMP
#pragma omp atomic
#endif
	++(_funcCallInfo[CALCULATEOBSERVABLES].nmbCalls);
	TStopwatch timerTot;
	timerTot.Start();

	rpwa::nBodyPhaseSpaceKinematics nBodyPhaseSpace;
	if (not initializeNBodyPhaseSpace(nBodyPhaseSpace, parameters, true)) {
		printErr << "error while initializing n-body phase space. Aborting..." << std::endl;
		throw;
	}

	const double mX = parameters[3*_nPart-6];
	nBodyPhaseSpace.calcBreakupMomenta();
	nBodyPhaseSpace.calcEventKinematics(TLorentzVector(0., 0., 0., mX));

	// calculate the squared masses of all possible two-body sub-systems
	// (used by BAT to create plots)
	size_t countPair = 0;
	for (size_t partA = 0; partA < _nPart; ++partA) {
		for (size_t partB = 0; partB < partA; ++partB) {
			GetObservable(countPair) = (nBodyPhaseSpace.daughter(partA) + nBodyPhaseSpace.daughter(partB)).M2();
			++countPair;
		}
	}

	// only store events if the file writer has been initialised and BAT is
	// no longer in the pre-run phase
	if (_fileWriter.initialized() and GetPhase() > 0) {
		// apply lag:
		// the 'CalculateObservables' method might be called more often
		// than the prior and log-likelihood functions in case points
		// outside the parameter limits are proposed inside BAT. to
		// protect from this, but also to reduce the impact of
		// correlations between two consecutive points, only store
		// every (fMCMCNLag)th point
		if (fMCMCCurrentIteration % fMCMCNLag != 0) {
			return;
		}

		// generate the production kinematics of the current event
		bool           validBeamAndTPrime;
		double         tPrime;
		TLorentzVector pBeam;
		TLorentzVector pX;
#ifdef _OPENMP
#pragma omp critical(prodKin)
#endif
		{
			boost::tuples::tie(validBeamAndTPrime, tPrime, pBeam, pX) = getProductionKinematics(mX);
		}
		if (not validBeamAndTPrime) {
			printErr << "vertex and beam generation or picking of tPrime failed. Aborting..." << std::endl;
			throw;
		}

		// boost the event into the laboratory system
		const TLorentzRotation toLab = rpwa::isobarAmplitude::gjTransform(pBeam, pX).Inverse();
		std::vector<TVector3> decayKinMomenta(_nPart);
		for (size_t part = 0; part < _nPart; ++part) {
			TLorentzVector p = nBodyPhaseSpace.daughter(part);
			p.Transform(toLab);
			decayKinMomenta[part] = p.Vect();
		}

		std::vector<double> additionalVars;
		if (_storeMassAndTPrime) {
			additionalVars.push_back(pX.M());
			additionalVars.push_back(tPrime);
		}
		std::vector<TVector3> prodKinMomenta(1, pBeam.Vect());
#ifdef _OPENMP
#pragma omp critical(fileWriter)
#endif
		{
			_fileWriter.addEvent(prodKinMomenta, decayKinMomenta, additionalVars);
		}
	}

	timerTot.Stop();
#ifdef _OPENMP
#pragma omp atomic
#endif
	_funcCallInfo[CALCULATEOBSERVABLES].totalTime += timerTot.RealTime();
}


bool
rpwa::importanceSampler::initializeFileWriter(TFile*             outFile,
                                              const std::string& auxString,
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

	multibinBoundariesType binning;
	binning[massVariableName]   = _massAndTPrimePicker->massRange();
	binning[tPrimeVariableName] = _massAndTPrimePicker->tPrimeRange();
	std::vector<std::string> additionalVarLabels;
	if (_storeMassAndTPrime) {
		additionalVarLabels.push_back(massVariableName);
		additionalVarLabels.push_back(tPrimeVariableName);
	}

	const bool valid = _fileWriter.initialize(*outFile,
	                                          auxString,
	                                          rpwa::eventMetadata::REAL,
	                                          prodKinParticleNames,
	                                          decayKinParticleNames,
	                                          binning,
	                                          additionalVarLabels);

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
		printWarn << "could not generate vertex and beam." << std::endl;
		return boost::tuples::make_tuple(false, 0., TLorentzVector(), TLorentzVector());
	}
	const TLorentzVector targetLab(0., 0., 0., _target.targetParticle.mass());
	const TLorentzVector beamLorentzVector = _beamAndVertexGenerator->getBeam();
	const TLorentzVector overallCm = beamLorentzVector + targetLab;  // beam-target center-of-mass system

	double tPrime;
	if (not _massAndTPrimePicker->pickTPrimeForMass(xMass, tPrime)) {
		printWarn << "t' pick failed." << std::endl;
		return boost::tuples::make_tuple(false, 0., TLorentzVector(), TLorentzVector());
	}

	// get Lorentz vector of X in lab system
	const TLorentzVector xSystemLab = rpwa::diffractivePhaseSpace::calculateXSystemLab(beamLorentzVector,
	                                                                                   _target.targetParticle.mass(),
	                                                                                   xMass,
	                                                                                   _target.recoilParticle.mass(),
	                                                                                   overallCm.M2(), tPrime);

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
