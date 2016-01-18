#include<cmath>
#include<iostream>
#include<sstream>

#include<TFile.h>

#include<BAT/BCMath.h>

#include"importanceSampler.h"
#include"physUtils.hpp"
#include"randomNumberGenerator.h"


size_t rpwa::importanceSampler::_nCalls = 0;


rpwa::importanceSampler::importanceSampler(const double            mMin,
                                           const double            mMax,
                                           rpwa::modelIntensityPtr model)
	: BCModel("phaseSpaceImportanceSampling"),
	  _phaseSpaceOnly(false),
	  _productionGeneratorInitialized(false),
	  _zeroBinWidth(false),
	  _nPart(0),
	  _mMin(mMin),
	  _mMax(mMax),
	  _masses(model->finalStateMasses()),
	  _massPrior(0),
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
	AddParameter("phi1", 0., 2*M_PI, "#phi_{1}");
	AddParameter("cosTheta1", -1., 1., "cos(#theta_{1})"); // sample cos(theta), because dOmega = sin(theta) dtheta dphi = dcos(theta) dphi
	double mSum = 0;
	for (size_t p = 1; p < _nPart; ++p) { // The first particle never appears in the mSum
		mSum += _masses[p];
	}
	for (size_t part = 2; part < _nPart; ++part) {
		std::stringstream nameM;
		std::stringstream namePhi;
		std::stringstream nameTheta;
		std::stringstream texM;
		std::stringstream texPhi;
		std::stringstream texTheta;

		nameM << "m";
		texM << "m_{";
		for (size_t p = part; p < _nPart+1; ++p) {
			nameM << p;
			texM << p;
		}
		texM << "}";
		namePhi << "phi" << part;
		texPhi << "#phi_{" << part << "}";
		nameTheta << "cosTheta" << part; // sample cos(theta), because dOmega = sin(theta) dtheta dphi = dcos(theta) dphi
		texTheta << "cos(#theta_{" << part << "})";
		AddParameter(nameM.str(), mSum, _mMax, texM.str());
		mSum -= _masses[part-1];
		AddParameter(namePhi.str(), 0., 2*M_PI, texPhi.str());
		AddParameter(nameTheta.str(), -1., 1., texTheta.str());
	}
	if (_mMin < _mMax) {
		AddParameter("mass", _mMin, _mMax, "m_{X}");
	} else if (_mMin == _mMax) {
		printInfo << "zeroBinWidth, do not sample mX" << std::endl;
		_zeroBinWidth = true;
	} else {
		printErr << "_mMin = " << _mMin << " > _mMax = " << _mMax << ". Aborting..." << std::endl;
		throw;
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
rpwa::importanceSampler::LogLikelihood(const std::vector<double>& parameters)
{
	++_nCalls;
	if (!_productionGeneratorInitialized) {
		printErr << "Production kinematics not initilized. Cannot generate beam." << std::endl;
		throw;
	}
	if (!_phaseSpaceOnly) {
		std::vector<TLorentzVector> fourVectors = getFinalStateMomenta(parameters);
		std::vector<TVector3> decayKinMomenta(_nPart);
		for (size_t part = 0; part < _nPart; ++part) {
			decayKinMomenta[part] = fourVectors[part].Vect();
		}
		const double intens = _model->getIntensity(decayKinMomenta);
		return std::log(intens);
	} else {
		return -1;
	}
}


double
rpwa::importanceSampler::LogAPrioriProbability(const std::vector<double>& parameters)
{
	// Just calculate non constant parts;
	double mass = _mMin;
	if (!_zeroBinWidth) {
		mass = parameters[parameters.size()-1];
		if (mass < _mMin || mass >= _mMax) {
			return -std::numeric_limits<double>::infinity();
		}
	}

	double phaseSpace = 1./pow(mass, 4); // Prior = 1/M^2 => calculate prior^2 -> 1/M^4
	if (_nPart == 2) { // Two particle phase space is constant. Handle this case special that part == 0 and part == _mParr-2 can appear laser.
		return std::log(phaseSpace)/2.;
	}

	for (size_t part = 0; part < _nPart - 1; ++part) {
		double M, m1, m2;
		if (part == 0) {
			M  = mass;
			m1 = _masses[0];
			m2 = parameters[2];
		} else if (part == _nPart-2) {
			M  = parameters[3*part-1];
			m1 = _masses[part];
			m2 = _masses[part+1];
		} else {
			M  = parameters[3*part-1];
			m1 = _masses[part];
			m2 = parameters[3*part+2];
		}
		if (M < m1 + m2) {
			return -std::numeric_limits<double>::infinity();
		}
		phaseSpace *= breakupMomentumSquared(M, m1, m2, true);
		if (phaseSpace < 0.) {
			return -std::numeric_limits<double>::infinity();
		}
	}
	if (_massPrior) {
		double massPrior = _massPrior->Eval(mass);
		phaseSpace *= massPrior*massPrior; // massPrior^2 since everything is sqared up to now
	}
	return std::log(phaseSpace)/2; // The product above is actually multiplied with the square of the phase space, positivity is ensured, so the missung sqrt becomes a factor 1/2 when pulling it out of the log
}


void
rpwa::importanceSampler::CalculateObservables(const std::vector<double>& parameters)
{
	std::vector<TLorentzVector> vectors = getFinalStateMomenta(parameters);
	size_t count = 0;
	for (size_t part = 0; part < _nPart; ++part) {
		for (size_t qart = 0; qart < part; ++qart) {
			GetObservable(count) = (vectors[part] + vectors[qart]).M2();
			++count;
		}
	}
	if (_fileWriter.initialized() && GetPhase() > 0) { // Put file writer stuff here
		if (fMCMCCurrentIteration % fMCMCNLag != 0) {
			return; // apply lag
		}
		std::vector<TVector3> decayKinMomenta(vectors.size());
		double mass = _mMin;
		if (!_zeroBinWidth) {
			mass = parameters[parameters.size()-1];
		}

		std::pair<std::pair<bool, double>, std::vector<TLorentzVector> > production = getProductionKinematics(mass);
		if (!production.first.first) {
			printErr << "generation of production failed" << std::endl;
			throw;
		}
		TLorentzVector pX;
		std::pair<bool, TLorentzRotation> inverseGJ = getInverseGJTransform(production.second[0], production.second[1]);
		for (size_t part = 0; part < vectors.size(); ++part) {
			vectors[part].Transform(inverseGJ.second);
			decayKinMomenta[part] = vectors[part].Vect();
			pX+=vectors[part];
		}
		std::vector<double>var(3+_nPart);
		var[0] = production.first.second;
		var[1] = pX.M();
		var[2] = pX.E();
		for (size_t part = 0; part < _nPart; ++part) {
			var[part + 3] = vectors[part].E();
		}
		std::vector<TVector3> prodKinMomenta(1, production.second[0].Vect());
		_fileWriter.addEvent(prodKinMomenta, decayKinMomenta, var);
	}
}


std::vector<TLorentzVector>
rpwa::importanceSampler::getFinalStateMomenta(const std::vector<double>& samplingParameters) const
{
	std::vector<TLorentzVector> retVal(_nPart);
	double MX = _mMin;
	if (!_zeroBinWidth) {
		MX = samplingParameters[samplingParameters.size()-1];
	}
	if (_nPart == 2) {
		double q = rpwa::breakupMomentum(MX, _masses[0], _masses[1]);
		double E1 = pow(q*q+_masses[0]*_masses[0], .5);
		double E2 = pow(q*q+_masses[1]*_masses[1], .5);
		double cosT = samplingParameters[1];
		double sinT = sqrt(1. - cosT*cosT); // Can be done, since theta in [0, pi] -> sin(theta) >= 0.
		retVal[0] = TLorentzVector(q*cos(samplingParameters[0])*sinT,
		                           q*sin(samplingParameters[0])*sinT,
		                           q*cosT,
		                           E1);
		retVal[1] = TLorentzVector(-q*cos(samplingParameters[0])*sinT,
		                           -q*sin(samplingParameters[0])*sinT,
		                           -q*cosT,
		                            E2);
		return retVal;
	}

	TLorentzVector isobarMomentum(0., 0., 0., MX);
	double q  = 0.;
	double M  = 0.;
	double m1 = 0.;
	double m2 = 0.;
	double E1 = 0.;
	double E2 = 0.;
	for (size_t part = 0; part < _nPart-1; ++part) {
		if (part == 0) {
			M  = MX;
			m1 = _masses[part];
			m2 = samplingParameters[3*part+2];
		} else if (part == _nPart-2) {
			M  = samplingParameters[3*part-1];
			m1 = _masses[part  ];
			m2 = _masses[part+1];
		} else {
			M  = samplingParameters[3*part-1];
			m1 = _masses[part];
			m2 = samplingParameters[3*part+2];
		}
		q  = rpwa::breakupMomentumSquared(M, m1, m2, true);
		if (q < 0.) {
			q = 0.;
		} else {
			q = sqrt(q);
		}
		E1 = pow(q*q+m1*m1, .5);
		E2 = pow(q*q+m2*m2, .5);

		TVector3 boostVector = isobarMomentum.BoostVector();
		double cosT = samplingParameters[3*part+1];
		double sinT = sqrt(1. - cosT*cosT); // Can be done, since theta in [0, pi] -> sin(theta) >= 0.
		retVal[part]         = TLorentzVector(q*cos(samplingParameters[3*part])*sinT,
		                                      q*sin(samplingParameters[3*part])*sinT,
		                                      q*cosT,
		                                      E1);
		retVal[part].Boost(boostVector);
		isobarMomentum = TLorentzVector(-q*cos(samplingParameters[3*part])*sinT,
		                                -q*sin(samplingParameters[3*part])*sinT,
		                                -q*cosT,
		                                 E2);
		isobarMomentum.Boost(boostVector);
		if (part == _nPart-2) {
			retVal[part+1] = isobarMomentum;
		}
	}
	return retVal;
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


std::pair<bool, TLorentzRotation>
rpwa::importanceSampler::getInverseGJTransform(const TLorentzVector& beamLv,
                                               const TLorentzVector& XLv)
{
	TLorentzVector beam = beamLv;
	TLorentzVector X    = XLv;
	const TVector3 yGjAxis = beam.Vect().Cross(X.Vect()); // y-axis of Gottfried-Jackson frame
	// rotate so that yGjAxis becomes parallel to y-axis and beam momentum ends up in (x, z)-plane
	TRotation rot1;
	rot1.RotateZ(piHalf - yGjAxis.Phi());
	rot1.RotateX(yGjAxis.Theta() - piHalf);
	beam *= rot1;
	X    *= rot1;
	// boost to X RF
	TLorentzRotation boost;
	boost.Boost(-X.BoostVector());
	beam *= boost;
	// rotate about yGjAxis so that beam momentum is along z-axis
	TRotation rot2;
	rot2.RotateY(-signum(beam.X()) * beam.Theta());
	// construct total transformation
	TLorentzRotation gjTransform(rot1);
	gjTransform.Transform(boost);
	gjTransform.Transform(rot2);
	return std::pair<bool, TLorentzRotation>(true, gjTransform.Inverse());
}
