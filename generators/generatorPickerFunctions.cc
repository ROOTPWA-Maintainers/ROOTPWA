
#include<assert.h>
#include<map>

#include<boost/assign/std/vector.hpp>

#include "randomNumberGenerator.h"
#include "libConfigUtils.hpp"
#include "generatorPickerFunctions.h"


using namespace boost::assign;
using namespace libconfig;
using namespace std;
using namespace rpwa;


void massAndTPrimePicker::overrideMassRange(double lowerLimit, double upperLimit) {
	if(not _initialized) {
		printErr << "cannot override massRange on uninitialized massAndTPrimePicker." << endl;
		throw;
	}
	_massRange.first = lowerLimit;
	_massRange.second = upperLimit;
	printInfo << "mass range overriden: new range ]" << _massRange.first << ", " << _massRange.second << "]." << endl;
}


const std::pair<double, double>& massAndTPrimePicker::massRange() const {
	if(not _initialized) {
		printErr << "cannot call massRange on uninitialized massAndTPrimePicker." << endl;
		throw;
	}
	return _massRange;
}


bool massAndTPrimePicker::initTPrimeAndMassRanges(const libconfig::Setting& setting) {
	if(_initialized) {
		printErr << "trying to initialize a massAndTPrimePicker class twice." << endl;
		return false;
	}
	map < string, Setting::Type > mandatoryArguments;
	insert (mandatoryArguments)
		("massMin", Setting::TypeFloat)
		("massMax", Setting::TypeFloat);
	if(not checkIfAllVariablesAreThere(&setting, mandatoryArguments)) {
		printErr << "found invalid settings for the mass range of a mass and t' picker." << endl;
		return false;
	}
	_massRange.first = setting["massMin"];
	_massRange.second = setting["massMax"];
	if(_massRange.second < _massRange.first) {
		printErr << "'massMax' must not be smaller than 'massMin'." << endl;
		return false;
	}
	if(setting.exists("tPrimeMin")) {
		if(not setting.lookupValue("tPrimeMin", _tPrimeRange.first)) {
			printWarn << "'tPrimeMin' setting is invalid. Setting 'tPrimeMin' to "
			          << _tPrimeRange.first << "." << endl;
		}
	} else {
		printInfo << "'tPrimeMin' not specified. Setting it to "
		          << _tPrimeRange.first << "." << endl;
	}
	if(setting.exists("tPrimeMax")) {
		if(not setting.lookupValue("tPrimeMax", _tPrimeRange.second)) {
			printWarn << "'tPrimeMax' setting is invalid. Setting 'tPrimeMax' to "
			          << _tPrimeRange.second << "." << endl;
		}
	} else {
		printInfo << "'tPrimeMax' not specified. Setting it to "
		          << _tPrimeRange.second << "." << endl;
	}
	if(_tPrimeRange.second < _tPrimeRange.first) {
		printErr << "'tPrimeMax' must not be smaller than 'tPrimeMax'."
		         << endl;
		return false;
	}
	return true;
}


uniformMassExponentialTPicker::uniformMassExponentialTPicker()
	: massAndTPrimePicker(),
	  _nExponential(0) { }


uniformMassExponentialTPicker::uniformMassExponentialTPicker(const uniformMassExponentialTPicker& picker)
	: massAndTPrimePicker(picker),
	  _tSlopesForMassBins(picker._tSlopesForMassBins),
	  _nExponential(picker._nExponential) { }


bool uniformMassExponentialTPicker::init(const Setting& setting) {
	if(not massAndTPrimePicker::initTPrimeAndMassRanges(setting)) {
		printErr << "could not initialize t' or mass range settings in 'uniformMassExponentialT'." << endl;
		return false;
	}
	map < string, Setting::Type > mandatoryArguments;
	insert (mandatoryArguments)
		("invariantMasses", Setting::TypeArray)
		("tSlopes", Setting::TypeList);
	if(not checkIfAllVariablesAreThere(&setting, mandatoryArguments)) {
		printErr << "found an invalid setting for function 'uniformMassExponentialT'." << endl;
		return false;
	}
	if(setting["invariantMasses"].getLength() != setting["tSlopes"].getLength()) {
		printErr << "'invariantMasses' and 'tSlopes' have to be arrays of the same length." << endl;
		return false;
	}
	if(not setting["invariantMasses"][0].isNumber()) {
		printErr << "'invariantMasses' has to be array of numbers." << endl;
		return false;
	}
	if(not (setting["tSlopes"][0].isNumber()
			|| (setting["tSlopes"][0].isArray() && setting["tSlopes"][0][0].isNumber()))) {
		printErr << "'tSlopes' has to be number or array of numbers." << endl;
		return false;
	}
	bool tArray = false;
	if(setting["tSlopes"][0].isNumber()) {
		_nExponential = 1;
	} else {
		if(setting["tSlopes"][0].getLength()%2 != 1) {
			printErr << "number of parameters in arrays in 'tSlopes' has to be odd." << endl;
			return false;
		}
		_nExponential = (setting["tSlopes"][0].getLength()+1) / 2;
		tArray = true;
	}
	for(int i = 0; i < setting["invariantMasses"].getLength(); ++i) {
		vector<double> param;
		if(tArray && not setting["tSlopes"][i].isArray()) {
			printErr << "all entries of 'tSlopes' have to be either numbers or arrays of numbers." << endl;
			return false;
		}
		if(tArray) {
			unsigned int nExponential = (setting["tSlopes"][i].getLength()+1) / 2;
			if(nExponential != _nExponential) {
				printErr << "all entries of 'tSlopes' have to be arrays of the same length" << endl;
				return false;
			}
			for(int j = 0; j < setting["tSlopes"][i].getLength(); ++j) {
				param.push_back(setting["tSlopes"][i][j]);
			}
		} else {
			param.push_back(setting["tSlopes"][i]);
		}
		// check parameters for correctness
		for(unsigned int j = 0; j < _nExponential; ++j) {
			if(param[2*j] >= 0.) {
				const double mass = setting["invariantMasses"][i];
				printErr << "while reading 'tSlopes' for invariant mass " << mass << ": "
				         << "parameter " << 2*j << " (current value = " << param[2*j] << ") has to be smaller than zero." << endl;
				return false;
			}
			if(j > 0) {
				if(param[2*j-1] < 0.) {
					const double mass = setting["invariantMasses"][i];
					printErr << "while reading 'tSlopes' for invariant mass " << mass << ": "
					         << "parameter " << 2*j-1 << " (current value = " << param[2*j-1] << ") has to be larger than or equal to zero." << endl;
					return false;
				}
			}
		}
		_tSlopesForMassBins.insert(pair<double, vector<double> >(setting["invariantMasses"][i], param));
	}
	_initialized = true;
	return true;
}


bool uniformMassExponentialTPicker::operator() (double& invariantMass, double& tPrime) {
	if(not _initialized) {
		printErr << "trying to use an uninitialized massAndTPrimePicker." << endl;
		return false;
	}
	if(_tSlopesForMassBins.empty()) {
		printErr << "no t' slopes to generate t' from." << endl;
		return false;
	}
	TRandom3* randomNumbers = randomNumberGenerator::instance()->getGenerator();
	invariantMass = randomNumbers->Uniform(_massRange.first, _massRange.second);
	vector<double> param;
	if(_tSlopesForMassBins.size() == 1) {
		param = _tSlopesForMassBins.begin()->second;
	} else {
		if(invariantMass > _tSlopesForMassBins.rbegin()->first) {
			param = _tSlopesForMassBins.rbegin()->second;
		} else if(invariantMass < _tSlopesForMassBins.begin()->first) {
			param = _tSlopesForMassBins.begin()->second;
		} else {
			map<double, vector<double> >::const_iterator itUpper = _tSlopesForMassBins.lower_bound(invariantMass);
			map<double, vector<double> >::const_iterator itLower = itUpper; --itLower;
			const double mass1 = itLower->first;
			const double mass2 = itUpper->first;
			assert(itLower->second.size() == itUpper->second.size());
			for(size_t i = 0; i < itLower->second.size(); i++) {
				const double param1 = itLower->second[i];
				const double param2 = itUpper->second[i];
				param.push_back(param1 + ((param2-param1) / (mass2-mass1) * (invariantMass-mass1)));
			}
		}
	}
	if(((param.size()+1) / 2) != _nExponential) {
		printErr << "error when calculating the parameters for t'-slope." << endl;
		return false;
	}
        // short-cut for one exponential, then t can analytically be calculated
	if (_nExponential == 1) {
		const double r = randomNumbers->Uniform();
		const double Fmin = exp(param[0] * _tPrimeRange.first);
		const double Fmax = exp(param[0] * _tPrimeRange.second);

		tPrime = log(r*(Fmax-Fmin) + Fmin) / param[0];

		return true;
	}
	// search for the t' using the Newton-Raphson method
	// way faster than using ROOT TF1s
	// the function to find the root for is
	//          Int(f(t'),t'=(tmin..t))
	// g(t) = -------------------------- - Rndm
	//        Int(f(t'),t'=(tmin..tmax))
	// with
	// f(t) = exp(p[0]*t) + Sum(p[2*i-1] * exp(p[2*i]*t),i=1.._nExponential)
	//
	// The idea behind this is the usual ansatz of getting one random
	// number r between 0 and 1, and then find the t for which the integral
	// from tmin to t is equal to the integral from tmin to tmax times r.
	// t is then distributed according to the function.
	//
	// The procedure here works only if the function is strictly monotonic
	// decreasing, maybe that should be checked in init().
	double Fmin = exp(param[0] * _tPrimeRange.first) / param[0];
	double Fmax = exp(param[0] * _tPrimeRange.second) / param[0];
	for(unsigned int i = 1; i < _nExponential; ++i) {
		Fmin += param[2*i-1] * exp(param[2*i] * _tPrimeRange.first) / param[2*i];
		Fmax += param[2*i-1] * exp(param[2*i] * _tPrimeRange.second) / param[2*i];
	}
	bool done = false;
	unsigned int restarts = 0;
	while (not done) {
		tPrime = 0.;
		const double r = randomNumbers->Uniform();
		unsigned int steps = 0;
		while(not done) {
			// calculate g (here: F) and g' (here: f)
			double Ft = exp(param[0] * tPrime) / param[0];
			double ft = exp(param[0] * tPrime);
			for(unsigned int i = 1; i < _nExponential; ++i) {
				Ft += param[2*i-1] * exp(param[2*i] * tPrime) / param[2*i];
				ft += param[2*i-1] * exp(param[2*i] * tPrime);
			}
			// check whether the current value for tPrime is already okay
			const double errF = r*(Fmax-Fmin)*numeric_limits<double>::epsilon();
			const double diffF = (Ft-Fmin) - r*(Fmax-Fmin);
			if(abs(diffF) < 10.*errF) {
				done = true;
				break;
			}
			const double errU = tPrime*numeric_limits<double>::epsilon();
			const double update = (Ft - r*Fmax - (1-r)*Fmin) / ft;
			if(abs(update) < 10.*errU) {
				done = true;
				break;
			}
			if(update > 0.) {
				// for this strictly monotonic decreasing
				// function with a strictly monotonic decreaing
				// derivative, a step should never be done
				// "backwards" (if it happens it should
				// typically be a sign for numeric limitations)
				done = true;
				break;
			}
			if(steps >= 50) {
				printWarn << "t' could not be found in " << steps << " steps:" << endl
				          << "\t\t r*(Fmax-Fmin) = " << r*(Fmax-Fmin) << endl
				          << "\t\t     (Ft-Fmin) = " << (Ft-Fmin) << endl
				          << "\t\t   actual diff = " << diffF << endl
				          << "\t\t  allowed diff = " << errF << endl
				          << "\t\t actual update = " << update << endl
				          << "\t\tminimal update = " << errU << endl
				          << "\t\t        tPrime = " << tPrime << endl
				          << "restarting with a new random number." << endl;
				done = false;
				break;
			}
			// do the step to the next iteration
			tPrime -= update;
			steps++;
		}
		restarts++;
		if(restarts >= 5 && not done) {
                            printErr << restarts << " attempts to generate t' failed." << endl;
			return false;
		}
	}
	return true;
}


ostream& uniformMassExponentialTPicker::print(ostream& out) const {
	out << "'uniformMassExponentialT' weighter parameters:" << endl;
	out << "    minimum Mass ... " << _massRange.first << endl;
	out << "    maximum Mass ... " << _massRange.second << endl;
	out << "    minimum t' ..... " << _tPrimeRange.first << endl;
	out << "    maximum t' ..... " << _tPrimeRange.second << endl;
	if(_tSlopesForMassBins.size() == 0) {
		out << "    t' slope ....... []";
	} else if(_tSlopesForMassBins.size() == 1) {
		out << "    t' slope ....... ";
		printSlice(out, _tSlopesForMassBins.begin()->second);
		out << endl;
	} else {
		map<double, vector<double> >::const_iterator it = _tSlopesForMassBins.begin();
		const map<double, vector<double> >::const_iterator itLast = --(_tSlopesForMassBins.end());
		out << "    t' slopes for interpolation:" << endl;
		out << "        mass < " << setw(7) << left << it->first
		    << " -> t' slope = ";
		printSlice(out, it->second);
		out << endl;
		while(++it != itLast) {
			stringstream strStr;
			strStr << "        t' slope(mass = " << it->first << ") ";
			out << setw(35) << left << strStr.str()
			    << setw(0) << "= ";
			printSlice(out, it->second);
			out << endl;
		}
		out << "        mass >= " << setw(6) << left << it->first
		    << " -> t' slope = ";
		printSlice(out, it->second);
		out << endl;
	}
	return out;
}


ostream& uniformMassExponentialTPicker::printSlice(ostream& out, const vector<double>& param) const {
	out << "[";
	for(size_t i = 0; i < param.size(); ++i) {
		if(i != 0) {
			out << ",";
		}
		out << param[i];
	}
	out << "]";
	return out;
}


polynomialMassAndTPrimeSlopePicker::polynomialMassAndTPrimeSlopePicker()
	: massAndTPrimePicker() { }


polynomialMassAndTPrimeSlopePicker::polynomialMassAndTPrimeSlopePicker(const polynomialMassAndTPrimeSlopePicker& picker)
	: massAndTPrimePicker(picker),
	  _massPolynomial(picker._massPolynomial),
	  _tPrimeSlopePolynomial(picker._tPrimeSlopePolynomial) { }


bool polynomialMassAndTPrimeSlopePicker::init(const Setting& setting) {
	if(not massAndTPrimePicker::initTPrimeAndMassRanges(setting)) {
		printErr << "could not initialize t' or mass range settings in 'polynomialMassAndTPrime'." << endl;
		return false;
	}
	map < string, Setting::Type > mandatoryArguments;
	insert (mandatoryArguments)
		("coeffsMass", Setting::TypeArray)
		("coeffsTSlopes", Setting::TypeArray);
	if(not checkIfAllVariablesAreThere(&setting, mandatoryArguments)) {
		printErr << "found an invalid setting for function 'polynomialMassAndTPrime'." << endl;
		return false;
	}
	const Setting& configCoeffsMass = setting["coeffsMass"];
	unsigned int numberOfMassCoeffs = configCoeffsMass.getLength();
	if(numberOfMassCoeffs < 1) {
		printErr << "'coeffsMass' array must not be empty." << endl;
		return false;
	}
	if(not configCoeffsMass[0].isNumber()) {
		printErr << "'coeffsMass' array has to be made up of numbers." << endl;
		return false;
	}
	stringstream strStr;
	strStr << "pol" << (numberOfMassCoeffs - 1);
	_massPolynomial = TF1("generatorPickerMassPolynomial", strStr.str().c_str(), _massRange.first, _massRange.second);
	for(unsigned int i = 0; i < numberOfMassCoeffs; ++i) {
		_massPolynomial.SetParameter(i, configCoeffsMass[i]);
	}
	const Setting& configCoeffsTSlopes = setting["coeffsTSlopes"];
	unsigned int numberOfTSlopeCoeffs = configCoeffsTSlopes.getLength();
	if(numberOfTSlopeCoeffs < 1) {
		printErr << "'coeffsTSlopes' array must not be empty." << endl;
		return false;
	}
	if(not configCoeffsTSlopes[0].isNumber()) {
		printErr << "'coeffsTSlopes' array has to be made up of numbers." << endl;
		return false;
	}
	strStr.str("");
	strStr << "pol" << (numberOfTSlopeCoeffs - 1);
	_tPrimeSlopePolynomial = TF1("generatorPickerMassPolynomial", strStr.str().c_str());
	for(unsigned int i = 0; i < numberOfTSlopeCoeffs; ++i) {
		_tPrimeSlopePolynomial.SetParameter(i, configCoeffsTSlopes[i]);
	}
	_initialized = true;
	return true;
}


bool polynomialMassAndTPrimeSlopePicker::operator()(double& invariantMass, double& tPrime) {
	if(not _initialized) {
		printErr << "trying to use an uninitialized massAndTPrimePicker." << endl;
		return false;
	}
	TRandom3* randomNumbers = randomNumberGenerator::instance()->getGenerator();
	invariantMass = _massPolynomial.GetRandom(_massRange.first, _massRange.second);
	double tPrimeSlope = _tPrimeSlopePolynomial.Eval(invariantMass);
	do {
		tPrime = randomNumbers->Exp(1. / tPrimeSlope);
	} while(tPrime < _tPrimeRange.first || tPrime > _tPrimeRange.second);
	return true;
}


ostream& polynomialMassAndTPrimeSlopePicker::print(ostream& out) const {
	out << "'polynomialMassAndTPrime' weighter parameters:" << endl;
	out << "    minimum Mass ........... " << _massRange.first << endl;
	out << "    maximum Mass ........... " << _massRange.second << endl;
	out << "    mass polynomial ........ " << _massPolynomial.GetExpFormula() <<endl;
	out << "    t' slopes polynomial ... " << _tPrimeSlopePolynomial.GetExpFormula() <<endl;
	return out;
}
