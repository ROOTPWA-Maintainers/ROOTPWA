#ifndef GENERATORWEIGHTFUNCTIONS_HH_
#define GENERATORWEIGHTFUNCTIONS_HH_


#include<map>

#include<boost/assign/std/vector.hpp>

#include<libconfig.h++>

#include "randomNumberGenerator.h"
#include "libConfigUtils.hpp"


namespace rpwa {

	class massAndTPrimePicker {

	  public:

		massAndTPrimePicker() { };
		virtual ~massAndTPrimePicker() { };

		virtual bool init(const libconfig::Setting& setting) = 0;

		virtual bool operator() (double& invariantMass, double& tPrime) = 0;

		virtual std::ostream& print(std::ostream& out) = 0;

	};

	class uniformMassExponentialTPicker : public massAndTPrimePicker {

	  public:

		uniformMassExponentialTPicker()
		  : _initialized(false),
		    _minimumTPrime(0.),
		    _maximumTPrime(std::numeric_limits<double>::max()) { }
		virtual ~uniformMassExponentialTPicker() { }

		virtual bool init(const libconfig::Setting& setting) {
			if(_initialized) {
				printErr << "trying to initialize a massAndTPrimePicker class twice." << std::endl;
				return false;
			}
			std::map<std::string, libconfig::Setting::Type> mandatoryArguments;
			boost::assign::insert (mandatoryArguments)
				("mass_min", libconfig::Setting::TypeFloat)
			    ("mass_max", libconfig::Setting::TypeFloat)
			    ("t_slope", libconfig::Setting::TypeArray)
			    ("inv_m", libconfig::Setting::TypeArray);
			if(not checkIfAllVariablesAreThere(&setting, mandatoryArguments)) {
				printErr << "found an invalid settings for function 'uniformMassExponentialT'." << std::endl;
				return false;
			}
			_minimumMass = setting["mass_min"];
			_maximumMass = setting["mass_max"];
			if(_maximumMass < _minimumMass) {
				printErr << "'mass_max' must not be smaller than 'mass_min'." << std::endl;
				return false;
			}
			if(setting.exists("t_prime_min")) {
				if(not setting.lookupValue("t_prime_min", _minimumTPrime)) {
					printWarn << "'t_prime_min' setting is invalid. Setting 't_prime_min' to "
					          << _minimumTPrime << "." << std::endl;
				}
			} else {
				printInfo << "'t_prime_min' not specified. Setting it to "
				          << _minimumTPrime << "." << std::endl;
			}
			if(setting.exists("t_prime_max")) {
				if(not setting.lookupValue("t_prime_max", _maximumTPrime)) {
					printWarn << "'t_prime_max' setting is invalid. Setting 't_prime_max' to "
					          << _maximumTPrime << "." << std::endl;
				}
			} else {
				printInfo << "'t_prime_max' not specified. Setting it to "
				          << _maximumTPrime << "." << std::endl;
			}
			if(_maximumTPrime < _minimumTPrime) {
				printErr << "'t_prime_max' must not be smaller than 't_prime_min'." << std::endl;
				return false;
			}
			if(setting["inv_m"].getLength() != setting["t_slope"].getLength()) {
				printErr << "'inv_m' and 't_slope' have to be arrays of the same length." << std::endl;
				return false;
			}
			if(not (setting["inv_m"][0].isNumber()) && (setting["t_slope"][0].isNumber())) {
				printErr << "'inv_m' and 't_slope' have to be array of numbers." << std::endl;
				return false;
			}
			for(int i = 0; i < setting["inv_m"].getLength(); ++i) {
				_tSlopesForMassBins.push_back(std::pair<double, double>(setting["inv_m"][i], setting["t_slope"][i]));
			}
			_initialized = true;
			return true;
		}

		virtual bool operator() (double& invariantMass, double& tPrime) {
			if(not _initialized) {
				printErr << "trying to use an uninitialized massAndTPrimePicker." << std::endl;
				return false;
			}
			TRandom3* randomNumbers = rpwa::randomNumberGenerator::instance()->getGenerator();
			invariantMass = randomNumbers->Uniform(_minimumMass, _maximumMass);
			double tPrimeSlope = -1.;
			if(_tSlopesForMassBins.size() == 1) {
				tPrimeSlope = _tSlopesForMassBins[0].second;
			} else {
				if(invariantMass > _tSlopesForMassBins[_tSlopesForMassBins.size()-1].first) {
					tPrimeSlope = _tSlopesForMassBins[_tSlopesForMassBins.size()-1].second;
				} else {
					unsigned int i = 0;
					for(; (invariantMass > _tSlopesForMassBins[i].first); ++i);
					tPrimeSlope = _tSlopesForMassBins[i].second;
				}
			}
			if(tPrimeSlope < 0) {
				printErr << "error when calculating the t'-slope." << std::endl;
				return false;
			}
			tPrime = _minimumTPrime - 1.;
			while(tPrime < _minimumTPrime || tPrime > _maximumTPrime) {
				tPrime = randomNumbers->Exp(tPrimeSlope);
			}
			printDebug << "generated mass=" << invariantMass << " | t'Slope=" << tPrimeSlope << std::endl;
			return true;
		}

		virtual std::ostream& print(std::ostream& out) {
			out << "'uniformMassExponentialT' weighter parameters:" << std::endl;
			out << "    minimum Mass ... " << _minimumMass << std::endl;
			out << "    maximum Mass ... " << _maximumMass << std::endl;
			out << "    minimum t' ..... " << _minimumTPrime << std::endl;
			out << "    maximum t' ..... " << _maximumTPrime << std::endl;
			if(_tSlopesForMassBins.size() == 1) {
				out << "    t' slope ....... " << _tSlopesForMassBins[0].second << std::endl;
			} else {
				unsigned int nSlots = _tSlopesForMassBins.size();
				out << "    t' slopes:" << std::endl;
				out << "        mass < " << _tSlopesForMassBins[0].first << " -> t' slope = "
				    << _tSlopesForMassBins[0].second << std::endl;
				for(unsigned int i = 0; i < nSlots - 1; ++i) {
					out << "        mass in [" << _tSlopesForMassBins[i].first << ", "
					    << _tSlopesForMassBins[i+1].first << "] -> t' slope = "
					    << _tSlopesForMassBins[i+1].second << std::endl;
				}
				out << "        mass >= " <<  _tSlopesForMassBins[nSlots-1].first << " -> t' slope = "
				    << _tSlopesForMassBins[nSlots-1].second << std::endl;
			}
			return out;
		}

	  private:

		std::vector<std::pair<double, double> > _tSlopesForMassBins;

		bool _initialized;

		double _minimumMass;
		double _maximumMass;
		double _minimumTPrime;
		double _maximumTPrime;

	};

}

#endif
