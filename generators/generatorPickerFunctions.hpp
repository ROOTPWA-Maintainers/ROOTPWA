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

		uniformMassExponentialTPicker() { }
		virtual ~uniformMassExponentialTPicker() { }

		virtual bool init(const libconfig::Setting& setting) {
			std::map<std::string, libconfig::Setting::Type> mandatoryArguments;
			boost::assign::insert (mandatoryArguments)
				("mass_min", libconfig::Setting::TypeFloat)
			    ("mass_max", libconfig::Setting::TypeFloat)
			    ("t_slope", libconfig::Setting::TypeArray)
			    ("inv_m", libconfig::Setting::TypeArray)
			    ("t_min", libconfig::Setting::TypeFloat);
			if(not checkIfAllVariablesAreThere(&setting, mandatoryArguments)) {
				printErr << "found an invalid settings for function 'uniformMassExponentialT'." << std::endl;
				return false;
			}
			_minimumMass = setting["mass_min"];
			_maximumMass = setting["mass_max"];
			_minimumTPrime = setting["t_min"];
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
			return true;
		}

		virtual bool operator() (double& invariantMass, double& tPrime) {
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
			while(tPrime <= _minimumTPrime) {
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

		double _minimumMass;
		double _maximumMass;
		double _minimumTPrime;

	};

}

#endif
