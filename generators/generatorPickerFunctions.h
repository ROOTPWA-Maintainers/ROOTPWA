#ifndef GENERATORWEIGHTFUNCTIONS_HH_
#define GENERATORWEIGHTFUNCTIONS_HH_


#include<libconfig.h++>


namespace rpwa {

	class massAndTPrimePicker {

	  public:

		massAndTPrimePicker()
			: _massRange(std::pair<double, double>(0., 0.)),
			  _initialized(false) { };

		virtual ~massAndTPrimePicker() { };

		virtual bool init(const libconfig::Setting& setting) = 0;

		virtual std::pair<double, double> massRange() {
			if(not _initialized) {
				printErr << "cannot call massRange on uninitialized massAndTPrimePicker." << std::endl;
				throw;
			}
			return _massRange;
		};

		virtual bool operator() (double& invariantMass, double& tPrime) = 0;

		virtual std::ostream& print(std::ostream& out) = 0;

	  protected:

		std::pair<double, double> _massRange;
		bool _initialized;

	};

	class uniformMassExponentialTPicker : public massAndTPrimePicker {

	  public:

		uniformMassExponentialTPicker();
		virtual ~uniformMassExponentialTPicker() { }

		virtual bool init(const libconfig::Setting& setting);

		virtual bool operator() (double& invariantMass, double& tPrime);

		virtual std::ostream& print(std::ostream& out);

	  private:

		std::vector<std::pair<double, double> > _tSlopesForMassBins;

		std::pair<double, double> _tPrimeRange;

	};

}

#endif
