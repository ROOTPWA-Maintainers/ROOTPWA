#ifndef GENERATORWEIGHTFUNCTIONS_HH_
#define GENERATORWEIGHTFUNCTIONS_HH_


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

		uniformMassExponentialTPicker();
		virtual ~uniformMassExponentialTPicker() { }

		virtual bool init(const libconfig::Setting& setting);

		virtual bool operator() (double& invariantMass, double& tPrime);

		virtual std::ostream& print(std::ostream& out);

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
