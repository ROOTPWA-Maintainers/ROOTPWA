#ifndef GENERATORWEIGHTFUNCTIONS_HH_
#define GENERATORWEIGHTFUNCTIONS_HH_


#include<limits>
#include<map>

#include<boost/shared_ptr.hpp>

#include<libconfig.h++>

#include<TF1.h>


namespace rpwa {

	class massAndTPrimePicker;
	typedef boost::shared_ptr<massAndTPrimePicker> massAndTPrimePickerPtr;

	class massAndTPrimePicker {

	  public:

		massAndTPrimePicker()
			: _massRange(std::pair<double, double>(0., 0.)),
			  _tPrimeRange(std::pair<double, double>(0., std::numeric_limits<double>::max())),
			  _initialized(false) { };

		massAndTPrimePicker(const massAndTPrimePicker& picker)
			: _massRange(picker._massRange),
			  _tPrimeRange(picker._tPrimeRange),
			  _initialized(picker._initialized) { };

		virtual ~massAndTPrimePicker() { };

		virtual bool init(const libconfig::Setting& setting) = 0;

		virtual void overrideMassRange(double lowerLimit, double upperLimit);
		virtual const std::pair<double, double>& massRange() const;

		virtual bool operator() (double& invariantMass, double& tPrime) = 0;

		virtual std::ostream& print(std::ostream& out) const = 0;

	  protected:

		virtual bool initTPrimeAndMassRanges(const libconfig::Setting& setting);

		std::pair<double, double> _massRange;
		std::pair<double, double> _tPrimeRange;
		bool _initialized;

	};

	inline std::ostream& operator<< (std::ostream& out, const massAndTPrimePicker& picker)
	{
		return picker.print(out);
	}

	class uniformMassExponentialTPicker : public massAndTPrimePicker {

	  public:

		uniformMassExponentialTPicker();
		uniformMassExponentialTPicker(const uniformMassExponentialTPicker& picker);
		virtual ~uniformMassExponentialTPicker() { }

		virtual bool init(const libconfig::Setting& setting);

		virtual bool operator() (double& invariantMass, double& tPrime);

		virtual std::ostream& print(std::ostream& out) const;

	  private:

		virtual std::ostream& printSlice(std::ostream& out, const std::vector<double>& param) const;

		std::map<double, std::vector<double> > _tSlopesForMassBins;
		unsigned int _nExponential;

	};

	// Class to pick a m and t' slope. First the mass is picked from a polynomial,
	// then the t' slope is determined in dependence of the mass (also a polynomial).
	class polynomialMassAndTPrimeSlopePicker : public massAndTPrimePicker {

	  public:

		polynomialMassAndTPrimeSlopePicker();
		polynomialMassAndTPrimeSlopePicker(const polynomialMassAndTPrimeSlopePicker& picker);
		virtual ~polynomialMassAndTPrimeSlopePicker() { }

		virtual bool init(const libconfig::Setting& setting);

		virtual bool operator() (double& invariantMass, double& tPrime);

		virtual std::ostream& print(std::ostream& out) const;

	  private:

		TF1 _massPolynomial;
		TF1 _tPrimeSlopePolynomial;

	};

}

#endif
