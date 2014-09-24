
#ifndef HASHCALCULATOR_H
#define HASHCALCULATOR_H

#include <complex>

#include <TMD5.h>

class TVector3;


namespace rpwa {

	class hashCalculator : public TMD5 {

	  public:

		hashCalculator()
			: TMD5() { }

		void Update(const double& value);
		void Update(const std::complex<double>& value);
		void Update(const TVector3& vector);

		std::string hash() {
			TMD5::Final();
			return TMD5::AsString();
		}

		const bool& debug() { return _debug; }
		void setDebug(const bool& debug = true) { _debug = debug; }

	  private:

		static bool _debug;

	}; // class hashCalculator

} // namespace rpwa

#endif
