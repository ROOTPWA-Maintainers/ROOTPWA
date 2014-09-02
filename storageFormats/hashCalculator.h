
#ifndef HASHCALCULATOR_H
#define HASHCALCULATOR_H

#include <TMD5.h>

class TVector3;


namespace rpwa {

	class hashCalculator : public TMD5 {

	  public:

		hashCalculator()
			: TMD5() { }

		void Update(const double& value);
		void Update(const TVector3& vector);

		std::string hash() {
			TMD5::Final();
			return TMD5::AsString();
		}

	}; // class hashCalculator

} // namespace rpwa

#endif
