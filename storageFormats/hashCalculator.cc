
#include "hashCalculator.h"

#include <TVector3.h>

#include "reportingUtils.hpp"


using namespace std;
using namespace rpwa;


bool rpwa::hashCalculator::_debug = false;


void rpwa::hashCalculator::Update(const double& value) {
	if(_debug) {
		printDebug << "updating with " << value << "." << endl;
	}
	TMD5::Update((UChar_t*)&value, 8);
}


void rpwa::hashCalculator::Update(const complex<double>& value)
{
	Update(value.real());
	Update(value.imag());
}


void rpwa::hashCalculator::Update(const TVector3& vector) {
	if(_debug) {
		printDebug << "updating with TVector3(" << vector.X() << ", " << vector.Y() << ", " << vector.Z() << ")." << endl;
	}
	Update(vector.X());
	Update(vector.Y());
	Update(vector.Z());
}
