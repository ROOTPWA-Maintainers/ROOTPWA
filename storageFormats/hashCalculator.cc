#include "TVector3.h"

#include "hashCalculator.h"
#include "reportingUtils.hpp"
#include "reportingUtilsRoot.hpp"


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
	if(_debug) {
		printDebug << "updating with complex value " << value << "." << endl;
	}
	Update(value.real());
	Update(value.imag());
}


void rpwa::hashCalculator::Update(const TVector3& vector) {
	if(_debug) {
		printDebug << "updating with TVector3" << vector << "." << endl;
	}
	Update(vector.X());
	Update(vector.Y());
	Update(vector.Z());
}


void rpwa::hashCalculator::Update(const string& value) {
	if(_debug) {
		printDebug << "updating with string '" << value << "'." << endl;
	}
	TMD5::Update(reinterpret_cast<const UChar_t*>(value.data()), value.size());
}
