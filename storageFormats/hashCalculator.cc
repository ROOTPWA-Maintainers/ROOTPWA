
#include "hashCalculator.h"

#include <TVector3.h>


using namespace std;
using namespace rpwa;


void rpwa::hashCalculator::Update(const double& value) {
	TMD5::Update((UChar_t*)&value, 8);
}

void rpwa::hashCalculator::Update(const TVector3& vector) {
	Update(vector.X());
	Update(vector.Y());
	Update(vector.Z());
}
