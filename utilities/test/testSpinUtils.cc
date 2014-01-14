///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
//
//    This file is part of rootpwa
//
//    rootpwa is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    rootpwa is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with rootpwa. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
//
// Description:
//      test some functions
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include "reportingUtils.hpp"
#include "conversionUtils.hpp"
#include "spinUtils.hpp"


using namespace std;
using namespace rpwa;
using namespace boost::tuples;


int
main(int    argc,
     char** argv)
{
	// test spin range
	tuple<int, int> demandedRange = make_tuple(0, 10);
	printInfo << "demanded range: [" << spinQn(get<0>(demandedRange)) << ", "
	          << spinQn(get<1>(demandedRange)) << "]" << endl;
	for (int s1 = 0; s1 < 5; ++s1)
		for (int s2 = 0; s2 < 5; ++s2) {
			bool valid;
			const tuple<int, int> allowedRange = getSpinRange(s1, s2, demandedRange, &valid);
			printInfo << "s1 = " << spinQn(s1) << ", s2 = " << spinQn(s2) << ": "
			          << "[" << spinQn(get<0>(allowedRange)) << ", " << spinQn(get<1>(allowedRange)) << "] "
			          << "valid = " << trueFalse(valid) << endl;
		}
	cout << endl;

	demandedRange = make_tuple(1, 2);
	printInfo << "demanded range: [" << spinQn(get<0>(demandedRange)) << ", "
	          << spinQn(get<1>(demandedRange)) << "]" << endl;
	for (int s1 = 0; s1 < 5; ++s1)
		for (int s2 = 0; s2 < 5; ++s2) {
			bool valid;
			const tuple<int, int> allowedRange = getSpinRange(s1, s2, demandedRange, &valid);
			printInfo << "s1 = " << spinQn(s1) << ", s2 = " << spinQn(s2) << ": "
			          << "[" << spinQn(get<0>(allowedRange)) << ", " << spinQn(get<1>(allowedRange)) << "] "
			          << "valid = " << trueFalse(valid) << endl;
		}
}
