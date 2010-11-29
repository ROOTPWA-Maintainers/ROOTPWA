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
// File and Version Information:
// $Rev:: 471                         $: revision of last commit
// $Author:: bgrube                   $: author of last commit
// $Date:: 2010-11-11 19:33:24 +0100 #$: date of last commit
//
// Description:
//      basic test program for integral
//
//
// Author List:
//      Boris Grube    TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include "normalizationIntegral.h"


using namespace std;
using namespace rpwa;


int
main(int argc, char** argv)
{
	// switch on debug output
	normalizationIntegral::setDebug(true);

	if (1) {
		normalizationIntegral integral;
		
		integral.readAscii("testIntegral.int");
		normalizationIntegral integral2(integral);
		integral2.writeAscii("testIntegral2.int");
	}
}
