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
//      test file functions
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include "reportingUtils.hpp"
#include "fileUtils.hpp"


using namespace std;
using namespace rpwa;


int
main(int    argc,
     char** argv)
{

	if (1) {
		const string somePath = "/local/data/compass/hadronData/massBins/2004/Q3PiData/r481.trunk/1260.1300/PSPAMPS/1-4++1+rho770_41_pi-.amp";
		printInfo << somePath << endl
		          << directoryFromPath(somePath) << endl
		          << fileNameFromPath(somePath) << endl
		          << fileNameNoExtFromPath(somePath) << endl
		          << extensionFromPath(somePath) << endl
		          << changeFileExtension(somePath, ".root") << endl
		          << changeFileExtension(somePath, "root") << endl
		          << changeFileExtension(somePath) << endl;
	}

}
