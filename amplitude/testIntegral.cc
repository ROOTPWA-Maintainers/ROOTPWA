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
// $Rev::                             $: revision of last commit
// $Author::                          $: author of last commit
// $Date::                            $: date of last commit
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


#include <algorithm>

#include <boost/assign/list_of.hpp>

#include "TROOT.h"
#include "TFile.h"

#include "reportingUtils.hpp"
#include "fileUtils.hpp"
#include "normalizationIntegral.h"


using namespace std;
using namespace boost::assign;
using namespace rpwa;


int
main(int argc, char** argv)
{
	// switch on debug output
	normalizationIntegral::setDebug(true);

	if (0) {
		const string somePath = "/local/data/compass/hadronData/massBins/2004/Q3PiData/r481.trunk/1260.1300/PSPAMPS/1-4++1+rho770_41_pi-.amp";
		cout << somePath << endl
		     << directoryFromPath(somePath) << endl
		     << fileNameFromPath(somePath) << endl
		     << fileNameNoExtFromPath(somePath) << endl
		     << extensionFromPath(somePath) << endl
		     << changeFileExtension(somePath, ".root") << endl
		     << changeFileExtension(somePath, "root") << endl
		     << changeFileExtension(somePath) << endl;
	}


	if (0) {
		// test integral calculation
		// get file list
		// vector<string> binAmpFileNames = filesMatchingGlobPattern("/data/compass/hadronData/massBins/2004/Q3PiData/r481.trunk/1260.1300/PSPAMPS/SYM/*.amp", true);
		// sort(binAmpFileNames.begin(), binAmpFileNames.end());
		// make sure order of waves is the same as used by int
		vector<string> binAmpFileNames = list_of
			("1-0-+0+f0980_00_pi-.amp"   )
			("1-0-+0+rho770_11_pi-.amp"  )
			("1-0-+0+sigma_00_pi-.amp"   )
			("1-1++0+f21270_12_pi-.amp"  )
			("1-1++0+rho770_01_pi-.amp"  )
			("1-1-+0-rho770_11_pi-.amp"  )
			("1-1++0+rho770_21_pi-.amp"  )
			("1-1++0+sigma_10_pi-.amp"   )
			("1-1++1+f21270_12_pi-.amp"  )
			("1-1++1-rho770_01_pi-.amp"  )
			("1-1++1+rho770_01_pi-.amp"  )
			("1-1-+1-rho770_11_pi-.amp"  )
			("1-1-+1+rho770_11_pi-.amp"  )
			("1-1++1+rho770_21_pi-.amp"  )
			("1-1++1+sigma_10_pi-.amp"   )
			("1-2-+0+f21270_02_pi-.amp"  )
			("1-2++0-f21270_12_pi-.amp"  )
			("1-2-+0+f21270_22_pi-.amp"  )
			("1-2-+0+rho770_11_pi-.amp"  )
			("1-2++0-rho770_21_pi-.amp"  )
			("1-2-+0+rho770_31_pi-.amp"  )
			("1-2-+0+sigma_20_pi-.amp"   )
			("1-2-+1-f21270_02_pi-.amp"  )
			("1-2-+1+f21270_02_pi-.amp"  )
			("1-2++1-f21270_12_pi-.amp"  )
			("1-2++1+f21270_12_pi-.amp"  )
			("1-2-+1+f21270_22_pi-.amp"  )
			("1-2-+1+rho770_11_pi-.amp"  )
			("1-2++1+rho770_21_pi-.amp"  )
			("1-2-+1+rho770_31_pi-.amp"  )
			("1-2-+1+sigma_20_pi-.amp"   )
			("1-3++0+f21270_12_pi-.amp"  )
			("1-3++0+rho31690_03_pi-.amp")
			("1-3++0+rho770_21_pi-.amp"  )
			("1-3++1+f21270_12_pi-.amp"  )
			("1-3++1+rho31690_03_pi-.amp")
			("1-3++1+rho770_21_pi-.amp"  )
			("1-4-+0+rho770_31_pi-.amp"  )
			("1-4++1+f21270_32_pi-.amp"  )
			("1-4-+1+rho770_31_pi-.amp"  )
			("1-4++1+rho770_41_pi-.amp"  );
		for (size_t i = 0; i < binAmpFileNames.size(); ++i)
			binAmpFileNames[i] = "/data/compass/hadronData/massBins/2004/Q3PiData/r481.trunk/1260.1300/PSPAMPS/SYM/" + binAmpFileNames[i];
		vector<string> rootAmpFileNames;
		normalizationIntegral integral;
		integral.integrate(binAmpFileNames, rootAmpFileNames);
		integral.writeAscii("testIntegral2.int");
	}


	if (1) {
		// test I/O and copying
		normalizationIntegral integral;
		// ascii I/O
		integral.readAscii("testIntegral.int");
		normalizationIntegral integral2(integral);
		integral2.writeAscii("testIntegral2.int");
		// root I/O
		// force loading predefined std::complex dictionary
#ifdef USE_STD_COMPLEX_TREE_LEAFS
		gROOT->ProcessLine("#include <complex>");
		{
			TFile* outFile = TFile::Open("testIntegral.root", "RECREATE");
			integral.Write("integral");
			outFile->Close();
		}
		{
			TFile*                 inFile    = TFile::Open("testIntegral.root", "READ");
			normalizationIntegral* integral3 = 0;
			inFile->GetObject("integral", integral3);
			if (not integral3)
				printErr << "cannot find integral 'integral'" << endl;
			else
				integral3->writeAscii("testIntegral3.int");
			inFile->Close();
		}
#endif  // USE_STD_COMPLEX_TREE_LEAFS
	}


	if (0) {
		// test renormalization
		normalizationIntegral integral;
		integral.readAscii("testIntegral.int");
		const unsigned int nmbEvents = integral.nmbEvents();
		integral.renormalize(10000);
		integral.renormalize(nmbEvents);  // scale back to original
		integral.writeAscii("testIntegral2.int");
		cout << integral;
		//integral.print(cout, true);
	}


}
